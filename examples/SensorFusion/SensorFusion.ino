/* SensorFusion: Sensor fusion on Arduino using TinyEKF.  Currently just uses sine and cosine
 *  waves as measurement proxy.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * This code is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this code.  If not, see <http:#www.gnu.org/licenses/>.
 */

// Un-comment this to continue working on GPS / Barometer fusion
//#define ACTUAL

#ifdef ACTUAL

#include <SFE_BMP180.h>
#include <UBX_Parser.h>
#include <Wire.h>

SFE_BMP180 baro;
double pressureBaseline; 

class POSLLH_Parser : public UBX_Parser {

    void handle_NAV_POSLLH(unsigned long iTOW, 
            long lon, 
            long lat, 
            long height, 
            long hMSL, 
            unsigned long hAcc, 
            unsigned long vAcc) {

        // Grab starting altitude first time around
        if (!this->ready) {
          this->hMSL_Basline = hMSL;
        }
        this->ready = true;

        this->lat = lat / 1e7;
        this->lon = lon / 1e7;
        this->alt = (hMSL-this->hMSL_Basline)/1e3;

    }  

  public:

   double lat;
   double lon;
   double alt;

 private:
   double hMSL_Basline;
   bool ready;
    
};


POSLLH_Parser parser;


bool success;

#else

#define _N 2
#define _M 2

#include <Arduino.h>
#include "TinyEKF.hpp"



class Fuser : public TinyEKF {

    public:

    Fuser()
    {            
        for (int i=0; i<2; ++i) {
            this->setQ(i, i, .01);
            this->setR(i, i, .01);
        }
     }

    protected:

        void model(double fx[_N], double F[_N][_N], double hx[_N], double H[_M][_N])
        {
            fx[0] = this->x[0];
            fx[1] = this->x[1];

            hx[0] = fx[0];
            hx[1] = fx[1];

            for (int i=0; i<2; ++i) {
                F[i][i] = 1;
                H[i][i] = 1;
            }
        }
};

Fuser ekf;
#endif

void setup() {

    Serial.begin(9600);

#ifdef ACTUAL
  // Start reading from BMP180
  baro.begin();

  // Start reading from UBlox
  Serial1.begin(57600);

  pressureBaseline = getPressure();
#else
    ekf.setX(0, 0);
    ekf.setX(1, 0);
#endif
}

void loop() {

#ifdef ACTUAL

   static double baroRelativeAltitude;

   // If byte available from GPS, parse it
    if (Serial1.available()) {
        parser.parse(Serial1.read());
    }

    // Otherwise, read from baro
    else {

        baroRelativeAltitude = baro.altitude(getPressure(), pressureBaseline);

        // tiny pause to avoid contention
        delay(100);
    }

    Serial.print(parser.lat, 7);
    Serial.print(" ");
    Serial.print(parser.lon, 7);
    Serial.print(" ");
    Serial.print(parser.alt);
    Serial.print(" ");
    Serial.print(baroRelativeAltitude);
    Serial.println();
#else

    static int count;
    const int LOOPSIZE = 100;

    double z[2];

    double t = 2*M_PI*count/(LOOPSIZE-1);

    z[0] = sin(t);
    z[1] = cos(t);

    ekf.step(z);

    Serial.print(ekf.getX(0));
    Serial.print(" ");
    Serial.print(ekf.getX(1));
    Serial.println();

    count = (count + 1) % LOOPSIZE;
    
#endif

    
}

#ifdef ACTUAL
double getPressure()
{
  char status;
  double T,P;
  
  // You must first get a temperature measurement to perform a pressure reading.
  
  // Start a temperature measurement:
  // If request is successful, the number of ms to wait is returned.
  // If request is unsuccessful, 0 is returned.

  status = baro.startTemperature();

  if (status != 0)
  {
    // Wait for the measurement to complete:
    delay(status);

    // Retrieve the completed temperature measurement:
    // Note that the measurement is stored in the variable T.
    // Use '&T' to provide the address of T to the function.
    // Function returns 1 if successful, 0 if failure.

    status = baro.getTemperature(T);
    if (status != 0)
    {
      // Start a pressure measurement:
      // The parameter is the oversampling setting, from 0 to 3 (highest res, longest wait).
      // If request is successful, the number of ms to wait is returned.
      // If request is unsuccessful, 0 is returned.

      status = baro.startPressure(3);
      if (status != 0)
      {
   
        // Wait for the measurement to complete:
        delay(status);

        // Retrieve the completed pressure measurement:
        // Note that the measurement is stored in the variable P.
        // Use '&P' to provide the address of P.
        // Note also that the function requires the previous temperature measurement (T).
        // (If temperature is stable, you can do one temperature measurement for a number of pressure measurements.)
        // Function returns 1 if successful, 0 if failure.

        status = baro.getPressure(P,T);
        if (status != 0)
        {
          return(P);
        }
        else Serial.println("error retrieving pressure measurement\n");
      }
      else Serial.println("error starting pressure measurement\n");
    }
    else Serial.println("error retrieving temperature measurement\n");
  }
  else Serial.println("error starting temperature measurement\n");

  return -1;
}
#endif
