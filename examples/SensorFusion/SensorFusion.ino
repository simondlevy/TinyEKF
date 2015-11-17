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

// These must be defined before including TinyEKF.h
#define N 1
#define M 1

#include <TinyEKF.h>

#include <SFE_BMP180.h>
#include <Wire.h>

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

        void model(double fx[N], double F[N][N], double hx[N], double H[M][N])
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
SFE_BMP180 baro;

void setup() {

    Serial.begin(9600);

    // Start reading from BMP180
    baro.begin();
}

void loop() {

    double temperature, pressure;
    getReadings(temperature, pressure);
    //Serial.print(temperature);
    //Serial.print(" ");
    //Serial.println(pressure);
}

void getReadings(double & T, double & P)
{
  char status;
  
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

      Serial.println(T);

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

        if (status == 0)
            Serial.println("error retrieving pressure measurement\n");
      }
      else Serial.println("error starting pressure measurement\n");
    }
    else Serial.println("error retrieving temperature measurement\n");
  }
  else Serial.println("error starting temperature measurement\n");
}
