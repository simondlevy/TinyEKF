/* SensorFusion: Sensor fusion on Arduino using TinyEKF.  
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
#define N 1     // Two state values: pressure, temperature
#define M 2     // Three measurements: baro pressure, baro temperature, LM35 temperature

#define LM35_PIN 0

#include <TinyEKF.h>
#include <SFE_BMP180.h>
#include <Wire.h>

class Fuser : public TinyEKF {

    public:

        Fuser()
        {            
            this->setQ(0, 0, 1);

            this->setR(0, 0, 1);
            this->setR(1, 1, 1);
        }

    protected:

        void model(double fx[N], double F[N][N], double hx[M], double H[M][N])
        {
            fx[0] = this->x[0];

            F[0][0] = 1;

            hx[0] = this->x[0];
            hx[1] = this->x[0];

            H[0][0] = 1;
            H[1][0] = 1;
        }
};


Fuser ekf;
SFE_BMP180 baro;


void setup() {

    Serial.begin(9600);

    // Start reading from baro
    baro.begin();

    // Set up to read from LM35
    analogReference(INTERNAL);
}

void loop() {

    float lm35Temperature = analogRead(LM35_PIN) / 9.31;

    Serial.println(lm35Temperature); return;

    double baroTemperature, baroPressure;
    getBaroReadings(baroTemperature, baroPressure);

    double z[2] = {baroTemperature, lm35Temperature};

    ekf.step(z);

    Serial.println(ekf.getX(0));
}

void getBaroReadings(double & T, double & P)
{
    char status = baro.startTemperature();
   
    if (status != 0) {
        delay(status);
        status = baro.getTemperature(T);
        if (status != 0) {
            status = baro.startPressure(3);
            if (status != 0) {
                delay(status);
                status = baro.getPressure(P,T);
                if (status == 0)
                    Serial.println("error retrieving pressure measurement");
            }
            else Serial.println("error starting pressure measurement");
        }
        else Serial.println("error retrieving temperature measurement");
    }
    else Serial.println("error starting temperature measurement");
}
