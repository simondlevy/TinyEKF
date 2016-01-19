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
#define N 2     // Two state values: pressure, temperature
#define M 3     // Three measurements: baro pressure, baro temperature, LM35 temperature

#define LM35_PIN 0

#include <TinyEKF.h>
#include <SFE_BMP180.h>
#include <Wire.h>

class Fuser : public TinyEKF {

    public:

        Fuser()
        {            
            // We approximate the process noise using a small constant
            this->setQ(0, 0, .0001);
            this->setQ(1, 1, .0001);

            // Same for measurement noise
            this->setR(0, 0, .0001);
            this->setR(1, 1, .0001);
            this->setR(2, 2, .0001);
        }

    protected:

        void model(double fx[N], double F[N][N], double hx[M], double H[M][N])
        {
            // Process model is f(x) = x
            fx[0] = this->x[0];
            fx[1] = this->x[1];

            // So process model Jacobian is identity matrix
            F[0][0] = 1;
            F[1][1] = 1;

            // Measurement function
            hx[0] = this->x[0]; // Barometric pressure from previous state
            hx[1] = this->x[1]; // Baro temperature from previous state
            hx[2] = this->x[1]; // LM35 temperature from previous state

            // Jacobian of measurement function
            H[0][0] = 1;        // Barometric pressure from previous state
            H[1][1] = 1 ;       // Baro temperature from previous state
            H[2][1] = 1 ;       // LM35 temperature from previous state
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

    // Read pressure, temperature from BMP180
    double baroTemperature, baroPressure;
    getBaroReadings(baroTemperature, baroPressure);

    // Read temperature from LM35
    float lm35Temperature = analogRead(LM35_PIN) / 9.31;

    // Send these measurements to the EKF
    double z[3] = {baroPressure, baroTemperature, lm35Temperature};
    ekf.step(z);

    // Report measured and predicte/fused values
    Serial.print(z[0]);
    Serial.print(" ");
    Serial.print(z[1]);
    Serial.print(" ");
    Serial.print(z[2]);
    Serial.print(" ");
    Serial.print(ekf.getX(0));
    Serial.print(" ");
    Serial.println(ekf.getX(1));
}


// Adapted from https://github.com/sparkfun/BMP180_Breakout
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
