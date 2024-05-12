/* SensorFusion: Sensor fusion on Arduino using TinyEKF.  
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */


// These must be defined before including TinyEKF.h
#define EKF_N 2     // pressure, temperature
#define EKF_M 3     // baro pressure, baro temperature, LM35 temperature

static const uint8_t LM35_PIN = 0;

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

        void model(
                double fx[EKF_N], 
                double F[EKF_N][EKF_N], 
                double hx[EKF_M], 
                double H[EKF_M][EKF_N])
        {
            // Process model is f(x) = x
            fx[0] = this->x[0];
            fx[1] = this->x[1];

            // So process model Jacobian is identity matrix
            F[0][0] = 1;
            F[1][1] = 1;

            // Measurement function simplifies the relationship between state
            // and sensor readings for convenience.  A more realistic
            // measurement function would distinguish between state value and
            // measured value; e.g.:
            //   hx[0] = pow(this->x[0], 1.03);
            //   hx[1] = 1.005 * this->x[1];
            //   hx[2] = .9987 * this->x[1] + .001;
            hx[0] = this->x[0]; // Barometric pressure from previous state
            hx[1] = this->x[1]; // Baro temperature from previous state
            hx[2] = this->x[1]; // LM35 temperature from previous state

            // Jacobian of measurement function
            H[0][0] = 1;        // Barometric pressure from previous state
            H[1][1] = 1 ;       // Baro temperature from previous state
            H[2][1] = 1 ;       // LM35 temperature from previous state
        }
};

static Fuser ekf;

static SFE_BMP180 baro;

// Adapted from https://github.com/sparkfun/BMP180_Breakout
static void getBaroReadings(double & T, double & P)
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

void setup() {

    Serial.begin(115200);

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
    Serial.print("BMP180Press:");
    Serial.print(z[0]);
    Serial.print(" ");
    Serial.print(" BMP180Temp:");
    Serial.print(z[1]);
    Serial.print(" LM35Temp:");
    Serial.print(z[2]);
    Serial.print(" EKFPress:");
    Serial.print(ekf.getX(0));
    Serial.print(" EKFTemp:");
    Serial.println(ekf.getX(1));
}



