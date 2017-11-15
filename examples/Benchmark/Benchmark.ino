/* Benchmark: Benchmark TinyEKF on Arduino/Teensy.  
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */


#define Nsta 5     // states, will also be measurement values M
#define Mobs Nsta

#include <TinyEKF.h>
#include <elapsedMillis.h>

//#define DEBUG   // comment this out to do timing

class Fuser : public TinyEKF {

    public:

        Fuser()
        {            
            // We approximate the process noise using a small constant
            for (int j=0; j<Nsta; ++j)
                this->setQ(j, j, .0001);

            // Same for measurement noise
            for (int j=0; j<Nsta; ++j)
                this->setR(j, j, .0001);
        }

    protected:

        void model(double fx[Nsta], double F[Nsta][Nsta], double hx[Mobs], double H[Mobs][Nsta])
        {
            
            for (int j=0; j<Nsta; ++j) {

                // Process model is f(x) = x
                fx[j] = x[j];

                // So process model Jacobian is identity matrix
                F[j][j] = 1;

                // Mobseasurement function
                hx[j] = this->x[j]; 

                // Jacobian of measurement function
                H[j][j] = 1;  
            }
        }
};

Fuser ekf;
unsigned long count;
elapsedMillis timeElapsed; 
double z[Mobs];

void setup() {

    Serial.begin(9600);

    // Fake up some sensor readings for later
    for (int j=0; j<Nsta; ++j)
        z[j] = j;

    count = 0;
}


void loop() {

    ekf.step(z);

#ifdef DEBUG
    for (int j=0; j<Nsta; ++j) {
        Serial.print(ekf.getX(j));
        Serial.print(" ");
    }
    Serial.println();
#else
    count++;

    if (timeElapsed > 1000) {
        Serial.print("N = M = ");
        Serial.print(Nsta);
        Serial.print(" : ");
        Serial.print(count);
        Serial.println(" updates per second");
        timeElapsed = 0;
        count = 0;
    }
#endif
}
