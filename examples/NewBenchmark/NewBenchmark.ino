/* Benchmark: Benchmark TinyEKF on Arduino/Teensy.  
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */


#define EKF_N 5     // states, will also be measurement values M
#define EKF_M EKF_N

#include <NewTinyEKF.h>

class Fuser : public TinyEKF {

    public:

        Fuser()
        {            
            // We approximate the process noise using a small constant
            for (int j=0; j<EKF_N; ++j)
                this->setQ(j, j, .0001);

            // Same for measurement noise
            for (int j=0; j<EKF_N; ++j)
                this->setR(j, j, .0001);
        }

    protected:

        void model(
                float fx[EKF_N], 
                float F[EKF_N][EKF_N], 
                float hx[EKF_M], 
                float H[EKF_M][EKF_N])
        {
            
            for (int j=0; j<EKF_N; ++j) {

                // Process model is f(x) = x
                fx[j] = x[j];

                // So process model Jacobian is identity matrix
                F[j][j] = 1;

                // EKF_Measurement function
                hx[j] = this->x[j]; 

                // Jacobian of measurement function
                H[j][j] = 1;  
            }
        }
};

static Fuser ekf;
static unsigned long count;
static float z[EKF_M];

void setup() {

    Serial.begin(115200);

    // Fake up some sensor readings for later
    for (int j=0; j<EKF_N; ++j)
        z[j] = j;

    count = 0;
}


void loop() {

    ekf.step(z);

    count++;

    const auto msec = millis();

    static uint32_t msec_prev;

    if (msec - msec_prev > 1000) {
        Serial.print("N = M = ");
        Serial.print(EKF_N);
        Serial.print(" : ");
        Serial.print(count);
        Serial.println(" updates per second");
        count = 0;
        msec_prev = msec;
    }
}
