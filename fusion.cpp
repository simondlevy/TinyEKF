/*
References:
1. R G Brown, P Y C Hwang, "Introduction to random signals and applied 
Kalman filtering : with MATLAB exercises and solutions",1996
2. Pratap Misra, Per Enge, "Global Positioning System Signals, 
Measurements, and Performance(Second Edition)",2006
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

#include "tinyekf.hpp"

class Fuser : public TinyEKF {

    public:

        // Eight state values, four measurement values
        Fuser() : TinyEKF(1, 1)
        {            
            this->setP(0, 0, .01);
            this->setQ(0, 0, .01);
            this->setR(0, 0, .01);

        }

    protected:

        void model(double * fx, double * F, double * hx, double * H)
        {
        }
};

int main(int argc, char ** argv)
{    
    // Create the EKF
    Fuser ekf;

    static const int steps = 1000;

    // Loop till no more data
    for (int i=0; i<steps; ++i) {

        double z[1] = {0};

        ekf.step(z);

        printf("%f\n", ekf.getX(0));
    }
}
