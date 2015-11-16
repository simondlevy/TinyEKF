#define _N 2
#define _M 2

#include "TinyEKF.hpp"
#include <math.h>
#include <stdio.h>

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

int main(int argc, char ** argv)
{
    Fuser ekf;

    ekf.setX(0, 0);
    ekf.setX(1, 0);

    const int STEPS = 10;

    for (int i=0; i<STEPS; ++i) {

        double z[2];

        double t = M_PI*i/(STEPS-1);

        z[0] = sin(t);
        z[1] = cos(t);

        ekf.step(z);
    }

    return 0;
}
