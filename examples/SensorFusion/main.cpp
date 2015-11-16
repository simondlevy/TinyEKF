/*
 * Command-line test program for TinyEKF SensorFusion demo.
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

    const int PERIOD = 10;

    int i = 0;

    while (true) {

        double z[2];

        double t = M_PI*i/(PERIOD-1);

        z[0] = sin(t);
        z[1] = cos(t);

        ekf.step(z);

        printf("%f %f\n", ekf.getX(0), ekf.getX(1));

        i = (i+1) % PERIOD;
    }

    return 0;
}
