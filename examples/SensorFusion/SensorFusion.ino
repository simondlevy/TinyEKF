/* SensorFusion: Sensor fusion on Arduino using TinyEKF
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

class Fuser : public TinyEKF {

    public:

    Fuser()
    {            
        for (int i=0; i<2; ++i)
            for (int j=0; j<2; ++j) {
                this->setQ(i, j, .01);
                this->setR(i, j, .01);
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

void setup() {

    Serial.begin(9600);

    ekf.setX(0, 0);
    ekf.setX(1, 0);
}

void loop() {

    static int count;
    const int LOOPSIZE = 1000;

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
}
