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

#include "TinyEKF.h"

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
            fx[0] = this->x[0];
            hx[0] = fx[0];

            set(F, 0, 0, 1);
            set(H, 0, 0, 1);
        }
};

Fuser ekf;

void setup() {

    Serial.begin(9600);
}


void loop() {

    static int count;
    const int LOOPSIZE = 1000;

    double z[1];

    z[0] = sin(2*M_PI*count/LOOPSIZE);

    ekf.step(z);

    Serial.println(ekf.getX(0));

    count = (count + 1) % LOOPSIZE;
}
