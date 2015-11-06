/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
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

class TinyEKF {

    private:

        int n;          // number of states
        int m;          // number of measurements

        double * P;     // prediction error covariance
        double * Q;     // process noise covariance
        double * R;     // measurement error covariance

        double * G;     // Kalman gain; a.k.a. K
        double * F;     // Jacobian of process model
        double * H;     // Jacobian of measurement model

        double * Ht;    // transpose of measurement Jacobian
        double * Ft;    // transpose of process Jacobian
        double * Pp;    // P, post-prediction, pre-update

        double * fx;    // output of user defined f() state-transition function
        double * hx;    // output of user defined h() measurement function

        // temporary storage
        double * tmp1;
        double * tmp2;
        double * tmp3;
        double * tmp4;
        double * tmp5;

    protected:

        double * x;     // state

        TinyEKF(int n, int m);

        ~TinyEKF();

        virtual void model(double * fx, double * F, double * hx, double * H) = 0;

        void set(double * A, int i, int j, double value);

        void setP(int i, int j, double value);

        void setQ(int i, int j, double value);

        void setR(int i, int j, double value);

    public:

        double getX(int i);

        /**
          Returns true on success, false on failure caused by non-positive-definite matrix.
         */
        bool step(double * z);
};

void error(const char * msg);
