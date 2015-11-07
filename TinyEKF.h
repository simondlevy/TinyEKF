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

/**
  * A class for the Extended Kalman Filter.
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

        /**
        * State vector, length <i>n</i>.
        */
        double * x;

        /**
          * Constructs a TinyEKF object.
          * @param n number of states
          * @param m number of observable
          */
        TinyEKF(int n, int m);

        /**
          * Deallocates memory for a TinyEKF object.
          */
         ~TinyEKF();

        virtual void model(double * fx, double * F, double * hx, double * H) = 0;

        /**
          * A convience function for setting values in a matrix: <i>A<sub>i,j</sub> = value</i>
          @param A the matrix
          @param i row index (first = 0)
          @param j row index (first = 0)
          @param value value to set
          */
        void set(double * A, int i, int j, double value);

        /**
          * Sets a the value in the state noise covariance matrix <i>P</i>.
          @param i row index (first = 0)
          @param j row index (first = 0)
          @param value value to set
          */
         void setP(int i, int j, double value);

        /**
          * Sets a the value in the process noise covariance matrix <i>Q</i>.
          @param i row index (first = 0)
          @param j row index (first = 0)
          @param value value to set
          */
         void setQ(int i, int j, double value);

        /**
          * Sets a the value in the measurement noise covariance matrix <i>R</i>.
          @param i row index (first = 0)
          @param j row index (first = 0)
          @param value value to set
          */
         void setR(int i, int j, double value);

    public:

        /**
          * Returns the state element at a given index.
          * @param i the index (at least 0 and less than <i>n</i>
          * @return state value at index
          */
        double getX(int i);

        /**
          Performs one step of the prediction and update.
          @param z observation vector, length <i>m</i>
          @return true on success, false on failure caused by non-positive-definite matrix.
         */
        bool step(double * z);
};
