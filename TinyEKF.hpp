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

#include <tinyekf.h>

class TinyEKF {

    protected:

        typedef struct {

            int n;          /* number of state values */
            int m;          /* number of observables */

            double x[_N];    /* state vector */

            double P[_N][_N];  /* prediction error covariance */
            double Q[_N][_N];  /* process noise covariance */
            double R[_M][_M];  /* measurement error covariance */

            double G[_N][_M];  /* Kalman gain; a.k.a. K */

            double F[_N][_N];  /* Jacobian of process model */
            double H[_M][_N];  /* Jacobian of measurement model */

            double Ht[_N][_M]; /* transpose of measurement Jacobian */
            double Ft[_N][_N]; /* transpose of process Jacobian */
            double Pp[_N][_N]; /* P, post-prediction, pre-update */

            double fx[_N];   /* output of user defined f() state-transition function */
            double hx[_N];   /* output of user defined h() measurement function */

            /* temporary storage */
            double tmp1[_N][_N];
            double tmp2[_M][_N];
            double tmp3[_M][_M];
            double tmp4[_M][_M];
            double tmp5[_M]; 

        } ekf_t;        


        /**
         * Initializes a TinyEKF object.
         */
        TinyEKF() { }

        /**
         * Deallocates memory for a TinyEKF object.
         */
        ~TinyEKF() { }

        /**
         * Implement this function for your EKF model.
         * @param fx gets output of state-transition function <i>f(x<sub>0 .. n-1</sub>)</i>
         * @param F gets <i>n &times; n</i> Jacobian of <i>f(x)</i>
         * @param hx gets output of observation function <i>h(x<sub>0 .. n-1</sub>)</i>
         * @param H gets <i>m &times; n</i> Jacobian of <i>h(x)</i>
         */
        virtual void model(double fx[_N], double F[_N][_N], double hx[_N], double H[_M][_N]) = 0;

    public:

        /**
          * Returns the state element at a given index.
          * @param i the index (at least 0 and less than <i>n</i>
          * @return state value at index
          */
        double getX(int i) { }

        /**
          Performs one step of the prediction and update.
          * @param z observation vector, length <i>m</i>
          * @return true on success, false on failure caused by non-positive-definite matrix.
         */
        bool step(double * z) { return true;}
};
