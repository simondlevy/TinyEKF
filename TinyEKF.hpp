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

#include <stdio.h>
#include <stdlib.h>

//extern "C" {
    void ekf_init(void *, int, int);
    int ekf_step(void *, double *);
//}

static void dump(double * a, int m, int n, const char * fmt)
{
    int i,j;

    char f[100];
    sprintf(f, "%s ", fmt);
    for(i=0; i<m; ++i) {
        for(j=0; j<n; ++j)
            printf(f, a[i*n+j]);
        printf("\n");
    }
}


/**
 * A header-only class for the Extended Kalman Filter.  Your implementing class should #define the constant _N and 
 * _M and then #include <TinyKalman.hpp>
 */
class TinyEKF {

    private:

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

        ekf_t ekf;

    protected:

        /**
          * The current state.
          */
        double * x;

        /**
         * Initializes a TinyEKF object.
         */
        TinyEKF() { 
            ekf_init(&this->ekf, _N, _M); 
            this->x = this->ekf.x; 
            //printf("TinyE:  x=%p\n", this->x);
            //exit(0);
        }

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

        /**
         * Sets the specified value of the prediction error covariance. <i>P<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setP(int i, int j, double value) 
        { 
            this->ekf.P[i][j] = value; 
        }

        /**
         * Sets the specified value of the process noise covariance. <i>Q<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
         void setQ(int i, int j, double value) 
        { 
            this->ekf.Q[i][j] = value; 
        }

        /**
         * Sets the specified value of the observation noise covariance. <i>R<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
         void setR(int i, int j, double value) 
        { 
            this->ekf.R[i][j] = value; 
        }

    public:

        /**
          * Returns the state element at a given index.
          * @param i the index (at least 0 and less than <i>n</i>
          * @return state value at index
          */
        double getX(int i) 
        { 
            return this->ekf.x[0]; 
        }

        double setX(int i, double value) 
        { 
            this->ekf.x[i] = value; 
        }

        /**
          Performs one step of the prediction and update.
          * @param z observation vector, length <i>m</i>
          * @return true on success, false on failure caused by non-positive-definite matrix.
         */
        bool step(double * z) 
        { 
            this->model(this->ekf.fx, this->ekf.F, this->ekf.hx, this->ekf.H); 

            return ekf_step(&this->ekf, z) ? false : true;
        }
};
