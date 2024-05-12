/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include <stdio.h>
#include <stdlib.h>

typedef struct {

    int n;          /* number of state values */
    int m;          /* number of observables */

    float x[EKF_N];    /* state vector */

    float P[EKF_N][EKF_N];  /* prediction error covariance */
    float Q[EKF_N][EKF_N];  /* process noise covariance */
    float R[EKF_M][EKF_M];  /* measurement error covariance */

    float G[EKF_N][EKF_M];  /* Kalman gain; a.k.a. K */

    float F[EKF_N][EKF_N];  /* Jacobian of process model */
    float H[EKF_M][EKF_N];  /* Jacobian of measurement model */

    float Ht[EKF_N][EKF_M]; /* transpose of measurement Jacobian */
    float Ft[EKF_N][EKF_N]; /* transpose of process Jacobian */
    float Pp[EKF_N][EKF_N]; /* P, post-prediction, pre-update */

    float fx[EKF_N];   /* output of user defined f() state-transition function */
    float hx[EKF_M];   /* output of user defined h() measurement function */

    /* temporary storage */
    float tmp0[EKF_N][EKF_N];
    float tmp1[EKF_N][EKF_M];
    float tmp2[EKF_M][EKF_N];
    float tmp3[EKF_M][EKF_M];
    float tmp4[EKF_M][EKF_M];
    float tmp5[EKF_M]; 

} ekf_t;        

// Support both Arduino and command-line versions
#ifndef MAIN
extern "C" {
#endif
    void ekf_init(void *, int, int);
    int ekf_step(void *, float *);
#ifndef MAIN
}
#endif

/**
 * A header-only class for the Extended Kalman Filter.  Your implementing class should #define the constant N and 
 * and then #include <TinyEKF.h>  You will also need to implement a model() method for your application.
 */
class TinyEKF {

    private:

        ekf_t ekf;

    protected:

        /**
          * The current state.
          */
        float * x;

        /**
         * Initializes a TinyEKF object.
         */
        TinyEKF() { 
            ekf_init(&this->ekf, EKF_N, EKF_M); 
            this->x = this->ekf.x; 
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
        virtual void model(float fx[EKF_N], float F[EKF_N][EKF_N], float hx[EKF_M], float H[EKF_M][EKF_N]) = 0;

        /**
         * Sets the specified value of the prediction error covariance. <i>P<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setP(int i, int j, float value) 
        { 
            this->ekf.P[i][j] = value; 
        }

        /**
         * Sets the specified value of the process noise covariance. <i>Q<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setQ(int i, int j, float value) 
        { 
            this->ekf.Q[i][j] = value; 
        }

        /**
         * Sets the specified value of the observation noise covariance. <i>R<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setR(int i, int j, float value) 
        { 
            this->ekf.R[i][j] = value; 
        }

    public:

        /**
         * Returns the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @return state value at index
         */
        float getX(int i) 
        { 
            return this->ekf.x[i]; 
        }

        /**
         * Sets the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @param value value to set
         */
        void setX(int i, float value) 
        { 
            this->ekf.x[i] = value; 
        }

        /**
          Performs one step of the prediction and update.
         * @param z observation vector, length <i>m</i>
         * @return true on success, false on failure caused by non-positive-definite matrix.
         */
        bool step(float * z) 
        { 
            this->model(this->ekf.fx, this->ekf.F, this->ekf.hx, this->ekf.H); 

            return ekf_step(&this->ekf, z) ? false : true;
        }
};
