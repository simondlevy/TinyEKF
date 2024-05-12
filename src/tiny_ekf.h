/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */


/**
  * Initializes an EKF structure.
  * @param ekf pointer to EKF structure to initialize
  * @param n number of state variables
  * @param m number of observables
  *
  * <tt>ekf</tt> should be a pointer to a structure defined as follows, where <tt>N</tt> and </tt>M</tt> are 
  * constants:
  * <pre>
        int n;           // number of state values 
        int m;           // number of observables 

        float x[N];     // state vector

        float P[N][N];  // prediction error covariance
        float Q[N][N];  // process noise covariance 
        float R[M][M];  // measurement error covariance

        float G[N][M];  // Kalman gain; a.k.a. K

        float F[N][N];  // Jacobian of process model
        float H[M][N];  // Jacobian of measurement model

        float Ht[N][M]; // transpose of measurement Jacobian
        float Ft[N][N]; // transpose of process Jacobian
        float Pp[N][N]; // P, post-prediction, pre-update

        float fx[N];   // output of user defined f() state-transition function
        float hx[M];   // output of user defined h() measurement function

      &nbsp; // temporary storage
        float tmp0[N][N];
        float tmp1[N][Msta];
        float tmp2[M][N];
        float tmp3[M][M];
        float tmp4[M][M];
        float tmp5[M]; 
    * </pre>
  */
void ekf_init(void * ekf, int n, int m);

/**
  * Runs one step of EKF prediction and update. Your code should first build a model, setting
  * the contents of <tt>ekf.fx</tt>, <tt>ekf.F</tt>, <tt>ekf.hx</tt>, and <tt>ekf.H</tt> to appropriate values.
  * @param ekf pointer to structure EKF 
  * @param z array of measurement (observation) values
  * @return 0 on success, 1 on failure caused by non-positive-definite matrix.
  */
int ekf_step(void * ekf, float * z);
