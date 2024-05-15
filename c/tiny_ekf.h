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

        double x[N];     // state vector

        double P[N][N];  // prediction error covariance
        double Q[N][N];  // process noise covariance 
        double R[M][M];  // measurement error covariance

        double G[N][M];  // Kalman gain; a.k.a. K

        double F[N][N];  // Jacobian of process model
        double H[M][N];  // Jacobian of measurement model

        double Ht[N][M]; // transpose of measurement Jacobian
        double Ft[N][N]; // transpose of process Jacobian
        double Pp[N][N]; // P, post-prediction, pre-update

        double fx[N];   // output of user defined f() state-transition function
        double hx[M];   // output of user defined h() measurement function

      &nbsp; // temporary storage
        double tmp0[N][N];
        double tmp1[N][Msta];
        double tmp2[M][N];
        double tmp3[M][M];
        double tmp4[M][M];
        double tmp5[M]; 
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
int ekf_step(void * ekf, double * z);
