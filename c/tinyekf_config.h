/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
 *
 * tinyekf_config.h: static configuration parameters
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

typedef struct {

    int n;           /* number of state values */
    int m;           /* number of observables */

    double x[EKF_N];     /* state vector */

    double P[EKF_N][EKF_N];  /* prediction error covariance */
    double Q[EKF_N][EKF_N];  /* process noise covariance */
    double R[EKF_M][EKF_M];  /* measurement error covariance */

    double G[EKF_N][EKF_M];  /* Kalman gain; a.k.a. K */

    double F[EKF_N][EKF_N];  /* Jacobian of process model */
    double H[EKF_M][EKF_N];  /* Jacobian of measurement model */

    double Ht[EKF_N][EKF_M]; /* transpose of measurement Jacobian */
    double Ft[EKF_N][EKF_N]; /* transpose of process Jacobian */
    double Pp[EKF_N][EKF_N]; /* P, post-prediction, pre-update */

    double fx[EKF_N];   /* output of user defined f() state-transition function */
    double hx[EKF_M];   /* output of user defined h() measurement function */

    /* temporary storage */
    double tmp0[EKF_N][EKF_N];
    double tmp1[EKF_N][EKF_M];
    double tmp2[EKF_M][EKF_N];
    double tmp3[EKF_M][EKF_M];
    double tmp4[EKF_M][EKF_M];
    double tmp5[EKF_M]; 

} ekf_t;
