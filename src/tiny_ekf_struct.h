/*
 * tiny_ekf_struct.h: common data structure for TinyEKF
 *
 * You should #include this file after using #define for N (states) and M
*  (observations)
 *
 * Copyright (C) 2016 Simon D. Levy
 *
 * MIT License
 */

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
