/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
 *
 * tinyekf_config.h: static configuration parameters
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include <stdbool.h>

typedef struct {

    double x[EKF_N];     // state vector

    double P[EKF_N*EKF_N];  // prediction error covariance
    double Q[EKF_N*EKF_N];  // process noise covariance
    double R[EKF_M*EKF_M];  // measurement error covariance

    double G[EKF_N*EKF_M];  // Kalman gain; a.k.a. K

    double F[EKF_N*EKF_N];  // Jacobian of process model
    double H[EKF_M*EKF_N];  // Jacobian of measurement model

    double fx[EKF_N];   // output of user defined f() state-transition function
    double hx[EKF_M];   // output of user defined h() measurement function

    double Pp[EKF_N*EKF_N]; // P, post-prediction, pre-update

} ekf_t;

void ekf_init(ekf_t * ekf);

bool ekf_step(ekf_t * ekf, double * z);
