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

#ifndef _float_t
typedef float _float_t;
#endif

typedef struct {

    _float_t x[EKF_N];     // state vector

    _float_t P[EKF_N*EKF_N];  // prediction error covariance
    _float_t Q[EKF_N*EKF_N];  // process noise covariance
    _float_t R[EKF_M*EKF_M];  // measurement error covariance

    _float_t G[EKF_N*EKF_M];  // Kalman gain; a.k.a. K

    _float_t F[EKF_N*EKF_N];  // Jacobian of process model
    _float_t H[EKF_M*EKF_N];  // Jacobian of measurement model

    _float_t fx[EKF_N];   // output of user defined f() state-transition function
    _float_t hx[EKF_M];   // output of user defined h() measurement function

    _float_t Pp[EKF_N*EKF_N]; // P, post-prediction, pre-update

} ekf_t;

void ekf_initialize(ekf_t * ekf);

void ekf_predict(ekf_t * ekf);

// Returns false iff matrix inversion fails
bool ekf_update(ekf_t * ekf, _float_t * z);
