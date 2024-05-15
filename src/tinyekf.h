/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
 *
 * tinyekf_config.h: static configuration parameters
 *
 * Copyright (C) 2024 Simon D. Levy
 *
 * MIT License
 */

#include <stdbool.h>

#ifndef _float_t
typedef float _float_t;
#endif

typedef struct {

    _float_t x[EKF_N];         // state vector
    _float_t P[EKF_N*EKF_N];  // prediction error covariance
    _float_t Pp[EKF_N*EKF_N]; // P, post-prediction, pre-update

} ekf_t;

void ekf_initialize(ekf_t * ekf);

void ekf_predict(
        ekf_t * ekf, 
        const _float_t fx[EKF_N],
        const _float_t F[EKF_N*EKF_N],
        const _float_t Q[EKF_N*EKF_N]);

// Returns false iff matrix inversion fails
bool ekf_update(
        ekf_t * ekf, 
        const _float_t z[EKF_M], 
        const _float_t hx[EKF_N],
        const _float_t H[EKF_M*EKF_N],
        const _float_t R[EKF_M*EKF_M]);