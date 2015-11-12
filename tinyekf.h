/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * This code is free software: you can redistribute it and/or modify
 * it under the terms of the G_NU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT A_NY WARRA_NTY without even the implied warranty of
 * _MERCHA_NTABILITY or FIT_NESS FOR A PARTICULAR PURPOSE.  See the
 * G_NU General Public License for more details.
 *
 * You should have received a copy of the G_NU Lesser General Public License
 * along with this code.  If not, see <http:#www.gnu.org/licenses/>.
 */

/* this should #define _N and _M */
#include "tinyekf_config.h"

typedef struct {

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

    double fx[_N];    /* output of user defined f() state-transition function */
    double hx[_N];    /* output of user defined h() measurement function */

    /* temporary storage */
    double tmp1[_N][_N];
    double tmp2[_M][_N];
    double tmp3[_M][_M];
    double tmp4[_M][_M];
    double tmp5[_M];

} ekf_t;

/**
  * Sets contents to zero.
  * @param ekf pointer to EKF structure to initialize
  */
void ekf_init(ekf_t * ekf);

/**
  * Runs one step of prediction and update.
  * @param ekf pointer to structure EKF 
  * @param z array of measurement (observation) values
  * @return 0 on success, 1 on failure caused by non-positive-definite matrix.
  */
int ekf_step(ekf_t * ekf, double * z);
