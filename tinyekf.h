/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
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

#include "tinyekf_config.h"

typedef struct {

    int n;          /* number of states */
    int m;          /* number of measurements */

    double x[N];    /* state vector */

    double P[N*N];  /* prediction error covariance */
    double Q[N*N];  /* process noise covariance */
    double R[M*M];  /* measurement error covariance */

    double G[N*M];  /* Kalman gain; a.k.a. K */
    double * F;     /* Jacobian of process model */
    double * H;     /* Jacobian of measurement model */

    double Ht[N*M]; /* transpose of measurement Jacobian */
    double Ft[N*N];    /* transpose of process Jacobian */
    double * Pp;    /* P, post-prediction, pre-update */

    double fx[N];    /* output of user defined f() state-transition function */
    double hx[N];    /* output of user defined h() measurement function */

    /* temporary storage */
    double * tmp1;
    double * tmp2;
    double * tmp3;
    double * tmp4;
    double * tmp5;

} ekf_t;

void   ekf_init(ekf_t * ekf, int n, int m);

void   ekf_free(ekf_t * ekf);

/**
  * Write this function for you application.
  */
void ekf_model(double * x, double * fx, double * F, double * hx, double * H);

void ekf_set(ekf_t * ekf, double * A, int i, int j, double value);

void ekf_setP(ekf_t * ekf, int i, int j, double value);

void ekf_setQ(ekf_t * ekf, int i, int j, double value);

void ekf_setR(ekf_t * ekf, int i, int j, double value);

double ekf_getX(ekf_t * ekf, int i);

/**
  * @return 0 on success, 1 on failure caused by non-positive-definite matrix.
  */
int ekf_step(ekf_t * ekf, double * z);
