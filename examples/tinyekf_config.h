/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
 *
 * tinyekf_config.h: static configuration parameters
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

/* states */
#define N 8

/* observables */
#define M 4

typedef struct {

    int n;           /* number of state values */
    int m;           /* number of observables */

    double x[N];     /* state vector */

    double P[N][N];  /* prediction error covariance */
    double Q[N][N];  /* process noise covariance */
    double R[M][M];  /* measurement error covariance */

    double G[N][M];  /* Kalman gain; a.k.a. K */

    double F[N][N];  /* Jacobian of process model */
    double H[M][N];  /* Jacobian of measurement model */

    double Ht[N][M]; /* transpose of measurement Jacobian */
    double Ft[N][N]; /* transpose of process Jacobian */
    double Pp[N][N]; /* P, post-prediction, pre-update */

    double fx[N];   /* output of user defined f() state-transition function */
    double hx[M];   /* output of user defined h() measurement function */

    /* temporary storage */
    double tmp1[N][N];
    double tmp2[M][N];
    double tmp3[M][M];
    double tmp4[M][M];
    double tmp5[M]; 

} ekf_t;
