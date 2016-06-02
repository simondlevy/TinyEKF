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
        double tmp1[N][N];
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
