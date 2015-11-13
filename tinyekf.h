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
  * Sets contents to zero.
  * @param ekf pointer to EKF structure to initialize
  * @param n number of state variables
  * @param m number of observables
  */
void ekf_init(void * ekf, int n, int m);

/**
  * Runs one step of prediction and update.
  * @param ekf pointer to structure EKF 
  * @param z array of measurement (observation) values
  * @return 0 on success, 1 on failure caused by non-positive-definite matrix.
  */
int ekf_step(void * ekf, double * z);
