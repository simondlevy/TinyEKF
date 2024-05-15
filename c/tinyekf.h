/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */


void ekf_init(void * ekf);

int ekf_step(void * ekf, double * z);
