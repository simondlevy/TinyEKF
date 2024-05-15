/*
 * static configuration parameters for GPS examples
 *
 * Copyright (C) 2024 Simon D. Levy
 *
 * MIT License
 */

// Size of state space
#define EKF_N 8

// Size of observation (measurement) space
#define EKF_M 4

// We need double-precision to replicate the published results
typedef double _float_t;
