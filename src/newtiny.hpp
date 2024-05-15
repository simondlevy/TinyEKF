/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include <tinyekf.h>

class NewTinyEKF {

    public:

        void initialize(
                const float pdiag[EKF_N],
                const float min_covariance=0, 
                const float max_covariance=0)
        {
        }

        void predict(
                const float fx[EKF_N],
                const float F[EKF_N*EKF_N],
                const float Q[EKF_N*EKF_N])
        { 
        }

        bool update(
                const float z[EKF_M], 
                const float hx[EKF_M],
                const float H[EKF_M*EKF_N],
                const float R[EKF_M*EKF_M])
        { 
            // success
            return true;
        }

         void update_with_scalar(
                const float h[EKF_N], 
                const float z,
                const float hx,
                const float r)
        {
        }


        void finalize(
                const float newx[EKF_N], 
                const float A[EKF_N*EKF_N],
                const bool isErrorSufficient) 
        {
        }

        float * get(void) 
        { 
            return NULL;
        }
};
