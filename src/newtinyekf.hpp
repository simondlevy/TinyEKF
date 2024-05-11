/*
 * Extended Kalman filter with separate prediction and update
 *
 * Copyright (C) 2024 Simon D. Levy
 *
 * MIT License
 */

#pragma once

#include <stdint.h>
#include <string.h>

class TinyEkf {

    public:

        void initialize(
                const float  pdiag[EKF_N],
                const uint32_t nowMsec,
                const float min_covariance, 
                const float max_covariance)
        {
            _isUpdated = false;
            _min_covariance = min_covariance;
            _max_covariance = max_covariance;

            for (uint8_t i=0; i<EKF_N; ++i) {
                for (uint8_t j=0; j<EKF_N; ++j) {
                    _p[i][j] = i==j ? pdiag[i] : 0;
                }
                _x[i] = 0;
            }
        }

        void predict(
                const float x[EKF_N],
                const float F[EKF_N][EKF_N], 
                const float Q[EKF_N][EKF_N])
        {
            _isUpdated = true;

            // $\hat{x}_k = f(\hat{x}_{k-1})$
            for (uint8_t i=0; i<EKF_N; ++i) {
                _x[i] = x[i];
            }

            // # $P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}$
            multiplyCovariance(F);
            add(_p, Q, _p);

            cleanupCovariance();
        }

        void update(
                const float h[EKF_N], 
                const float error, 
                const float stdMeasNoise)
        {

            float ph[EKF_N] = {};
            multiply(_p, h, ph);
            const auto r = stdMeasNoise * stdMeasNoise;
            const auto hphr = r + dot(h, ph); // HPH' + R

            float g[EKF_N] = {};
            for (uint8_t i=0; i<EKF_N; ++i) {
                g[i] = ph[i] / hphr;
            }

            // $\hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))$
            for (uint8_t i=0; i<EKF_N; ++i) {
                _x[i] += g[i] * error;
            }

            float GH[EKF_N][EKF_N] = {};
            outer(g, h, GH); 

            for (int i=0; i<EKF_N; i++) { 
                GH[i][i] -= 1;
            }

            // $P_k = (I - G_k H_k) P_k$
            multiplyCovariance(GH);

            // Add the measurement variance 
            for (int i=0; i<EKF_N; i++) {
                for (int j=0; j<EKF_N; j++) {
                    _p[i][j] += j < i ? 0 : r * g[i] * g[j];
                }
            }

            cleanupCovariance();

            _isUpdated = true;
        }

        void finalize(
                const float x[EKF_N], 
                const float A[EKF_N][EKF_N], 
                const bool isErrorSufficient)
        {
            if (_isUpdated) {

                if (isErrorSufficient) {

                    multiplyCovariance(A);

                    cleanupCovariance();
                }

                _isUpdated = false;
            }
        }

        float * getState(void)
        {
            static float x[EKF_N];

            for (uint8_t i=0; i<EKF_N; ++i) {
                x[i] = _x[i];
            }

            return x;
        }
};
