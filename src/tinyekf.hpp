/*
 * Extended Kalman filter with separate prediction and update
 *
 * Copyright (C) 2024 Simon D. Levy
 *
 * MIT License
 */

#pragma once

#include <stdint.h>

class TinyEkf {

    public:

        void initialize(
                const float pdiag[EKF_N],
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

        // # $P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}$
        void predict(
            const float xnew[EKF_N],
            const float F[EKF_N][EKF_N],
            const float Q[EKF_N][EKF_N])
        {
            // $\hat{x}_k = f(\hat{x}_{k-1})$
            for (uint8_t i=0; i<EKF_N; ++i) {
                _x[i] = xnew[i];
            }

            multiplyCovariance(F);
            add(_p, Q, _p);
            cleanupCovariance();
         }

        void update(
                const float h[EKF_M], 
                const float z[EKF_N],
                const float hx[EKF_N],
                const float r)
        {
        }

         void update(
                const float h[EKF_N], 
                const float z,
                const float hx,
                const float r)
        {
            float ph[EKF_N] = {};
            multiply(_p, h, ph);
            const auto hphr = r + dot(h, ph); // HPH' + R

            float g[EKF_N] = {};
            for (uint8_t i=0; i<EKF_N; ++i) {
                g[i] = ph[i] / hphr;
            }

            // $\hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))$
            for (uint8_t i=0; i<EKF_N; ++i) {
                _x[i] += g[i] * (z - hx);
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
                const float newx[EKF_N], 
                const float A[EKF_N][EKF_N],
                const bool isErrorSufficient) 
        {
            if (_isUpdated) {

                for (uint8_t i=0; i<EKF_N; ++i) {
                    _x[i] = newx[i];
                }

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

    private:

        static void transpose(
                const float a[EKF_N][EKF_N], float at[EKF_N][EKF_N]) 
        {
            for (uint8_t i=0; i<EKF_N; ++i) {
                for (uint8_t j=0; j<EKF_N; ++j) {
                    auto tmp = a[i][j];
                    at[i][j] = a[j][i];
                    at[j][i] = tmp;
                }
            }
        }

        static float dot(const float x[EKF_N], const float y[EKF_N]) 
        {
            float d = 0;

            for (uint8_t k=0; k<EKF_N; k++) {
                d += x[k] * y[k];
            }

            return d;
        }

        // Matrix + Matrix
        static void add(
                const float a[EKF_N][EKF_N], 
                const float b[EKF_N][EKF_N], 
                float c[EKF_N][EKF_N]) 
        {
            for (uint8_t i=0; i<EKF_N; i++) {

                for (uint8_t j=0; j<EKF_N; j++) {

                    c[i][j] = a[i][j] + b[i][j];
                }
            }
        }

        // Matrix * Matrix
        static void multiply(
                const float a[EKF_N][EKF_N], 
                const float b[EKF_N][EKF_N], 
                float c[EKF_N][EKF_N]) 
        {
            for (uint8_t i=0; i<EKF_N; i++) {

                for (uint8_t j=0; j<EKF_N; j++) {

                    c[i][j] = 0;

                    for (uint8_t k=0; k<EKF_N; k++) {

                        c[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        }

        // Matrix * Vector
        static void multiply(
                const float a[EKF_N][EKF_N], 
                const float x[EKF_N],
                float y[EKF_N])

        {
            for (uint8_t i=0; i<EKF_N; i++) {
                y[i] = 0;
                for (uint8_t j=0; j<EKF_N; j++) {
                    y[i] += a[i][j] * x[j];
                }
            }
        }

        // Outer product
        static void outer(
                const float x[EKF_N],
                const float y[EKF_N],
                float a[EKF_N][EKF_N]) 
        {
            for (uint8_t i=0; i<EKF_N; i++) {
                for (uint8_t j=0; j<EKF_N; j++) {
                    a[i][j] = x[i] * y[j];
                }
            }
        }


        float _min_covariance;
        float _max_covariance;

        float _p[EKF_N][EKF_N];

        float _x[EKF_N];

        bool _isUpdated;

        void multiplyCovariance(const float a[EKF_N][EKF_N])
        {
            float at[EKF_N][EKF_N] = {};
            transpose(a, at);  
            float ap[EKF_N][EKF_N] = {};
            multiply(a, _p, ap);
            multiply(ap, at, _p);
        }


        void cleanupCovariance(void)
        {
            // Enforce symmetry of the covariance matrix, and ensure the
            // values stay bounded
            for (int i=0; i<EKF_N; i++) {

                for (int j=i; j<EKF_N; j++) {

                    const auto pval = (_p[i][j] + _p[j][i]) / 2;

                    _p[i][j] = _p[j][i] = 
                        pval > _max_covariance ?  _max_covariance :
                        (i==j && pval < _min_covariance) ?  _min_covariance :
                        pval;
                }
            }
        }
};
