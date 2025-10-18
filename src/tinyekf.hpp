/*
 * Extended Kalman Filter class for embedded processors
 *
 * Copyright (C) 2025 Simon D. Levy
 *
 * MIT License
 */

#include <matrix_typedef.h>

class TinyEkf {

    public:

        void init(
                const float pinit[TINYEKF_N],
                const float min_covariance,
                const float max_covariance)
        {
            for (uint8_t k=0; k<TINYEKF_N; ++k) {
                _p[k][k] = pinit[k] * pinit[k];
            }

            _min_covariance = min_covariance;
            _max_covariance = max_covariance;

            _p_m.numRows = TINYEKF_N;
            _p_m.numCols = TINYEKF_N;
            _p_m.pData = (float*)_p;
        }

        void addProcessNoise(const float pnoise[TINYEKF_N])
        {
            for (uint8_t k=0; k<TINYEKF_N; ++k) {
                _p[k][k] += pnoise[k] * pnoise[k];
            }

            enforceSymmetry();
        }

    private:

        // State vector
        __attribute__((aligned(4))) float _x[TINYEKF_N];

        // Covariance matrix
        __attribute__((aligned(4))) float _p[TINYEKF_N][TINYEKF_N];

        // Workspace
        matrix_t _p_m;

        float _min_covariance;
        float _max_covariance;

        void enforceSymmetry()
        {
            for (int i=0; i<TINYEKF_N; i++) {
                for (int j=i; j<TINYEKF_N; j++) {
                    float p = 0.5f*_p[i][j] + 0.5f*_p[j][i];
                    if (isnan(p) || p > _max_covariance) {
                        _p[i][j] = _p[j][i] = _max_covariance;
                    } else if ( i==j && p < _min_covariance ) {
                        _p[i][j] = _p[j][i] = _min_covariance;
                    } else {
                        _p[i][j] = _p[j][i] = p;
                    }
                }
            }

        }

};

