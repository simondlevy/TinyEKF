/*
 * Extended Kalman Filter for embedded processors
 *
 * Copyright (C) 2024 Simon D. Levy
 *
 * MIT License
 */

#include <math.h>
#include <stdbool.h>
#include <string.h>

/**
  * Floating-point precision defaults to single but can be made double via
    <tt><b>#define _float_t double</b></tt> before <tt>#include <tinyekf.h></tt>
  */
#ifndef _float_t
#define _float_t float
#endif

// Linear alegbra ////////////////////////////////////////////////////////////

/// @private
static void _mulmat(
        const _float_t * a, 
        const _float_t * b, 
        _float_t * c, 
        const int arows, 
        const int acols, 
        const int bcols)
{
    for (int i=0; i<arows; ++i) {
        for (int j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for (int k=0; k<acols; ++k) {
                c[i*bcols+j] += a[i*acols+k] * b[k*bcols+j];
            }
        }
    }
}

/// @private
static void _mulvec(
        const _float_t * a, 
        const _float_t * x, 
        _float_t * y, 
        const int m, 
        const int n)
{
    for (int i=0; i<m; ++i) {
        y[i] = 0;
        for (int j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

/// @private
static void _transpose(
        const _float_t * a, _float_t * at, const int m, const int n)
{
    for (int i=0; i<m; ++i)
        for (int j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

/// @private
static void _addmat(
        const _float_t * a, const _float_t * b, _float_t * c, 
        const int m, const int n)
{
    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            c[i*n+j] = a[i*n+j] + b[i*n+j];
        }
    }
}

/// @private
static void _negate(_float_t * a, const int m, const int n)
{        
    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            a[i*n+j] = -a[i*n+j];
        }
    }
}

/// @private
static void _addeye(_float_t * a, const int n)
{
    for (int i=0; i<n; ++i) {
        a[i*n+i] += 1;
    }
}


/* Cholesky-decomposition matrix-inversion code, adapated from
http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/_choles_cpp.txt */

/// @private
static int _choldc1(_float_t * a, _float_t * p, const int n) 
{
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            _float_t sum = a[i*n+j];
            for (int k = i - 1; k >= 0; k--) {
                sum -= a[i*n+k] * a[j*n+k];
            }
            if (i == j) {
                if (sum <= 0) {
                    return 1; /* error */
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j*n+i] = sum / p[i];
            }
        }
    }

    return 0; // success:w
}

/// @private
static int _choldcsl(const _float_t * A, _float_t * a, _float_t * p, const int n) 
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i*n+j] = A[i*n+j];
        }
    }
    if (_choldc1(a, p, n)) {
        return 1;
    }
    for (int i = 0; i < n; i++) {
        a[i*n+i] = 1 / p[i];
        for (int j = i + 1; j < n; j++) {
            _float_t sum = 0;
            for (int k = i; k < j; k++) {
                sum -= a[j*n+k] * a[k*n+i];
            }
            a[j*n+i] = sum / p[j];
        }
    }

    return 0; // success
}

/// @private
static int _cholsl(const _float_t * A, _float_t * a, _float_t * p, const int n) 
{
    if (_choldcsl(A,a,p,n)) {
        return 1;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            a[i*n+j] = 0.0;
        }
    }
    for (int i = 0; i < n; i++) {
        a[i*n+i] *= a[i*n+i];
        for (int k = i + 1; k < n; k++) {
            a[i*n+i] += a[k*n+i] * a[k*n+i];
        }
        for (int j = i + 1; j < n; j++) {
            for (int k = j; k < n; k++) {
                a[i*n+j] += a[k*n+i] * a[k*n+j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            a[i*n+j] = a[j*n+i];
        }
    }

    return 0; // success
}

/// @private
static void _addvec(
        const _float_t * a, const _float_t * b, _float_t * c, const int n)
{
    for (int j=0; j<n; ++j) {
        c[j] = a[j] + b[j];
    }
}

/// @private
static void _sub(
        const _float_t * a, const _float_t * b, _float_t * c, const int n)
{
    for (int j=0; j<n; ++j) {
        c[j] = a[j] - b[j];
    }
}

/// @private
static bool invert(const _float_t * a, _float_t * ainv)
{
    _float_t tmp[EKF_M];

    return _cholsl(a, ainv, tmp, EKF_M) == 0;
}


// EKF ///////////////////////////////////////////////////////////////////////

typedef struct {

    /** State vector **/
    _float_t x[EKF_N];        

    /** Prediction error covariance **/
    _float_t P[EKF_N*EKF_N];  

} ekf_t;

/**
 * Initializes the EKF
 * @param ekf pointer to an ekf_t structure
 * @param pdiag a vector of length EKF_N containing the initial values for the
 * covariance matrix diagonal
 */
static void ekf_initialize(ekf_t * ekf, const _float_t pdiag[EKF_N])
{
    for (int i=0; i<EKF_N; ++i) {

        for (int j=0; j<EKF_N; ++j) {

            ekf->P[i*EKF_N+j] = i==j ? pdiag[i] : 0;
        }

        ekf->x[i] = 0;
    }
}

/**
  * Runs the EKF prediction step
  * @param ekf pointer to an ekf_t structure
  * @param fx predicted values
  * @param F Jacobian of state-transition function
  * @param Q process noise matrix
  * 
  */static void ekf_predict(
        ekf_t * ekf, 
        const _float_t fx[EKF_N],
        const _float_t F[EKF_N*EKF_N],
        const _float_t Q[EKF_N*EKF_N])
{        
    // \hat{x}_k = f(\hat{x}_{k-1}, u_k)
    memcpy(ekf->x, fx, EKF_N*sizeof(_float_t));

    // P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}

    _float_t FP[EKF_N*EKF_N] = {};
    _mulmat(F, ekf->P,  FP, EKF_N, EKF_N, EKF_N);

    _float_t Ft[EKF_N*EKF_N] = {};
    _transpose(F, Ft, EKF_N, EKF_N);

    _float_t FPFt[EKF_N*EKF_N] = {};
    _mulmat(FP, Ft, FPFt, EKF_N, EKF_N, EKF_N);

    _addmat(FPFt, Q, ekf->P, EKF_N, EKF_N);
}

/// @private
static void ekf_update_step3(ekf_t * ekf, _float_t GH[EKF_N*EKF_N])
{
    _negate(GH, EKF_N, EKF_N);
    _addeye(GH, EKF_N);
    _float_t GHP[EKF_N*EKF_N];
    _mulmat(GH, ekf->P, GHP, EKF_N, EKF_N, EKF_N);
    memcpy(ekf->P, GHP, EKF_N*EKF_N*sizeof(_float_t));
}

/**
  * Runs the EKF update step
  * @param ekf pointer to an ekf_t structure
  * @param z observations
  * @param hx predicted values
  * @param H sensor-function Jacobian matrix
  * @param R measurement-noise matrix
  * 
  */
static bool ekf_update(
        ekf_t * ekf, 
        const _float_t z[EKF_M], 
        const _float_t hx[EKF_N],
        const _float_t H[EKF_M*EKF_N],
        const _float_t R[EKF_M*EKF_M])
{        
    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    _float_t G[EKF_N*EKF_M];
    _float_t Ht[EKF_N*EKF_M];
    _transpose(H, Ht, EKF_M, EKF_N);
    _float_t PHt[EKF_N*EKF_M];
    _mulmat(ekf->P, Ht, PHt, EKF_N, EKF_N, EKF_M);
    _float_t HP[EKF_M*EKF_N];
    _mulmat(H, ekf->P, HP, EKF_M, EKF_N, EKF_N);
    _float_t HpHt[EKF_M*EKF_M];
    _mulmat(HP, Ht, HpHt, EKF_M, EKF_N, EKF_M);
    _float_t HpHtR[EKF_M*EKF_M];
    _addmat(HpHt, R, HpHtR, EKF_M, EKF_M);
    _float_t HPHtRinv[EKF_M*EKF_M];
    if (!invert(HpHtR, HPHtRinv)) {
        return false;
    }
    _mulmat(PHt, HPHtRinv, G, EKF_N, EKF_M, EKF_M);

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
    _float_t z_hx[EKF_M];
    _sub(z, hx, z_hx, EKF_M);
    _float_t Gz_hx[EKF_M*EKF_N];
    _mulvec(G, z_hx, Gz_hx, EKF_N, EKF_M);
    _addvec(ekf->x, Gz_hx, ekf->x, EKF_N);

    // P_k = (I - G_k H_k) P_k
    _float_t GH[EKF_N*EKF_N];
    _mulmat(G, H, GH, EKF_N, EKF_M, EKF_N);
    ekf_update_step3(ekf, GH);

    // success
    return true;
}
