/*
 * TinyEKF: Extended Kalman Filter for embedded processors
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

/* C <- A * B */
static void mulmat(
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

static void mulvec(
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

static void transpose(
        const _float_t * a, _float_t * at, const int m, const int n)
{
    for (int i=0; i<m; ++i)
        for (int j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

typedef struct {

    _float_t x[EKF_N];         // state vector
    _float_t P[EKF_N*EKF_N];  // prediction error covariance
    _float_t Pp[EKF_N*EKF_N]; // P, post-prediction, pre-update

} ekf_t;

static void ekf_initialize(ekf_t * ekf, const _float_t pdiag[EKF_N])
{
    memset(ekf->P, 0, EKF_M*EKF_N*sizeof(_float_t));

    for (int i=0; i<EKF_N; ++i) {
        ekf->P[i*EKF_N+i] = pdiag[i];
    }
}


//////////////////////////////////////////////////////////////////////////////

#ifndef _FOO

/* C <- A + B */
static void addvec(
        const _float_t * a, const _float_t * b, _float_t * c, const int n)
{
    for (int j=0; j<n; ++j) {
        c[j] = a[j] + b[j];
    }
}

static void addmat(
        const _float_t * a, const _float_t * b, _float_t * c, 
        const int m, const int n)
{
    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            c[i*n+j] = a[i*n+j] + b[i*n+j];
        }
    }
}


/* Cholesky-decomposition matrix-inversion code, adapated from
http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */

static int choldc1(_float_t * a, _float_t * p, const int n) 
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

static int choldcsl(const _float_t * A, _float_t * a, _float_t * p, const int n) 
{
    int i,j,k; _float_t sum;
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            a[i*n+j] = A[i*n+j];
    if (choldc1(a, p, n)) return 1;
    for (i = 0; i < n; i++) {
        a[i*n+i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j*n+k] * a[k*n+i];
            }
            a[j*n+i] = sum / p[j];
        }
    }

    return 0; // success
}


static int cholsl(const _float_t * A, _float_t * a, _float_t * p, const int n) 
{
    int i,j,k;
    if (choldcsl(A,a,p,n)) return 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i*n+j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i*n+i] *= a[i*n+i];
        for (k = i + 1; k < n; k++) {
            a[i*n+i] += a[k*n+i] * a[k*n+i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i*n+j] += a[k*n+i] * a[k*n+j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i*n+j] = a[j*n+i];
        }
    }

    return 0; // success
}

/* C <- A - B */
static void sub(
        const _float_t * a, const _float_t * b, _float_t * c, const int n)
{
    for (int j=0; j<n; ++j) {
        c[j] = a[j] - b[j];
    }
}

static void negate(_float_t * a, const int m, const int n)
{        
    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            a[i*n+j] = -a[i*n+j];
        }
    }
}

static void mat_addeye(_float_t * a, const int n)
{
    for (int i=0; i<n; ++i) {
        a[i*n+i] += 1;
    }
}

//////////////////////////////////////////////////////////////////////////////


static void ekf_predict(
        ekf_t * ekf, 
        const _float_t fx[EKF_N],
        const _float_t F[EKF_N*EKF_N],
        const _float_t Q[EKF_N*EKF_N])
{        
    // \hat{x}_k = f(\hat{x}_{k-1}, u_k)
    memcpy(ekf->x, fx, EKF_N*sizeof(_float_t));

    // temporary storage
    _float_t tmp0[EKF_N*EKF_N];

    _float_t Ft[EKF_N*EKF_N]; // transpose of process Jacobian

    // P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}
    mulmat(F, ekf->P, tmp0, EKF_N, EKF_N, EKF_N);
    transpose(F, Ft, EKF_N, EKF_N);
    mulmat(tmp0, Ft, ekf->Pp, EKF_N, EKF_N, EKF_N);
    addmat(ekf->Pp, Q, ekf->Pp, EKF_N, EKF_N);
}

static bool ekf_update(
        ekf_t * ekf, 
        const _float_t z[EKF_M], 
        const _float_t hx[EKF_N],
        const _float_t H[EKF_M*EKF_N],
        const _float_t R[EKF_M*EKF_M])
{        
    /* temporary storage */
    _float_t tmp0[EKF_N*EKF_N];
    _float_t tmp1[EKF_N*EKF_M];
    _float_t tmp2[EKF_M*EKF_N];
    _float_t tmp3[EKF_M*EKF_M];
    _float_t tmp4[EKF_M*EKF_M];
    _float_t tmp5[EKF_M]; 

    _float_t G[EKF_N*EKF_M];  // Kalman gain; a.k.a. K
    _float_t Ht[EKF_N*EKF_M]; // transpose of measurement Jacobian

    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    transpose(H, Ht, EKF_M, EKF_N);
    mulmat(ekf->Pp, Ht, tmp1, EKF_N, EKF_N, EKF_M);
    mulmat(H, ekf->Pp, tmp2, EKF_M, EKF_N, EKF_N);
    mulmat(tmp2, Ht, tmp3, EKF_M, EKF_N, EKF_M);
    addmat(tmp3, R, tmp3, EKF_M, EKF_M);
    if (cholsl(tmp3, tmp4, tmp5, EKF_M)) return false;
    mulmat(tmp1, tmp4, G, EKF_N, EKF_M, EKF_M);

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
    sub(z, hx, tmp5, EKF_M);
    mulvec(G, tmp5, tmp2, EKF_N, EKF_M);
    addvec(ekf->x, tmp2, ekf->x, EKF_N);

    // P_k = (I - G_k H_k) P_k
    mulmat(G, H, tmp0, EKF_N, EKF_M, EKF_N);
    negate(tmp0, EKF_N, EKF_N);
    mat_addeye(tmp0, EKF_N);
    mulmat(tmp0, ekf->Pp, ekf->P, EKF_N, EKF_N, EKF_N);

    // success
    return true;
}
#endif
