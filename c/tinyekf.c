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

#include <tinyekf.h>

/* Cholesky-decomposition matrix-inversion code, adapated from
   http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */

static int choldc1(double * a, double * p, int n) {
    int i,j,k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i*n+j];
            for (k = i - 1; k >= 0; k--) {
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

    return 0; /* success */
}

static int choldcsl(double * A, double * a, double * p, int n) 
{
    int i,j,k; double sum;
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

    return 0; /* success */
}


static int cholsl(double * A, double * a, double * p, int n) 
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

    return 0; /* success */
}

static void zeros(double * a, int m, int n)
{
    int j;
    for (j=0; j<m*n; ++j)
        a[j] = 0;
}

#ifdef DEBUG
static void dump(double * a, int m, int n, const char * fmt)
{
    int i,j;

    char f[100];
    sprintf(f, "%s ", fmt);
    for(i=0; i<m; ++i) {
        for(j=0; j<n; ++j)
            printf(f, a[i*n+j]);
        printf("\n");
    }
}
#endif

/* C <- A * B */
static void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols)
{
    int i, j,l;

    for(i=0; i<arows; ++i)
        for(j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for(l=0; l<acols; ++l)
                c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
        }
}

static void mulvec(double * a, double * x, double * y, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i) {
        y[i] = 0;
        for(j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

static void transpose(double * a, double * at, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

/* A <- A + B */
static void accum(double * a, double * b, int m, int n)
{        
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] += b[i*n+j];
}

/* C <- A + B */
static void add(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] + b[j];
}


/* C <- A - B */
static void sub(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] - b[j];
}

static void negate(double * a, int m, int n)
{        
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] = -a[i*n+j];
}

static void mat_addeye(double * a, int n)
{
    int i;
    for (i=0; i<n; ++i)
        a[i*n+i] += 1;
}

void ekf_init(ekf_t * ekf)
{

    /* zero-out matrices */
    zeros(ekf->P, EKF_N, EKF_N);
    zeros(ekf->Q, EKF_N, EKF_N);
    zeros(ekf->R, EKF_M, EKF_M);
    zeros(ekf->G, EKF_N, EKF_M);
    zeros(ekf->F, EKF_N, EKF_N);
    zeros(ekf->H, EKF_M, EKF_N);
}

int ekf_step(ekf_t * ekf, double * z)
{        
    /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
    mulmat(ekf->F, ekf->P, ekf->tmp0, EKF_N, EKF_N, EKF_N);
    transpose(ekf->F, ekf->Ft, EKF_N, EKF_N);
    mulmat(ekf->tmp0, ekf->Ft, ekf->Pp, EKF_N, EKF_N, EKF_N);
    accum(ekf->Pp, ekf->Q, EKF_N, EKF_N);

    /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
    transpose(ekf->H, ekf->Ht, EKF_M, EKF_N);
    mulmat(ekf->Pp, ekf->Ht, ekf->tmp1, EKF_N, EKF_N, EKF_M);
    mulmat(ekf->H, ekf->Pp, ekf->tmp2, EKF_M, EKF_N, EKF_N);
    mulmat(ekf->tmp2, ekf->Ht, ekf->tmp3, EKF_M, EKF_N, EKF_M);
    accum(ekf->tmp3, ekf->R, EKF_M, EKF_M);
    if (cholsl(ekf->tmp3, ekf->tmp4, ekf->tmp5, EKF_M)) return 1;
    mulmat(ekf->tmp1, ekf->tmp4, ekf->G, EKF_N, EKF_M, EKF_M);

    /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
    sub(z, ekf->hx, ekf->tmp5, EKF_M);
    mulvec(ekf->G, ekf->tmp5, ekf->tmp2, EKF_N, EKF_M);
    add(ekf->fx, ekf->tmp2, ekf->x, EKF_N);

    /* P_k = (I - G_k H_k) P_k */
    mulmat(ekf->G, ekf->H, ekf->tmp0, EKF_N, EKF_M, EKF_N);
    negate(ekf->tmp0, EKF_N, EKF_N);
    mat_addeye(ekf->tmp0, EKF_N);
    mulmat(ekf->tmp0, ekf->Pp, ekf->P, EKF_N, EKF_N, EKF_N);

    /* success */
    return 0;
}
