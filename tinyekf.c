/*
 * TinyEKF: Extended Kalman Filter for embedded processors
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * This code is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ekf-> code.  If not, see <http:#www.gnu.org/licenses/>.
 */

#include "tinyekf.h"

#include <math.h>
#include <stdlib.h>

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

/*
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
*/

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
        for(j=0; j<n; ++j) 
            at[j*m+i] = a[i*n+j];
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

/* TinyEKF code ------------------------------------------------------------------- */

void ekf_init(ekf_t * ekf)
{
    zeros(&ekf->P[0][0], N, N);
    zeros(&ekf->Q[0][0], N, N);
    zeros(&ekf->R[0][0], M, M);
    zeros(&ekf->G[0][0], N, M);
    zeros(&ekf->F[0][0], N, N);
    zeros(&ekf->H[0][0], M, N);
}

double ekf_getX(ekf_t * ekf, int i)
{
    return ekf->x[i];
}

int ekf_step(ekf_t * ekf, double * z)
{        
    /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
    mulmat(&ekf->F[0][0], &ekf->P[0][0], &ekf->tmp1[0][0], N, N, N);
    transpose(&ekf->F[0][0], &ekf->Ft[0][0], N, N);
    mulmat(&ekf->tmp1[0][0], &ekf->Ft[0][0], &ekf->Pp[0][0], N, N, N);
    accum(&ekf->Pp[0][0], &ekf->Q[0][0], N, N);

    /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
    transpose(&ekf->H[0][0], &ekf->Ht[0][0], M, N);
    mulmat(&ekf->Pp[0][0], &ekf->Ht[0][0], &ekf->tmp1[0][0], N, N, M);
    mulmat(&ekf->H[0][0], &ekf->Pp[0][0], &ekf->tmp2[0][0], M, N, N);
    mulmat(&ekf->tmp2[0][0], &ekf->Ht[0][0], &ekf->tmp3[0][0], M, N, M);
    accum(&ekf->tmp3[0][0], &ekf->R[0][0], M, M);
    if (cholsl(&ekf->tmp3[0][0], &ekf->tmp4[0][0], ekf->tmp5, M)) return 1;
    mulmat(&ekf->tmp1[0][0], &ekf->tmp4[0][0], &ekf->G[0][0], N, M, M);

    /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k */
    sub(z, ekf->hx, &ekf->tmp1[0][0], M);
    mulvec(&ekf->G[0][0], &ekf->tmp1[0][0], &ekf->tmp2[0][0], N, M);
    add(ekf->fx, &ekf->tmp2[0][0], ekf->x, N);

    /* P_k = (I - G_k H_k) P_k */
    mulmat(&ekf->G[0][0], &ekf->H[0][0], &ekf->tmp1[0][0], N, M, N);
    negate(&ekf->tmp1[0][0], N, N);
    mat_addeye(&ekf->tmp1[0][0], N);
    mulmat(&ekf->tmp1[0][0], &ekf->Pp[0][0], &ekf->P[0][0], N, N, N);

    /* success */
    return 0;
}
