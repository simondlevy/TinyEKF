#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tinyekf.h"

static void choldc1(number_t * a, number_t * p, int n) {
    int i,j,k;
    number_t sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i*n+j];
            for (k = i - 1; k >= 0; k--) {
                sum -= a[i*n+k] * a[j*n+k];
            }
            if (i == j) {
                if (sum <= 0) {
                    printf(" a is not positive definite!\n");
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j*n+i] = sum / p[i];
            }
        }
    }
}

static void choldcsl(number_t * A, number_t * a, number_t * p, int n) 
{
    int i,j,k; number_t sum;
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            a[i*n+j] = A[i*n+j];
    choldc1(a, p, n);
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
}


static void cholsl(number_t * A, number_t * a, number_t * p, int n) 
{
    int i,j,k;
    choldcsl(A,a,p,n);
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
}

/*
static void vec_dump(number_t * x, int n, const char * fmt)
{
    int j;
    char f[100];
    sprintf(f, "%s ", fmt);
    for(j=0; j<n; ++j)
        printf(f, x[j]);
    printf("\n");
}
*/

static void mat_dump(number_t * a, int m, int n, const char * fmt)
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

// C <- A * B
static void mulmat(number_t * a, number_t * b, number_t * c, int arows, int acols, int bcols)
{
    int i, j,l;

    for(i=0; i<arows; ++i)
        for(j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for(l=0; l<acols; ++l)
                c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
        }
}

static void mulvec(number_t * a, number_t * x, number_t * y, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i) {
        y[i] = 0;
        for(j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

static void transpose(number_t * a, number_t * at, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j) 
            at[j*m+i] = a[i*n+j];
}

// A <- A + B
static void add(number_t * a, number_t * b, int m, int n)
{        
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] += b[i*n+j];
}

// C <- A - B
static void sub(number_t * a, number_t * b, number_t * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] - b[j];
}

static void negate(number_t * a, int m, int n)
{        
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] = -a[i*n+j];
}

static void invert(number_t * a, number_t * at, number_t * p, int n)
{
    cholsl(a, at, p, n);
}

static void mat_addeye(number_t * a, int n)
{
    int i;
    for (i=0; i<n; ++i)
        a[i*n+i] += 1;
}


// ----------------------------------------------------------

void ekf_step(
        ekf_t * ekf, 
        number_t * Z, 
        void (*f)(number_t x[N], number_t F[N][N]), 
        void (*h)(number_t x[N], number_t hx[N], number_t H[M][N]))
{        
    // \hat{x}_k = f(\hat{x}_{k-1})
    f(ekf->x, ekf->F);

    // h(\hat{x}_k
    h(ekf->x, ekf->hx, ekf->H);     

    // Predict and update
    ekf_predict_and_update(ekf, Z);
}

void ekf_predict_and_update(ekf_t * ekf, number_t * Z)
{    
    // P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}
    mulmat(&ekf->F[0][0], &ekf->P[0][0], ekf->tmp1, N, N, N);
    transpose(&ekf->F[0][0], &ekf->Ft[0][0], N, N);
    mulmat(ekf->tmp1, &ekf->Ft[0][0], &ekf->Pp[0][0], N, N, N);
    add(&ekf->Pp[0][0], &ekf->Q[0][0], N, N);

    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    transpose(&ekf->H[0][0], &ekf->Ht[0][0], M, N);
    mulmat(&ekf->Pp[0][0], &ekf->Ht[0][0], ekf->tmp1, N, N, M);
    mulmat(&ekf->H[0][0], &ekf->Pp[0][0], ekf->tmp2, M, N, N);
    mulmat(ekf->tmp2, &ekf->Ht[0][0], &ekf->tmp2_m_m[0][0], M, N, M);
    add(&ekf->tmp2_m_m[0][0], &ekf->R[0][0], M, M);
    invert(&ekf->tmp2_m_m[0][0], &ekf->tmp_m_m[0][0], ekf->tmp_m, M);
    mulmat(ekf->tmp1, &ekf->tmp_m_m[0][0], &ekf->G[0][0], N, M, M);

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k
    sub(ekf->tmp1, ekf->hx, Z, M);
    mulvec(&ekf->G[0][0], ekf->tmp1, &ekf->x[0], N, M);

    // P_k = (I - G_k H_k) P_k
    mulmat(&ekf->G[0][0], &ekf->H[0][0], ekf->tmp1, N, M, N);
    negate(ekf->tmp1, N, N);
    mat_addeye(ekf->tmp1, N);
    mulmat(ekf->tmp1, &ekf->Pp[0][0], &ekf->P[0][0], N, N, N);

    mat_dump(&ekf->P[0][0], N, N, "%+10.4f"); exit(0);
}
