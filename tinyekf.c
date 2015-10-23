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

static void zeros(number_t * a, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] = 0;
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

static void mat_set(number_t * a, int i, int j, int n, number_t value)
{
    a[i*n+j] = value;
}

static void mat_addeye(number_t * a, int n)
{
    for (int i=0; i<n; ++i)
        a[i*n+i] += 1;
}


// ----------------------------------------------------------

void ekf_setP(ekf_t * ekf, int i, int j, number_t value)
{
    mat_set(ekf->P, i, j, N, value);
}

void ekf_setQ(ekf_t * ekf, int i, int j, number_t value)
{
    mat_set(ekf->Q, i, j, N, value);
}

void ekf_setR(ekf_t * ekf, int i, int j, number_t value)
{
    mat_set(ekf->R, i, j, M, value);
}

void ekf_setX(ekf_t * ekf, int i, number_t value)
{
    ekf->X[i] = value;
}

static void ekf_pre_update(
        ekf_t * ekf, 
        void (*f)(number_t *, number_t *, number_t *), 
        void (*g)(number_t *, number_t *, number_t *))
{
    // 1, 2
    zeros(ekf->fy, N, N);
    f(ekf->X, ekf->Xp, ekf->fy);

    // 3
    zeros(ekf->H, M, N);
    g(ekf->Xp, ekf->gXp, ekf->H);     
}

void ekf_update(
        ekf_t * ekf, 
        number_t * Z, 
        void (*f)(number_t *, number_t *, number_t *), 
        void (*g)(number_t *, number_t *, number_t *))
{        
    // 1,2,3
    ekf_pre_update(ekf, f, g);

    // 4,5,6,7
    ekf_post_update(ekf, Z);
}

void ekf_post_update(ekf_t * ekf, number_t * Z)
{    
    // 4
    mulmat(ekf->fy, ekf->P, ekf->tmp_n_n, N, N, N);
    transpose(ekf->fy, ekf->fyt, N, N);
    mulmat(ekf->tmp_n_n, ekf->fyt, ekf->Pp, N, N, N);
    add(ekf->Pp, ekf->Q, N, N);

    // 5
    transpose(ekf->H, ekf->Ht, M, N);
    mulmat(ekf->Pp, ekf->Ht, ekf->tmp_n_m, N, N, M);
    mulmat(ekf->H, ekf->Pp, ekf->tmp_m_n, M, N, N);
    mulmat(ekf->tmp_m_n, ekf->Ht, ekf->tmp2_m_m, M, N, M);
    add(ekf->tmp2_m_m, ekf->R, M, M);
    invert(ekf->tmp2_m_m, ekf->tmp_m_m, ekf->tmp_m, M);
    mulmat(ekf->tmp_n_m, ekf->tmp_m_m, ekf->G, N, M, M);

    // 6
    sub(ekf->tmp_m, ekf->gXp, Z, M);
    mulvec(ekf->G, ekf->tmp_m, ekf->X, N, M);

    // 7
    mulmat(ekf->G, ekf->H, ekf->tmp_n_n, N, M, N);
    negate(ekf->tmp_n_n, N, N);
    mat_addeye(ekf->tmp_n_n, N);
    mulmat(ekf->tmp_n_n, ekf->Pp, ekf->P, N, N, N);

    mat_dump(ekf->P, N, N, "%+10.4f"); exit(0);
}
