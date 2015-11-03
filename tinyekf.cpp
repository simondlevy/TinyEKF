#include "tinyekf.hpp"

#include <stdio.h>
#include <math.h>
#include <strings.h>

static void choldc1(double * a, double * p, int n) {
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

static void choldcsl(double * A, double * a, double * p, int n) 
{
    int i,j,k; double sum;
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


static void cholsl(double * A, double * a, double * p, int n) 
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

// C <- A * B
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

// A <- A + B
static void accum(double * a, double * b, int m, int n)
{        
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] += b[i*n+j];
}

// C <- A + B
static void add(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] + b[j];
}


// C <- A - B
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

static void invert(double * a, double * at, double * p, int n)
{
    cholsl(a, at, p, n);
}

static void mat_addeye(double * a, int n)
{
    int i;
    for (i=0; i<n; ++i)
        a[i*n+i] += 1;
}


// -------------------------------------------------------------------

TinyEKF::TinyEKF(int n, int m) 
{
    this->n = n;
    this->m = m;

    bzero(this->P, N*N*sizeof(double)); 
    bzero(this->Q, N*N*sizeof(double)); 
    bzero(this->R, M*M*sizeof(double)); 
    bzero(this->G, N*M*sizeof(double)); 
    bzero(this->F, N*N*sizeof(double)); 
    bzero(this->H, M*N*sizeof(double)); 
}

TinyEKF::~TinyEKF()
{
}

void TinyEKF::setP(int i, int j, double value)
{
    this->P[i][j] = value;
}

void TinyEKF::setQ(int i, int j, double value)
{
    this->Q[i][j] = value;
}

void TinyEKF::setR(int i, int j, double value)
{
    this->R[i][j] = value;
}

void TinyEKF::setX(int i, double value)
{
    this->x[i] = value;
}

double TinyEKF::getX(int i)
{
    return this->x[i];
}

void TinyEKF::step(double * Z)
{        
    // Model
    this->f(this->x, this->fx, this->F); 
    this->h(this->fx, this->hx, this->H);     

    // P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}
    mulmat(&this->F[0][0], &this->P[0][0], this->tmp1, N, N, N);
    transpose(&this->F[0][0], &this->Ft[0][0], N, N);
    mulmat(this->tmp1, &this->Ft[0][0], &this->Pp[0][0], N, N, N);
    accum(&this->Pp[0][0], &this->Q[0][0], N, N);

    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    transpose(&this->H[0][0], &this->Ht[0][0], M, N);
    mulmat(&this->Pp[0][0], &this->Ht[0][0], this->tmp1, N, N, M);
    mulmat(&this->H[0][0], &this->Pp[0][0], this->tmp2, M, N, N);
    mulmat(this->tmp2, &this->Ht[0][0], this->tmp3, M, N, M);
    accum(this->tmp3, &this->R[0][0], M, M);
    invert(this->tmp3, this->tmp4, this->tmp5, M);
    mulmat(this->tmp1, this->tmp4, &this->G[0][0], N, M, M);

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k
    sub(Z, this->hx, this->tmp1, M);
    mulvec(&this->G[0][0], this->tmp1, this->tmp2, N, M);
    add(this->fx, this->tmp2, this->x, N);

    // P_k = (I - G_k H_k) P_k
    mulmat(&this->G[0][0], &this->H[0][0], this->tmp1, N, M, N);
    negate(this->tmp1, N, N);
    mat_addeye(this->tmp1, N);
    mulmat(this->tmp1, &this->Pp[0][0], &this->P[0][0], N, N, N);
}


