#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

static void choldc1(double ** a, double * p, int n) {
    int i,j,k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i][j];
            for (k = i - 1; k >= 0; k--) {
                sum -= a[i][k] * a[j][k];
            }
            if (i == j) {
                if (sum <= 0) {
                    printf(" a is not positive definite!\n");
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j][i] = sum / p[i];
            }
        }
    }
}

static void choldcsl(double ** A, double ** a, double * p, int n) 
{
    int i,j,k; double sum;
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            a[i][j] = A[i][j];
    choldc1(a, p, n);
    for (i = 0; i < n; i++) {
        a[i][i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j][k] * a[k][i];
            }
            a[j][i] = sum / p[j];
        }
    }
}


static void cholsl(double ** A, double ** a, double * p, int n) 
{
    int i,j,k;
    choldcsl(A,a,p,n);
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i][j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i][i] *= a[i][i];
        for (k = i + 1; k < n; k++) {
            a[i][i] += a[k][i] * a[k][i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i][j] += a[k][i] * a[k][j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i][j] = a[j][i];
        }
    }
}

#include "linalg.h"

void mat_init(mat_t * mat, int m, int n)
{
    int i;

    mat->data = (double **)malloc(m*sizeof(double *));

    for (i=0; i<m; ++i)
        mat->data[i] = (double *)malloc(n*sizeof(double));

    mat->rows = m;
    mat->cols = n;

    mat->tmp = (double *)malloc(m*sizeof(double));
}

void mat_free(mat_t mat)
{
    int i;

    for (i=0; i<mat.rows; ++i)
        free(mat.data[i]);

    free(mat.data);

    free(mat.tmp);
}

void zeros(mat_t mat, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            mat.data[i][j] = 0;
}

void eye(mat_t mat, int n, double s)
{
    zeros(mat, mat.rows, mat.cols);

    int k;

    for(k=0; k<n; ++k)
        mat.data[k][k] = s;
}

void vec_dump(double * x, int n, const char * fmt)
{
    int j;
    char f[100];
    sprintf(f, "%s ", fmt);
    for(j=0; j<n; ++j)
        printf(f, x[j]);
    printf("\n");
}

void mat_dump(mat_t mat, int m, int n, const char * fmt)
{
    int i,j;

    char f[100];
    sprintf(f, "%s ", fmt);
     for(i=0; i<m; ++i) {
        for(j=0; j<n; ++j)
            printf(f, mat.data[i][j]);
        printf("\n");
    }
}

// C <- A * B
void mulmat(mat_t a, mat_t b, mat_t c)
{
    int i, j,l;

    for(i=0; i<a.rows; ++i)
        for(j=0; j<b.cols; ++j) {
            c.data[i][j] = 0;
            for(l=0; l<a.cols; ++l)
                c.data[i][j] += a.data[i][l] * b.data[l][j];
        }
}

void mulvec(mat_t a, double * x, double * y)
{
    int i,j;

    for(i=0; i<a.rows; ++i) {
        y[i] = 0;
        for(j=0; j<a.cols; ++j)
            y[i] += x[j] * a.data[i][j];
    }
}

void transpose(mat_t a, mat_t at, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j) 
            at.data[j][i] = a.data[i][j];
}

// A <- A + B
void add(mat_t a, mat_t b, int m, int n)
{        
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a.data[i][j] += b.data[i][j];
}

// C <- A - B
void sub(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] - b[j];
}

void negate(mat_t a, int m, int n)
{        
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a.data[i][j] = -a.data[i][j];
}

void invert(mat_t a, mat_t at, int n)
{
    cholsl(a.data, at.data, a.tmp, n);
}

void mat_set(mat_t a, int i, int j, double value)
{
    a.data[i][j] = value;
}
