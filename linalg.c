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

vec_t * newvec(int n)
{
    vec_t * vec = (vec_t *)malloc(sizeof(vec_t));

    vec->data = (double *)malloc(n*sizeof(double));

    vec->len = n;

    return vec;
}

mat_t * newmat(int m, int n)
{
    int i;

    mat_t * mat = (mat_t *)malloc(sizeof(mat_t));

    mat->data = (double **)malloc(m*sizeof(double *));

    for (i=0; i<m; ++i)
        mat->data[i] = (double *)malloc(n*sizeof(double));

    mat->rows = m;
    mat->cols = n;

    mat->tmp = (double *)malloc(m*sizeof(double));

    return mat;
}

void deletemat(mat_t * mat)
{
    int i;

    for (i=0; i<mat->rows; ++i)
        free(mat->data[i]);

    free(mat->data);

    free(mat->tmp);

    free(mat);
}

void deletevec(vec_t * vec)
{
    free(vec->data);

    free(vec);
}


void zeros(mat_t * mat)
{
    int i;

    for(i=0; i<mat->rows; ++i)
        bzero(mat->data[i], mat->cols*sizeof(double));
}

void eye(mat_t * mat, double s)
{
    zeros(mat);

    int k;

    for(k=0; k<mat->rows; ++k)
        mat->data[k][k] = s;
}

void dumpvec(vec_t * vec, const char * fmt)
{
    int j;
    char f[100];
    sprintf(f, "%s ", fmt);
    for(j=0; j<vec->len; ++j)
        printf(f, vec->data[j]);
    printf("\n");
}

void dumpmat(mat_t * mat, const char * fmt)
{
    int i,j;

    char f[100];
    sprintf(f, "%s ", fmt);
     for(i=0; i<mat->rows; ++i) {
        for(j=0; j<mat->cols; ++j)
            printf(f, mat->data[i][j]);
        printf("\n");
    }
}

// C <- A * B
void mulmat(mat_t * a, mat_t * b, mat_t * c)
{
    int i, j,l;

    for(i=0; i<a->rows; ++i)
        for(j=0; j<b->cols; ++j) {
            c->data[i][j] = 0;
            for(l=0; l<a->cols; ++l)
                c->data[i][j] += a->data[i][l] * b->data[l][j];
        }
}

void mulvec(mat_t * a, vec_t * x, vec_t * y)
{
    int i,j;

    for(i=0; i<a->rows; ++i) {
        y->data[i] = 0;
        for(j=0; j<a->cols; ++j)
            y->data[i] += x->data[j] * a->data[i][j];
    }
}

void transpose(mat_t * a, mat_t * at)
{
    int i,j;

    for(i=0; i<a->rows; ++i)
        for(j=0; j<a->cols; ++j) 
            at->data[j][i] = a->data[i][j];
}

// A <- A + B
void add(mat_t * a, mat_t * b)
{        
    int i,j;

    for(i=0; i<a->rows; ++i)
        for(j=0; j<a->cols; ++j)
            a->data[i][j] += b->data[i][j];
}

// A <- A - B
void sub(vec_t * a, vec_t * b)
{
    int j;

    for(j=0; j<a->len; ++j)
        a->data[j] -= b->data[j];
}

void negate(mat_t * a)
{        
    int i, j;

    for(i=0; i<a->rows; ++i)
        for(j=0; j<a->cols; ++j)
            a->data[i][j] = -a->data[i][j];
}

void invert(mat_t * a, mat_t * at)
{
    cholsl(a->data, at->data, a->tmp, a->rows);
}
