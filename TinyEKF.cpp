/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
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
 * along with this code.  If not, see <http:#www.gnu.org/licenses/>.
 */

#include "TinyEKF.h"

#include <math.h>

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
                    return 1; // error
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j*n+i] = sum / p[i];
            }
        }
    }

    return 0; // success
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

    return 0; // success
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

    return 0; // success
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

    this->x = new double [n];
    this->fx = new double [n];
    this->hx = new double [n];

    this->F = new double [n*n];
    this->H = new double [m*n];
    this->P = new double [n*n];
    this->Q = new double [n*n];
    this->R = new double [m*m];
    this->G = new double [n*m];
    this->Ft = new double [n*n];
    this->Ht = new double [n*m];
    this->Pp = new double [n*n];

    this->tmp1 = new double[n*n];
    this->tmp2 = new double[m*n];
    this->tmp3 = new double[m*m];
    this->tmp4 = new double[m*m];
    this->tmp5 = new double[m];

    zeros(this->P, n, n);
    zeros(this->Q, n, n);
    zeros(this->R, m, m);
    zeros(this->G, n, m);
    zeros(this->F, n, n);
    zeros(this->H, m, n);
}

TinyEKF::~TinyEKF()
{
    delete this->x;
    delete this->fx;
    delete this->hx;

    delete this->F;
    delete this->H;
    delete this->P;
    delete this->Q;
    delete this->R;
    delete this->G;
    delete this->Ft;
    delete this->Ht;
    delete this->Pp;

    delete this->tmp1;
    delete this->tmp2;
    delete this->tmp3;
    delete this->tmp4;
    delete this->tmp5;

}

void TinyEKF::setP(int i, int j, double value)
{
    this->P[i*this->n+j] = value;
}

void TinyEKF::setQ(int i, int j, double value)
{
    this->Q[i*this->n+j] = value;
}

void TinyEKF::setR(int i, int j, double value)
{
    this->R[i*this->m+j] = value;
}

double TinyEKF::getX(int i)
{
    return this->x[i];
}

bool TinyEKF::step(double * z)
{        
    // Model
    this->model(this->fx, this->F, this->hx, this->H);

    // P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}
    mulmat(this->F, this->P, this->tmp1, this->n, this->n, this->n);
    transpose(this->F, this->Ft, this->n, this->n);
    mulmat(this->tmp1, this->Ft, this->Pp, this->n, this->n, this->n);
    accum(this->Pp, this->Q, this->n, this->n);

    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    transpose(this->H, this->Ht, this->m, this->n);
    mulmat(this->Pp, this->Ht, this->tmp1, this->n, this->n, this->m);
    mulmat(this->H, this->Pp, this->tmp2, this->m, this->n, this->n);
    mulmat(this->tmp2, this->Ht, this->tmp3, this->m, this->n, this->m);
    accum(this->tmp3, this->R, this->m, this->m);
    if (cholsl(this->tmp3, this->tmp4, this->tmp5, this->m)) return false;
    mulmat(this->tmp1, this->tmp4, this->G, this->n, this->m, this->m);

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k
    sub(z, this->hx, this->tmp1, this->m);
    mulvec(this->G, this->tmp1, this->tmp2, this->n, this->m);
    add(this->fx, this->tmp2, this->x, this->n);

    // P_k = (I - G_k H_k) P_k
    mulmat(this->G, this->H, this->tmp1, this->n, this->m, this->n);
    negate(this->tmp1, this->n, this->n);
    mat_addeye(this->tmp1, this->n);
    mulmat(this->tmp1, this->Pp, this->P, this->n, this->n, this->n);

    // success
    return true;
}

void TinyEKF::set(double * A, int i, int j, double value)
{
    A[i*this->n+j] = value;
}

