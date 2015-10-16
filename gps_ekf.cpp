/*
 References:
 1. R G Brown, P Y C Hwang, "Introduction to random signals and applied 
   Kalman filtering : with MATLAB exercises and solutions",1996
 2. Pratap Misra, Per Enge, "Global Positioning System Signals, 
   Measurements, and Performance(Second Edition)",2006
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strings.h>

#include "tinyekf.hpp"

class GPS_EKF : TinyEKF {

    public:

        GPS_EKF(int m, int n) : TinyEKF(m, n) 
        {
        }

    protected:

        void f(float * x, float * fx, float * dfx)
        {
        }

        void g(float * x, float * fx, float * dfx)
        {
        }

};

static char * readline(char * line, FILE * fp)
{
    return fgets(line, 1000, fp);
}

static void fill(char * line, double * SV_Pos, double * SV_Rho)
{
    char * p = strtok(line, ",");

    for (int k=0; k<12; ++k) {
        SV_Pos[k] = atof(p);
        printf("%2d: %f\n", k, SV_Pos[k]);
        p = strtok(NULL, ",");
    }
    printf("\n");
    for (int k=0; k<4; ++k) {
        SV_Rho[k] = atof(p);
        printf("%2d: %f\n", k, SV_Rho[k]);
        p = strtok(NULL, ",");
    }
    printf("-----\n");
}

static void dumpmat(double * a, int m, int n)
{
    for (int j=0; j<m; ++j) {
        for (int k=0; k<n; ++k)
            printf("%10.6e ", a[j*n+k]);
        printf("\n");
    }
}

static void blkfill(double * out, const double * a, int off)
{
    int k = off*18;
    
    out[k]   = a[0];
    out[k+1] = a[1];
    out[k+8] = a[2];
    out[k+9] = a[3];    
}

static void blkdiag4(double * out, 
        const double * a, const double * b, 
        const double * c, const double * d)
{
    bzero(out, 64*sizeof(double));

    blkfill(out, a, 0);
    blkfill(out, b, 1);
    blkfill(out, c, 2);
    blkfill(out, d, 3);
}

int main(int argc, char ** argv)
{
    GPS_EKF ekf(3, 4);
    
    double T = 1; // positioning interval

    // Set Q, see [1]
    const double Sf    = 36;
    const double Sg    = 0.01;
    const double sigma = 5;         // state transition variance
    const double Qb[4] = {Sf*T+Sg*T*T*T/3, Sg*T*T/2, Sg*T*T/2, Sg*T};
    const double Qxyz[4] = {sigma*sigma*T*T*T/3, sigma*sigma*T*T/2, sigma*sigma*T*T/2, sigma*sigma*T};
    double Q[64];
    blkdiag4(Q, Qxyz, Qxyz, Qxyz, Qb);
    
    // Initial state X
    double X[8];
    
    // position
    X[0] = -2.168816181271560e+006;
    X[2] =  4.386648549091666e+006;
    X[4] =  4.077161596428751e+006;
    
    // velocity
    X[1] = 0;
    X[3] = 0;
    X[5] = 0;
    
    // clock bias
    X[6] = 3.575261153706439e+006;
    
    // clock drift
    X[7] = 4.549246345845814e+001;         
    
    dumpmat(X, 1, 8);
   
    // Initial prediction covariance
    //P = eye(8)*10;
    
    
    return 0;

    FILE * fp = fopen("gps.csv", "r");
    char line[1000];

    // skip header
    readline(line, fp);

    while (true) {

        if (!readline(line, fp))
            break;

        double SV_Pos[12];
        double SV_Rho[4];

        fill(line, SV_Pos, SV_Rho);
   }

    fclose(fp);
}
