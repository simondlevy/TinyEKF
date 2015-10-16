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

static void zeros(double * a, int n)
{
    bzero(a, n*sizeof(double));
}

class GPS_EKF : public TinyEKF {

    public:

        // Eight state values, four measurement values
        GPS_EKF(double * X, double T) : TinyEKF(8, 4) 
        {
            this->T = T;
            memcpy(this->X, X, 8*sizeof(double));
        }
        
        void setPseudorange(double * SV)
        {
            memcpy(this->SV, SV, 12*sizeof(double));
        }
                
    protected:

        void f(double * x, double * fx)
        {
            zeros(fx, 8);
            
            for (int j=0; j<8; j+=2) {
                fx[j] = x[j] + this->T * x[j+1];
                fx[j+1] = x[j+1];
            }
        }
        
        void df(double * x, double * dfx)
        {
            zeros(dfx, 64);
            
            for (int j=0; j<8; ++j)
                dfx[j*8+j] = 1;
            
            for (int j=0; j<7; ++j)
                dfx[j*8+j+1] = this->T;
        }


        void g(double * xp, double * gx, double * dgx)
        {
            dump(this->SV, 4, 3);
            exit(0);
            
            //dX = bsxfun(@minus, X([1,3,5])', SV);% X - Xs
            //Val = sum(dX .^2, 2) .^0.5 + X(7);
        }
        
    private:
        
        double X[8];   // constant velocity
        double T;      // positioning interval
        double SV[12]; // pseudorange for g function
};

static char * readline(char * line, FILE * fp)
{
    return fgets(line, 1000, fp);
}

static bool readdata(FILE * fp, double * SV_Pos, double * SV_Rho)
{
    char line[1000];
    
    if (!readline(line, fp))
        return false;
    
    char * p = strtok(line, ",");

    for (int k=0; k<12; ++k) {
        SV_Pos[k] = atof(p);
        p = strtok(NULL, ",");
    }
    
    for (int k=0; k<4; ++k) {
        SV_Rho[k] = atof(p);
        p = strtok(NULL, ",");
    }
 
    return true;
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
    zeros(out, 64);

    blkfill(out, a, 0);
    blkfill(out, b, 1);
    blkfill(out, c, 2);
    blkfill(out, d, 3);
}

static void eye(double * a, int n, double s)
{
    zeros(a, n*n);
    
    for (int k=0; k<n; ++k)
        a[k*n+k] = s;
}

static void skipline(FILE * fp)
{
    char line[1000];
    readline(line, fp);
}


int main(int argc, char ** argv)
{
    // Positioning interval
    double T = 2; 
    
    // Set Q, see [1]
    const double Sf    = 36;
    const double Sg    = 0.01;
    const double sigma = 5;         // state transition variance
    const double Qb[4] = {Sf*T+Sg*T*T*T/3, Sg*T*T/2, Sg*T*T/2, Sg*T};
    const double Qxyz[4] = {sigma*sigma*T*T*T/3, sigma*sigma*T*T/2, 
                            sigma*sigma*T*T/2, sigma*sigma*T};
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
       
    // Initial prediction covariance
    double P[64];
    eye(P, 8, 10.0);
    
    // Variance of measurement error
    double R[16];
    eye(R, 4, 36);
    
    // Open data file
    FILE * fp = fopen("gps.csv", "r");
    
    // Skip CSV header
    skipline(fp);
    
    // Inititilize EKF with constant velocity and positioning interval
    GPS_EKF ekf(X, T);

    // Loop till no more data
    while (true) {
        
        double SV_Pos[12];
        double SV_Rho[4];

        if (!readdata(fp, SV_Pos, SV_Rho))
            break;
        
        ekf.setPseudorange(SV_Pos);
        
        ekf.update(Q, R, SV_Rho, X, P);
   }

    fclose(fp);
}
