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
#include <math.h>

#include "tinyekf.hpp"

class GPS_EKF : public TinyEKF {

    public:

        // Eight state values, four measurement values
        GPS_EKF(double * X, double T, double P0, double R0) : TinyEKF(8, 4) 
        {
            this->T = T;
            copy(this->X, X, 8);
            this->SV = newmat(4, 3);
            
            eye(this->P, 8, P0);
            eye(this->R, 4, R0);
        }
        
        void setPseudorange(double ** SV)
        {
            copy(this->SV, SV, 4, 3);
        }
                
    protected:
    
        void f(double * X, double * Xp, double ** fy)
        {
            zeros(Xp, 8);
            
            for (int j=0; j<8; j+=2) {
                Xp[j] = X[j] + this->T * X[j+1];
                Xp[j+1] = X[j+1];
            }
            
            zeros(fy, 8, 8);
            
            for (int j=0; j<8; ++j)
                fy[j][j] = 1;
                
            for (int j=0; j<4; ++j)
                fy[2*j][2*j+1] = this->T;
        }


        void g(double * Xp, double * gXp, double ** H)
        {
            double dx[12];
            
            for (int i=0; i<4; ++i) {
                gXp[i] = 0;
                for (int j=0; j<3; ++j) {
                    double d = this->X[j*2] - this->SV[i][j];
                    dx[i*3+j] = d;
                    gXp[i] += d*d;
                }
                gXp[i] = sqrt(gXp[i]) + this->X[6];
            }
            
            zeros(H, 4, 8);
            for (int i=0; i<4; ++i) {
                for (int j=0; j<3; ++j) {
                    H[i][j*2] = dx[i*3+j] / gXp[i];
                }
                H[i][6] = 1;
            }   
        }
        
    private:
        
        double    X[8]; // constant velocity
        double    T;    // positioning interval
        double ** SV;   // pseudorange for g function
};

static char * readline(char * line, FILE * fp)
{
    return fgets(line, 1000, fp);
}

static bool readdata(FILE * fp, double ** SV_Pos, double * SV_Rho)
{
    char line[1000];
    
    if (!readline(line, fp))
        return false;
    
    char * p = strtok(line, ",");

    for (int i=0; i<4; ++i)
        for (int j=0; j<3; ++j) {
            SV_Pos[i][j] = atof(p);
            p = strtok(NULL, ",");
        }
    
    for (int j=0; j<4; ++j) {
        SV_Rho[j] = atof(p);
        p = strtok(NULL, ",");
    }
 
    return true;
}

static void blkfill(double ** out, const double * a, int off)
{   
    off *= 2;
    
    out[off][off]     = a[0];
    out[off][off+1]   = a[1];
    out[off+1][off]   = a[2];
    out[off+1][off+1] = a[3];    
}


static void skipline(FILE * fp)
{
    char line[1000];
    readline(line, fp);
}


int main(int argc, char ** argv)
{
    // Positioning interval
    double T = 1; 
    
    // Set Q, see [1]
    const double Sf    = 36;
    const double Sg    = 0.01;
    const double sigma = 5;         // state transition variance
    const double Qb[4] = {Sf*T+Sg*T*T*T/3, Sg*T*T/2, Sg*T*T/2, Sg*T};
    const double Qxyz[4] = {sigma*sigma*T*T*T/3, sigma*sigma*T*T/2, 
                            sigma*sigma*T*T/2, sigma*sigma*T};
                            
    double ** Q = TinyEKF::newmat(8,8);
    blkfill(Q, Qxyz, 0);
    blkfill(Q, Qxyz, 1);
    blkfill(Q, Qxyz, 2);
    blkfill(Q, Qb,   3);    
    
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
       
    // Inititilize EKF with constant velocity and positioning interval
    // initial prediction covariance 10, initial measurement error 36
    GPS_EKF ekf(X, T, 10, 36);
    
    // Open data file
    FILE * fp = fopen("gps.csv", "r");
    
    // Skip CSV header
    skipline(fp);
        
    double ** SV_Pos = TinyEKF::newmat(4,3);
    double SV_Rho[4];
    
    // Loop till no more data
    while (true) {
        
        if (!readdata(fp, SV_Pos, SV_Rho))
            break;
                
        ekf.setPseudorange(SV_Pos);
        
        ekf.update(SV_Rho, X);
   }

    fclose(fp);
}
