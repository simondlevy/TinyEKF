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
        GPS_EKF() : TinyEKF(8, 4)
        {            
            // positioning interval
            this->T = 1; 

           // Set Q, see [1]
            const double Sf    = 36;
            const double Sg    = 0.01;
            const double sigma = 5;         // state transition variance
            const double Qb[4] = {Sf*T+Sg*T*T*T/3, Sg*T*T/2, Sg*T*T/2, Sg*T};
            const double Qxyz[4] = {sigma*sigma*T*T*T/3, sigma*sigma*T*T/2,
                sigma*sigma*T*T/2, sigma*sigma*T};

            this->blkfill(Qxyz, 0);
            this->blkfill(Qxyz, 1);
            this->blkfill(Qxyz, 2);
            this->blkfill(Qb,   3);

            // initial covariances of state, measurement noise 
            double P0 = 10;
            double R0 = 36;

            for (int i=0; i<8; ++i)
                this->setP(i, i, P0);

            for (int i=0; i<4; ++i)
                this->setR(i, i, R0);

            // position
            this->x[0] = -2.168816181271560e+006;
            this->x[2] =  4.386648549091666e+006;
            this->x[4] =  4.077161596428751e+006;

            // velocity
            this->x[1] = 0;
            this->x[3] = 0;
            this->x[5] = 0;

            // clock bias
            this->x[6] = 3.575261153706439e+006;

            // clock drift
            this->x[7] = 4.549246345845814e+001;
            }

        void setPseudorange(double  SV[4][3])
        {
            for (int i=0; i<4; ++i)
                for (int j=0; j<3; ++j)
                    this->SV[i][j] = SV[i][j];
        }

    protected:

        void f(double * x, double * fx, double * F)
        {
            for (int j=0; j<8; j+=2) {
                fx[j] = x[j] + this->T * x[j+1];
                fx[j+1] = x[j+1];
            }

            for (int j=0; j<8; ++j)
                this->set(F, j, j, 1);

            for (int j=0; j<4; ++j)
                this->set(F, 2*j, 2*j+1, this->T);
        }

        void h(double * fx, double * hx, double * H)
        {
            double dx[4][3];

            for (int i=0; i<4; ++i) {
                hx[i] = 0;
                for (int j=0; j<3; ++j) {
                    double d = fx[j*2] - this->SV[i][j];
                    dx[i][j] = d;
                    hx[i] += d*d;
                }
                hx[i] = pow(hx[i], 0.5) + fx[6];
            }

            for (int i=0; i<4; ++i) {
                for (int j=0; j<3; ++j) 
                    this->set(H, i, j*2, dx[i][j] / hx[i]);
                this->set(H, i, 6, 1);
            }   
        }

    private:

        void blkfill(const double * a, int off)
        {
            off *= 2;

            this->setQ(off, off,     a[0]); 
            this->setQ(off, off+1,   a[1]);
            this->setQ(off+1, off,   a[2]);
            this->setQ(off+1, off+1, a[3]);
        }

        double  T;          // positioning interval
        double  SV[4][3];   // pseudorange for g function
};

static char * readline(char * line, FILE * fp)
{
    return fgets(line, 1000, fp);
}

static bool readdata(FILE * fp, double SV_Pos[4][3], double SV_Rho[4])
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


static void skipline(FILE * fp)
{
    char line[1000];
    readline(line, fp);
}

int main(int argc, char ** argv)
{    
    // Create the EKF
    GPS_EKF ekf;

    // Open data file
    FILE * fp = fopen("gps.csv", "r");

    // Skip CSV header
    skipline(fp);

    // Make a place to store the data from the file
    double SV_Pos[4][3];
    double SV_Rho[4];

    // Loop till no more data
    while (true) {

        if (!readdata(fp, SV_Pos, SV_Rho))
            break;

        ekf.setPseudorange(SV_Pos);

        ekf.step(SV_Rho);

        printf("%f %f %f\n", ekf.getX(0), ekf.getX(2), ekf.getX(4));
    }

    fclose(fp);
}
