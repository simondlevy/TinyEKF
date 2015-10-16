#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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
            printf("%6.6f ", a[j*n+k]);
        printf("\n");
    }
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
    dumpmat((double *)Qxyz, 2, 2);
    //double Q[64];
    //blkdiag(Q, Qxyz,Qxyz,Qxyz,Qb);
    
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
