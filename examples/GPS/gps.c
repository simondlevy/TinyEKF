/* gps_ekf: TinyEKF test case using You Chong's GPS example:
 * 
 *   http://www.mathworks.com/matlabcentral/fileexchange/
 *     31487-extended-kalman-filter-ekf--for-gps
 * 
 * Reads file data.csv of satellite data and writes file ekf.csv of
 * mean-subtracted estimated positions.
 *
 *
 * References:
 *
 * 1. R G Brown, P Y C Hwang, "Introduction to random signals and applied 
 * Kalman filtering : with MATLAB exercises and solutions",1996
 *
 * 2. Pratap Misra, Per Enge, "Global Positioning System Signals, 
 * Measurements, and Performance(Second Edition)",2006
 * 
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

// Size of state space
#define EKF_N 8

// Size of observation (measurement) space
#define EKF_M 4

// We need double precision to replicate the published results
#define _float_t double

#include <tinyekf.h>

// positioning interval
static const double T = 1;

// initial covariances of state noise, measurement noise
static const double P0 = 10;
static const double R0 = 36;

// Set fixed process-noise covariance matrix Q, see [1]  ---------------------

static const double Sf    = 36;
static const double Sg    = 0.01;
static const double sigma = 5;         // state transition variance

static const double b0 = Sf * T+Sg * T * T * T/3; 
static const double b1 = Sg * T * T/2; 
static const double b2 = Sg * T * T/2; 
static const double b3 = Sg * T;

static const double xyz0 = sigma * sigma * T * T * T/3; 
static const double xyz1 = sigma * sigma * T * T/2; 
static const double xyz2 = sigma * sigma * T * T/2; 
static const double xyz3 = sigma * sigma * T;

static const double Q[8*8] = {

    xyz0, xyz1,  0,     0,     0,     0,     0,   0,
    xyz2, xyz3,  0,     0,     0,     0,     0,   0,
    0,     0,      xyz0, xyz1, 0,     0,     0,   0,
    0,     0,      xyz2, xyz3, 0,     0,     0,   0,
    0,     0,      0,     0,     xyz0, xyz1, 0,   0,
    0,     0,      0,     0,     xyz2, xyz3, 0,   0,
    0,     0,      0,     0,     0,     0,     b0, b1,
    0,     0,      0,     0,     0,     0,     b2, b3
};

// Set fixed measurement noise covariance matrix R ----------------------------

static const double R[4*4] = {
    R0, 0, 0, 0,
    0, R0, 0, 0,
    0, 0, R0, 0,
    0, 0, 0, R0

};

static void init(ekf_t * ekf)
{
    const double pdiag[8] = { P0, P0, P0, P0, P0, P0, P0, P0 };

    ekf_initialize(ekf, pdiag);

    // position
    ekf->x[0] = -2.168816181271560e+006;
    ekf->x[2] =  4.386648549091666e+006;
    ekf->x[4] =  4.077161596428751e+006;

    // velocity
    ekf->x[1] = 0;
    ekf->x[3] = 0;
    ekf->x[5] = 0;

    // clock bias
    ekf->x[6] = 3.575261153706439e+006;

    // clock drift
    ekf->x[7] = 4.549246345845814e+001;
}

static void run_model(
        ekf_t * ekf, 
        const double SV[4][3], 
        double fx[8],
        double F[8*8],
        double hx[4],
        double H[4*8])
{ 

    for (int j=0; j<8; j+=2) {
        fx[j] = ekf->x[j] + T * ekf->x[j+1];
        fx[j+1] = ekf->x[j+1];
    }

    for (int j=0; j<8; ++j) {
        F[j*8+j] = 1;
    }

    for (int j=0; j<4; ++j) {
        F[2*j*8+2*j+1] = T;
    }

    double dx[4][3];

    for (int i=0; i<4; ++i) {
        hx[i] = 0;
        for (int j=0; j<3; ++j) {
            double d = fx[j*2] - SV[i][j];
            dx[i][j] = d;
            hx[i] += d*d;
        }
        hx[i] = pow(hx[i], 0.5) + fx[6];
    }

    for (int i=0; i<4; ++i) {

        for (int j=0; j<3; ++j) {

            H[i*8+j*2] = dx[i][j] / hx[i];
        }

        H[i*8+6] = 1;
    }   
}

static void readline(char * line, FILE * fp)
{
    fgets(line, 1000, fp);
}

static void readdata(FILE * fp, double SV_Pos[4][3], double SV_Rho[4])
{
    char line[1000];

    readline(line, fp);

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
}


static void skipline(FILE * fp)
{
    char line[1000];
    readline(line, fp);
}

void error(const char * msg)
{
    fprintf(stderr, "%s\n", msg);
}

int main(int argc, char ** argv)
{    
    ekf_t ekf = {0};

    init(&ekf);

    // Open input data file
    FILE * ifp = fopen("data.csv", "r");

    // Skip CSV header
    skipline(ifp);

    // Make a place to store the data from the file and the output of the EKF
    double SV_Pos[4][3] = {0};
    double SV_Rho[4] = {0};
    double Pos_KF[25][3] = {0};

    // Open output CSV file and write header
    const char * OUTFILE = "ekf.csv";
    FILE * ofp = fopen(OUTFILE, "w");
    fprintf(ofp, "X,Y,Z\n");

    // Loop till no more data
    for (int j=0; j<25; ++j) {

        readdata(ifp, SV_Pos, SV_Rho);

        // -------------------------------------------------------------------

        // Run our model to get the EKF inputs
        double fx[8] = {0};
        double F[8*8] = {0};
        double hx[4] = {0};
        double H[4*8] = {0};
        run_model(&ekf, SV_Pos, fx, F, hx, H);

        // Run the EKF prediction step
        ekf_predict(&ekf, fx, F, Q);

        // Run the EKF update step
        ekf_update(&ekf, SV_Rho, hx, H, R);

        // grab positions, ignoring velocities
        for (int k=0; k<3; ++k) {
            Pos_KF[j][k] = ekf.x[2*k];
        }
    }

    // Compute means of filtered positions
    double mean_Pos_KF[3] = {0, 0, 0};
    for (int j=0; j<25; ++j) 
        for (int k=0; k<3; ++k)
            mean_Pos_KF[k] += Pos_KF[j][k];
    for (int k=0; k<3; ++k)
        mean_Pos_KF[k] /= 25;


    // Dump filtered positions minus their means
    for (int j=0; j<25; ++j) {
        fprintf(ofp, "%f,%f,%f\n", 
                Pos_KF[j][0]-mean_Pos_KF[0], 
                Pos_KF[j][1]-mean_Pos_KF[1], 
                Pos_KF[j][2]-mean_Pos_KF[2]);
        printf("%f %f %f\n", Pos_KF[j][0], Pos_KF[j][1], Pos_KF[j][2]);
    }

    // Done!
    fclose(ifp);
    fclose(ofp);
    printf("Wrote file %s\n", OUTFILE);
    return 0;
}
