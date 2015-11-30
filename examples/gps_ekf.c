/* gps_ekf: TinyEKF test case using You Chong's GPS example:
 * 
 *   http://www.mathworks.com/matlabcentral/fileexchange/31487-extended-kalman-filter-ekf--for-gps
 * 
 * Reads file gps.csv of satellite data and writes file ekf.csv of mean-subtracted estimated positions.
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

#include "tinyekf_config.h"
#include "tiny_ekf.h"

// positioning interval
static const double T = 1;

static void blkfill(ekf_t * ekf, const double * a, int off)
{
    off *= 2;

    ekf->Q[off]   [off]   = a[0]; 
    ekf->Q[off]   [off+1] = a[1];
    ekf->Q[off+1] [off]   = a[2];
    ekf->Q[off+1] [off+1] = a[3];
}


static void init(ekf_t * ekf)
{
    // Set Q, see [1]
    const double Sf    = 36;
    const double Sg    = 0.01;
    const double sigma = 5;         // state transition variance
    const double Qb[4] = {Sf*T+Sg*T*T*T/3, Sg*T*T/2, Sg*T*T/2, Sg*T};
    const double Qxyz[4] = {sigma*sigma*T*T*T/3, sigma*sigma*T*T/2, sigma*sigma*T*T/2, sigma*sigma*T};

    blkfill(ekf, Qxyz, 0);
    blkfill(ekf, Qxyz, 1);
    blkfill(ekf, Qxyz, 2);
    blkfill(ekf, Qb,   3);

    // initial covariances of state noise, measurement noise
    double P0 = 10;
    double R0 = 36;

    int i;

    for (i=0; i<8; ++i)
        ekf->P[i][i] = P0;

    for (i=0; i<4; ++i)
        ekf->R[i][i] = R0;

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

static void model(ekf_t * ekf, double SV[4][3])
{ 

    int i, j;

    for (j=0; j<8; j+=2) {
        ekf->fx[j] = ekf->x[j] + T * ekf->x[j+1];
        ekf->fx[j+1] = ekf->x[j+1];
    }

    for (j=0; j<8; ++j)
        ekf->F[j][j] = 1;

    for (j=0; j<4; ++j)
        ekf->F[2*j][2*j+1] = T;

    double dx[4][3];

    for (i=0; i<4; ++i) {
        ekf->hx[i] = 0;
        for (j=0; j<3; ++j) {
            double d = ekf->fx[j*2] - SV[i][j];
            dx[i][j] = d;
            ekf->hx[i] += d*d;
        }
        ekf->hx[i] = pow(ekf->hx[i], 0.5) + ekf->fx[6];
    }

    for (i=0; i<4; ++i) {
        for (j=0; j<3; ++j) 
            ekf->H[i][j*2]  = dx[i][j] / ekf->hx[i];
        ekf->H[i][6] = 1;
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

    int i, j;

    for (i=0; i<4; ++i)
        for (j=0; j<3; ++j) {
            SV_Pos[i][j] = atof(p);
            p = strtok(NULL, ",");
        }

    for (j=0; j<4; ++j) {
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
    // Do generic EKF initialization
    ekf_t ekf;
    ekf_init(&ekf, N, M);

    // Do local initialization
    init(&ekf);

    // Open input data file
    FILE * ifp = fopen("gps.csv", "r");

    // Skip CSV header
    skipline(ifp);

    // Make a place to store the data from the file and the output of the EKF
    double SV_Pos[4][3];
    double SV_Rho[4];
    double Pos_KF[25][3];

    // Open output CSV file and write header
    const char * OUTFILE = "ekf.csv";
    FILE * ofp = fopen(OUTFILE, "w");
    fprintf(ofp, "X,Y,Z\n");

    int j, k;

    // Loop till no more data
    for (j=0; j<25; ++j) {

        readdata(ifp, SV_Pos, SV_Rho);

        model(&ekf, SV_Pos);

        ekf_step(&ekf, SV_Rho);

        // grab positions, ignoring velocities
        for (k=0; k<3; ++k)
            Pos_KF[j][k] = ekf.x[2*k];
    }

    // Compute means of filtered positions
    double mean_Pos_KF[3] = {0, 0, 0};
    for (j=0; j<25; ++j) 
        for (k=0; k<3; ++k)
            mean_Pos_KF[k] += Pos_KF[j][k];
    for (k=0; k<3; ++k)
        mean_Pos_KF[k] /= 25;


    // Dump filtered positions minus their means
    for (j=0; j<25; ++j) {
        fprintf(ofp, "%f,%f,%f\n", 
                Pos_KF[j][0]-mean_Pos_KF[0], Pos_KF[j][1]-mean_Pos_KF[1], Pos_KF[j][2]-mean_Pos_KF[2]);
        printf("%f %f %f\n", Pos_KF[j][0], Pos_KF[j][1], Pos_KF[j][2]);
    }
    
    // Done!
    fclose(ifp);
    fclose(ofp);
    printf("Wrote file %s\n", OUTFILE);
    return 0;
}
