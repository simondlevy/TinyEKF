#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "tinyekf.hpp"

static char * readline(char * line, FILE * fp)
{
    return fgets(line, 1000, fp);
}

int main(int argc, char ** argv)
{
    FILE * fp = fopen("gps.csv", "r");
    char line[1000];

    // skip header
    readline(line, fp);

    while (true) {

        if (!readline(line, fp))
            break;

        double SV_Pos[12];
        double SV_Rho[4];

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

    fclose(fp);
}
