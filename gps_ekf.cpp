#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static char * readline(char * line, FILE * fp)
{
    return fgets(line, 1000, fp);
}

static void scan(char * s, double * dst, int off, int len)
{
    char * p = strtok(&s[off], ",");

    for (int k=0; k<len; ++k) {
        dst[k] = atof(p);
        printf("%2d: %f\n", k, dst[k]);
        p = strtok(NULL, ",");
    }

    exit(0);
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
        scan(line, SV_Pos, 0, 12);
        scan(line, SV_Rho, 12, 4);
    }

    fclose(fp);
}
