#include <stdio.h>

int main(int argc, char ** argv)
{
    FILE * fp = fopen("gps.csv", "r");

    fclose(fp);
}
