#include <stdio.h>

int main(int argc, char ** argv)
{
    FILE * fp = fopen("gps.csv", "r");


    while (true) {
        char line[1000];
        if (!fread(line, 1, 1000, fp))
            break;
        printf("%s", line);
    }

    fclose(fp);
}
