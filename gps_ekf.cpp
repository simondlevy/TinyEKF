#include <stdio.h>

int main(int argc, char ** argv)
{
    FILE * fp = fopen("gps.csv", "r");
    char line[1000];

    // skip header
    fread(line, 1, 1000, fp);

    while (true) {
        if (!fread(line, 1, 1000, fp))
            break;
        printf("%s", line);
    }

    fclose(fp);
}
