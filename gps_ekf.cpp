#include <stdio.h>

int main(int argc, char ** argv)
{
    FILE * fp = fopen("gps.csv", "r");
    char line[1000];

    // skip header
    fgets(line, 1000, fp);

    while (true) {
        if (!fgets(line, 1000, fp))
            break;
        printf("%s", line);
    }

    fclose(fp);
}
