#include <stdio.h>

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
        printf("%s", line);
    }

    fclose(fp);
}
