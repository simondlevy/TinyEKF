/* have to define N above here */

static void choldc1(double a[N][N], double p[N]) {
    int i,j,k;
    double sum;

    for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
            sum = a[i][j];
            for (k = i - 1; k >= 0; k--) {
                sum -= a[i][k] * a[j][k];
            }
            if (i == j) {
                if (sum <= 0) {
                    printf(" a is not positive definite!\n");
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j][i] = sum / p[i];
            }
        }
    }
}

static void choldcsl(double A[N][N], double a[N][N], double p[N]) 
{
    int i,j,k; double sum;
    for (i = 0; i < N; i++) 
        for (j = 0; j < N; j++) 
            a[i][j] = A[i][j];
    choldc1(a, p);
    for (i = 0; i < N; i++) {
        a[i][i] = 1 / p[i];
        for (j = i + 1; j < N; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j][k] * a[k][i];
            }
            a[j][i] = sum / p[j];
        }
    }
}


static void cholsl(double A[N][N], double a[N][N]) 
{
    int i,j,k;
    double p[N];
    choldcsl(A,a,p);
    for (i = 0; i < N; i++) {
        for (j = i + 1; j < N; j++) {
            a[i][j] = 0.0;
        }
    }
    for (i = 0; i < N; i++) {
        a[i][i] *= a[i][i];
        for (k = i + 1; k < N; k++) {
            a[i][i] += a[k][i] * a[k][i];
        }
        for (j = i + 1; j < N; j++) {
            for (k = j; k < N; k++) {
                a[i][j] += a[k][i] * a[k][j];
            }
        }
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < i; j++) {
            a[i][j] = a[j][i];
        }
    }
}

