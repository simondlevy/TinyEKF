static void choldc1(double ** a, double * p, int n) {
    int i,j,k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
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

static void choldcsl(double ** A, double ** a, double * p, int n) 
{
    int i,j,k; double sum;
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            a[i][j] = A[i][j];
    choldc1(a, p, n);
    for (i = 0; i < n; i++) {
        a[i][i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j][k] * a[k][i];
            }
            a[j][i] = sum / p[j];
        }
    }
}


static void invert(double ** A, double ** a, double * p, int n) 
{
    int i,j,k;
    choldcsl(A,a,p,n);
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i][j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i][i] *= a[i][i];
        for (k = i + 1; k < n; k++) {
            a[i][i] += a[k][i] * a[k][i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i][j] += a[k][i] * a[k][j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i][j] = a[j][i];
        }
    }
}

static double ** newmat(int m, int n)
{
    double ** a = new double * [m];

    for (int i=0; i<m; ++i)
        a[i] = new double [n];

    return a;
}

static double deletemat(double ** a, int m)
{
    for (int i=0; i<m; ++i)
        delete a[i];
}

static void zeros(double * x, int n)
{
    bzero(x, n*sizeof(double));
}

static void zeros(double ** a, int m, int n)
{
    for (int i=0; i<m; ++i)
        zeros(a[i], n);
}


static void eye(double ** a, int n, double s)
{
    zeros(a, n, n);

    for (int k=0; k<n; ++k)
        a[k][k] = s;
}

static void dump(double * x, int n, const char * fmt)
{
    char f[100];
    sprintf(f, "%s ", fmt);
    for (int j=0; j<n; ++j)
        printf(f, x[j]);
    printf("\n");
}

static void dump(double ** a, int m, int n, const char * fmt)
{
    for (int i=0; i<m; ++i) {
        dump(a[i], n, fmt);
    }
}

// C <- A * B
static void mul(double ** a, double ** b, double ** c, int m, int n, int p)
{
    for (int i=0; i<m; ++i)
        for (int j=0; j<n; ++j) {
            c[i][j] = 0;
            for (int l=0; l<p; ++l)
                c[i][j] += a[i][l] * b[l][j];
        }
}

static void mul(double ** a, double * x, double * y, int rows, int cols)
{
    for (int i=0; i<rows; ++i) {
        y[i] = 0;
        for (int j=0; j<cols; ++j)
            y[i] += x[j] * a[i][j];
    }
}

static void transpose(double ** a, double ** at, int rows, int cols)
{
    for (int i=0; i<rows; ++i)
        for (int j=0; j<cols; ++j) 
            at[j][i] = a[i][j];
}

// X <- X + Y
static void add(double * x, double * y, int n)
{        
   for (int j=0; j<n; ++j)
       x[j] += y[j];
}

// A <- A + B
static void add(double ** a, double ** b, int rows, int cols) 
{        
    for (int i=0; i<rows; ++i)
        for (int j=0; j<cols; ++j)
            a[i][j] += b[i][j];
}

// C <- A - B

static void sub(double * a, double * b, double * c, int n)
{
    for (int j=0; j<n; ++j)
        c[j] = a[j] - b[j];
}

static void sub(double ** a, double ** b, double ** c, int m, int n)
{
    for (int i=0; i<m; ++i)
        for (int j=0; j<n; ++j)
            c[i][j] = a[i][j] -  b[i][j];
}

typedef struct {

    double * data;
    int len;

} vec_t;

typedef struct {

    double ** data;
    int rows;
    int cols;
} mat_t;

static mat_t * newnewmat(int m, int n)
{
    mat_t * mat = new mat_t;

    mat->data = new double * [m];

    for (int i=0; i<m; ++i)
        mat->data[i] = new double [n];

    mat->rows = m;
    mat->cols = n;

    return mat;
}

static void deletemat(mat_t * mat)
{
    for (int i=0; i<mat->rows; ++i)
        delete mat->data[i];

    delete mat->data;

    delete mat;
}

static void zeros(vec_t * vec)
{
    bzero(vec->data, vec->len*sizeof(double));
}

static void zeros(mat_t * mat)
{
    for (int i=0; i<mat->rows; ++i)
        bzero(mat->data[i], mat->cols*sizeof(double));
}

static void eye(mat_t * mat, double s)
{
    zeros(mat);

    for (int k=0; k<mat->rows; ++k)
        mat->data[k][k] = s;
}

static void dump(vec_t * vec, const char * fmt)
{
    char f[100];
    sprintf(f, "%s ", fmt);
    for (int j=0; j<vec->len; ++j)
        printf(f, vec->data[j]);
    printf("\n");
}

static void dump(mat_t * mat, const char * fmt)
{
    char f[100];
    sprintf(f, "%s ", fmt);
     for (int i=0; i<mat->rows; ++i) {
        for (int j=0; j<mat->cols; ++j)
            printf(fmt, mat->data[i][j]);
        printf("\n");
    }
}

// C <- A * B
static void mul(mat_t * a, mat_t * b, mat_t * c)
{
    for (int i=0; i<a->rows; ++i)
        for (int j=0; j<b->cols; ++j) {
            c->data[i][j] = 0;
            for (int l=0; l<a->cols; ++l)
                c->data[i][j] += a->data[i][l] * b->data[l][j];
        }
}

static void mul(mat_t * a, vec_t * x, vec_t * y)
{
    for (int i=0; i<a->rows; ++i) {
        y->data[i] = 0;
        for (int j=0; j<a->cols; ++j)
            y->data[i] += x->data[j] * a->data[i][j];
    }
}

static void transpose(mat_t * a, mat_t * at)
{
    for (int i=0; i<a->rows; ++i)
        for (int j=0; j<a->cols; ++j) 
            at->data[j][i] = a->data[i][j];
}

// X <- X + Y
static void add(vec_t * x, vec_t * y)
{        
   for (int j=0; j<x->len; ++j)
       x->data[j] += y->data[j];
}

// A <- A + B
static void add(mat_t * a, mat_t * b)
{        
    for (int i=0; i<a->rows; ++i)
        for (int j=0; j<a->cols; ++j)
            a->data[i][j] += b->data[i][j];
}

// C <- A - B

static void sub(vec_t * a, vec_t * b, vec_t * c)
{
    for (int j=0; j<a->len; ++j)
        c->data[j] = a->data[j] - b->data[j];
}

static void sub(mat_t * a, mat_t * b, mat_t * c)
{
    for (int i=0; i<a->rows; ++i)
        for (int j=0; j<a->cols; ++j)
            c->data[i][j] = a->data[i][j] -  b->data[i][j];
}


