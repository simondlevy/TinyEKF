typedef struct {

    double ** data;
    int rows;
    int cols;

    double * tmp;

} mat_t;

void mat_init(mat_t * mat, int m, int n);

void mat_free(mat_t mat);

void zeros(mat_t mat);

void eye(mat_t mat, double s);

void vec_dump(double * x, int n, const char * fmt);

void mat_dump(mat_t mat, const char * fmt);

// C <- A * B
void mulmat(mat_t a, mat_t b, mat_t c);

void mulvec(mat_t a, double * x, double * y);

void transpose(mat_t a, mat_t at);

// A <- A + B
void add(mat_t a, mat_t b);

// C <- A - B
void sub(double * a, double * b, double * c, int n);

void negate(mat_t a);

void invert(mat_t a, mat_t at);

void mat_set(mat_t a, int i, int j, double value);
