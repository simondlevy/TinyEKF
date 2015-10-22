typedef struct {

    double * data;
    int len;

} vec_t;

typedef struct {

    double ** data;
    int rows;
    int cols;

    double * tmp;

} mat_t;

void vec_init(vec_t * vec, int n);

void vec_free(vec_t vec);

void mat_init(mat_t * mat, int m, int n);

void mat_free(mat_t mat);

void zeros(mat_t mat);

void eye(mat_t mat, double s);

void dumpvec(vec_t vec, const char * fmt);

void dumpmat(mat_t mat, const char * fmt);

// C <- A * B
void mulmat(mat_t a, mat_t b, mat_t c);

void mulvec(mat_t a, vec_t x, vec_t y);

void transpose(mat_t a, mat_t at);

// A <- A + B
void add(mat_t a, mat_t b);

// A <- A - B
void sub(vec_t a, vec_t b);

void negate(mat_t a);

void invert(mat_t a, mat_t at);
