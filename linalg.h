void mat_init(double * * mat, int m, int n);

void mat_free(double * mat);

void zeros(double * mat, int m, int n);

void eye(double * a, int n, double s);

void vec_dump(double * x, int n, const char * fmt);

void mat_dump(double * mat, int m, int n, const char * fmt);

// C <- A * B
void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols);

// Y <- A * X
void mulvec(double * a, double * x, double * y, int m, int n);

void transpose(double * a, double * at, int m, int n);

// A <- A + B
void add(double * a, double * b, int m, int n);

// C <- A - B
void sub(double * a, double * b, double * c, int n);

void negate(double * a, int m, int n);

void invert(double * a, double * at, double * p, int n);

void mat_set(double * a, int i, int j, int n, double value);

void mat_addeye(double * a, int n);
