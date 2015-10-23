#define N 8
#define M 4

typedef struct {

    double X[N];     // state

    double P[N*N];   // prediction error covariance
    double Q[N*N];   // process noise covariance
    double R[M*M];   // measurement error covariance

    double  G[N*M];  // Kalman gain; a.k.a. K

    double  Xp[N];   // output of state-transition function
    double  fy[N*N]; // Jacobean of process model
    double  H[M*N];  // Jacobean of measurement model
    double  gXp[N];

    double  Ht[N*M];
    double  fyt[N*N];
    double  Pp[N*N];

    // temporary storage
    double  tmp_n_m[N*M];
    double  tmp_n_n[N*N];
    double  tmp_m_n[M*N];
    double  tmp_m[M];
    double  tmp2_n_m[N*M];
    double  tmp_m_m[M*M];
    double  tmp2_m_m[M*M];

} ekf_t; 

void ekf_setP(ekf_t * ekf, int i, int j, double value);

void ekf_setQ(ekf_t * ekf, int i, int j, double value);

void ekf_setR(ekf_t * ekf, int i, int j, double value);

void ekf_setX(ekf_t * ekf, int i, double value);

void ekf_update(
        ekf_t * ekf, 
        double * Z, 
        void (*f)(double *, double *, double *), 
        void (*g)(double *, double *, double *));

void ekf_post_update(ekf_t * ekf, double * Z);
