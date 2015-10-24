/* ------------------------- */
#define N 8
#define M 4

typedef double number_t;
/* ------------------------- */

typedef struct {

    number_t X[N];     // state

    number_t P[N][N];   // prediction error covariance
    number_t Q[N][N];   // process noise covariance
    number_t R[M][M];   // measurement error covariance

    number_t  G[N][M];  // Kalman gain; a.k.a. K

    number_t  Xp[N];   // output of state-transition function
    number_t  fy[N][N]; // Jacobean of process model
    number_t  H[M][N];  // Jacobean of measurement model
    number_t  gXp[N];

    number_t  Ht[N][M];
    number_t  fyt[N][N];
    number_t  Pp[N][N];

    // temporary storage
    number_t  tmp_n_m[N][M];
    number_t  tmp_n_n[N][N];
    number_t  tmp_m_n[M][N];
    number_t  tmp_m[M];
    number_t  tmp2_n_m[N][M];
    number_t  tmp_m_m[M][M];
    number_t  tmp2_m_m[M][M];

} ekf_t; 

void ekf_update(
        ekf_t * ekf, 
        number_t * Z, 
        void (*f)(number_t X[N], number_t Xp[N], number_t fy[N][N]), 
        void (*g)(number_t Xp[N], number_t gXp[N], number_t H[M][N]));

void ekf_post_update(ekf_t * ekf, number_t * Z);
