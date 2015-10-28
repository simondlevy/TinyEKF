/* ------------------------- */
#define N 8
#define M 4

typedef double number_t;
/* ------------------------- */

typedef struct {

    number_t x[N];      // state

    number_t P[N][N];   // prediction error covariance
    number_t Q[N][N];   // process noise covariance
    number_t R[M][M];   // measurement error covariance

    number_t G[N][M];   // Kalman gain; a.k.a. K
    number_t F[N][N];   // Jacobian of process model
    number_t H[M][N];   // Jacobian of measurement model

    number_t Ht[N][M];  // transpose of measurement Jacobian
    number_t Ft[N][N];  // transpose of process Jacobian
    number_t Pp[N][N];  // P, post-prediction, pre-update

    number_t fx[N];     // f(x)
    number_t hx[N];     // h(x)

    // temporary storage
    number_t tmp1[N*N];
    number_t tmp2[M*N];
    number_t tmp3[M*M];
    number_t tmp4[M*M];
    number_t tmp5[M];

} ekf_t; 

void ekf_init(ekf_t * ekf);

void ekf_step(
        ekf_t * ekf, 
        number_t * Z, 
        void (*f)(number_t x[N],  number_t fx[N], number_t F[N][N]), 
        void (*h)(number_t fx[N], number_t hx[N], number_t H[M][N]));

void ekf_predict_and_update(ekf_t * ekf, number_t * Z);
