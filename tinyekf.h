typedef struct {

    vec_t * X;      // state
    mat_t * P;      // prediction error covariance
    mat_t * Q;      // process noise covariance
    mat_t * R;      // measurement error covariance

    mat_t  * G;    // Kalman gain; a.k.a. K

    vec_t  * Xp;   // output of state-transition function
    mat_t  * fy;   // Jacobean of process model
    mat_t  * H;    // Jacobean of measurement model
    vec_t  * gXp;

    mat_t  * Ht;
    mat_t  * fyt;
    mat_t  * Pp;

    mat_t  * eye_n_n;

    // temporary storage
    mat_t  * tmp_n_m;
    mat_t  * tmp_n_n;
    mat_t  * tmp_m_n;
    vec_t  * tmp_m;
    vec_t  * tmp_n;
    mat_t  * tmp2_n_m;
    mat_t  * tmp_m_m;
    mat_t  * tmp2_m_m;

} ekf_t; 
