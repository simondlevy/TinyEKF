/**
 * Custom EKF methods
 *
 * Copyright (C) 2024 Simon D. Levy
 *
 * MIT License
 */

/// @private
static void outer(
        const _float_t x[EKF_N],
        const _float_t y[EKF_N],
        _float_t a[EKF_N*EKF_N]) 
{
    for (int i=0; i<EKF_N; i++) {
        for (int j=0; j<EKF_N; j++) {
            a[i*EKF_N+j] = x[i] * y[j];
        }
    }
}

/// @private
static _float_t dot(const _float_t x[EKF_N], const _float_t y[EKF_N]) 
{
    _float_t d = 0;

    for (int k=0; k<EKF_N; k++) {
        d += x[k] * y[k];
    }

    return d;
}

/**
  * Runs a custom update on the covariance matrix
  * @param ekf pointer to an ekf_t structure
  * @param A from the update P <- A P A^T
  */
static void ekf_custom_multiply_covariance(
        ekf_t * ekf, const _float_t A[EKF_N*EKF_N]) 
{
    _float_t AP[EKF_N*EKF_N] = {};
    _mulmat(A, ekf->P,  AP, EKF_N, EKF_N, EKF_N);

    _float_t At[EKF_N*EKF_N] = {};
    _transpose(A, At, EKF_N, EKF_N);

    _mulmat(AP, At, ekf->P, EKF_N, EKF_N, EKF_N);
}

/**
  * Enforces symmetry of the covariance matrix and ensures that the its values stay bounded
  * @param ekf pointer to an ekf_t structure
  * @param minval minimum covariance bound
  * @param maxval maximum covariance bound
  * 
  */
static void ekf_custom_cleanup_covariance(
        ekf_t * ekf, const float minval, const float maxval)
{

    for (int i=0; i<EKF_N; i++) {

        for (int j=i; j<EKF_N; j++) {

            const _float_t pval = (ekf->P[i*EKF_N+j] + ekf->P[EKF_N*j+i]) / 2;

            ekf->P[i*EKF_N+j] = ekf->P[j*EKF_N+i] =
                pval > maxval ?  maxval :
                (i==j && pval < minval) ?  minval :
                pval;
        }
    }
}

/**
  * Updates the EKF with a single scalar observation
  * @param ekf pointer to an ekf_t structure
  * @param z the observation
  * @param hx the predicted value
  * @param h one column of the sensor-function Jacobian matrix H
  * @param r one entry in the measurement-noise matrix R
  * 
  */
static void ekf_custom_scalar_update(
        ekf_t * ekf,
        const _float_t z,
        const _float_t hx,
        const _float_t h[EKF_N], 
        const _float_t r)
{
    (void)ekf_update;

    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    _float_t ph[EKF_N] = {};
    _mulvec(ekf->P, h, ph, EKF_N, EKF_N);
    const _float_t hphtr_inv = 1 / (r + dot(h, ph)); 
    _float_t g[EKF_N] = {};
    for (int i=0; i<EKF_N; ++i) {
        g[i] = ph[i] * hphtr_inv;
    }

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
    for (int i=0; i<EKF_N; ++i) {
        ekf->x[i] += g[i] * (z - hx);
    }

    // P_k = (I - G_k H_k) P_k$
    _float_t GH[EKF_N*EKF_N];
    outer(g, h, GH); 
    ekf_update_step3(ekf, GH);

    // Does this belong here, or in caller?
    for (int i=0; i<EKF_N; i++) {
        for (int j=i; j<EKF_N; j++) {
            ekf->P[i*EKF_N+j] += r * g[i] * g[j];
        }
    }
}
