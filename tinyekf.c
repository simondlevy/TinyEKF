#include<stdlib.h>

#include "tinyekf.h"

void ekf_setP(ekf_t * ekf, int i, int j, double value)
{
    mat_set(ekf->P, i, j, N, value);
}

void ekf_setQ(ekf_t * ekf, int i, int j, double value)
{
    mat_set(ekf->Q, i, j, N, value);
}

void ekf_setR(ekf_t * ekf, int i, int j, double value)
{
    mat_set(ekf->R, i, j, M, value);
}

void ekf_setX(ekf_t * ekf, int i, double value)
{
    ekf->X[i] = value;
}

static void ekf_pre_update(
        ekf_t * ekf, 
        void (*f)(double *, double *, double *), 
        void (*g)(double *, double *, double *))
{
    // 1, 2
    zeros(ekf->fy, N, N);
    f(ekf->X, ekf->Xp, ekf->fy);

    // 3
    zeros(ekf->H, M, N);
    g(ekf->Xp, ekf->gXp, ekf->H);     
}

void ekf_update(
        ekf_t * ekf, 
        double * Z, 
        void (*f)(double *, double *, double *), 
        void (*g)(double *, double *, double *))
{        
    // 1,2,3
    ekf_pre_update(ekf, f, g);

    // 4,5,6,7
    ekf_post_update(ekf, Z);
}

void ekf_post_update(ekf_t * ekf, double * Z)
{    
    // 4
    mulmat(ekf->fy, ekf->P, ekf->tmp_n_n, N, N, N);
    transpose(ekf->fy, ekf->fyt, N, N);
    mulmat(ekf->tmp_n_n, ekf->fyt, ekf->Pp, N, N, N);
    add(ekf->Pp, ekf->Q, N, N);

    // 5
    transpose(ekf->H, ekf->Ht, M, N);
    mulmat(ekf->Pp, ekf->Ht, ekf->tmp_n_m, N, N, M);
    mulmat(ekf->H, ekf->Pp, ekf->tmp_m_n, M, N, N);
    mulmat(ekf->tmp_m_n, ekf->Ht, ekf->tmp2_m_m, M, N, M);
    add(ekf->tmp2_m_m, ekf->R, M, M);
    invert(ekf->tmp2_m_m, ekf->tmp_m_m, ekf->tmp_m, M);
    mulmat(ekf->tmp_n_m, ekf->tmp_m_m, ekf->G, N, M, M);

    // 6
    sub(ekf->tmp_m, ekf->gXp, Z, M);
    mulvec(ekf->G, ekf->tmp_m, ekf->X, N, M);

    // 7
    mulmat(ekf->G, ekf->H, ekf->tmp_n_n, N, M, N);
    negate(ekf->tmp_n_n, N, N);
    mat_addeye(ekf->tmp_n_n, N);
    mulmat(ekf->tmp_n_n, ekf->Pp, ekf->P, N, N, N);

    mat_dump(ekf->P, N, N, "%+10.4f"); exit(0);
}
