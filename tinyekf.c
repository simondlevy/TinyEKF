#include<stdlib.h>

#include "tinyekf.h"

void ekf_init(ekf_t * ekf, int n, int m)
{
    mat_init(&ekf->P,n, n);
    mat_init(&ekf->Q,n, n);
    mat_init(&ekf->R,m, m);
    mat_init(&ekf->G,n, m);

    mat_init(&ekf->H  ,m, n);
    mat_init(&ekf->fy ,n, n);

    mat_init(&ekf->Pp  ,n, n);

    mat_init(&ekf->Ht  ,n, m);
    mat_init(&ekf->fyt ,n, n);

    mat_init(&ekf->eye_n_n,n, n);
    eye(ekf->eye_n_n, N, 1);

    mat_init(&ekf->tmp_m_m, m, m);
    mat_init(&ekf->tmp2_m_m, m, m);
    mat_init(&ekf->tmp_m_n, m, n);
    mat_init(&ekf->tmp_n_m, n, m);
    mat_init(&ekf->tmp2_n_m, n, m);
    mat_init(&ekf->tmp_n_n, n, n);
}

void ekf_free(ekf_t ekf)
{
    mat_free(ekf.P);
    mat_free(ekf.Q);
    mat_free(ekf.R);
    mat_free(ekf.G);
    mat_free(ekf.H);

    mat_free(ekf.fy);
    mat_free(ekf.fyt);        

    mat_free(ekf.Pp);
    mat_free(ekf.Ht);

    mat_free(ekf.eye_n_n);

    mat_free(ekf.tmp_n_n);
    mat_free(ekf.tmp_m_n);

    mat_free(ekf.tmp_n_m);        
    mat_free(ekf.tmp2_n_m);
    mat_free(ekf.tmp_m_m);
}

void ekf_setP(ekf_t * ekf, int i, int j, double value)
{
    mat_set(ekf->P, i, j, value);
}

void ekf_setQ(ekf_t * ekf, int i, int j, double value)
{
    mat_set(ekf->Q, i, j, value);
}

void ekf_setR(ekf_t * ekf, int i, int j, double value)
{
    mat_set(ekf->R, i, j, value);
}

void ekf_setX(ekf_t * ekf, int i, double value)
{
    ekf->X[i] = value;
}

static void ekf_pre_update(
        ekf_t * ekf, 
        void (*f)(double *, double *, double **), 
        void (*g)(double *, double *, double **))
{
    // 1, 2
    zeros(ekf->fy, N, N);
    f(ekf->X, ekf->Xp, ekf->fy.data);

    // 3
    zeros(ekf->H, M, N);
    g(ekf->Xp, ekf->gXp, ekf->H.data);     
}

void ekf_update(
        ekf_t * ekf, 
        double * Z, 
        void (*f)(double *, double *, double **), 
        void (*g)(double *, double *, double **))
{        
    // 1,2,3
    ekf_pre_update(ekf, f, g);

    // 4,5,6,7
    ekf_post_update(ekf, Z);
}

void ekf_post_update(ekf_t * ekf, double * Z)
{    
    // 4
    mulmat(ekf->fy, ekf->P, ekf->tmp_n_n);
    transpose(ekf->fy, ekf->fyt);
    mulmat(ekf->tmp_n_n, ekf->fyt, ekf->Pp);
    add(ekf->Pp, ekf->Q);

    // 5
    transpose(ekf->H, ekf->Ht);
    mulmat(ekf->Pp, ekf->Ht, ekf->tmp_n_m);
    mulmat(ekf->H, ekf->Pp, ekf->tmp_m_n);
    mulmat(ekf->tmp_m_n, ekf->Ht, ekf->tmp2_m_m);
    add(ekf->tmp2_m_m, ekf->R);
    invert(ekf->tmp2_m_m, ekf->tmp_m_m);
    mulmat(ekf->tmp_n_m, ekf->tmp_m_m, ekf->G);

    // 6
    sub(ekf->tmp_m, ekf->gXp, Z, M);
    mulvec(ekf->G, ekf->tmp_m, ekf->X);

    // 7
    mulmat(ekf->G, ekf->H, ekf->tmp_n_n);
    negate(ekf->tmp_n_n);
    add(ekf->tmp_n_n, ekf->eye_n_n);
    mulmat(ekf->tmp_n_n, ekf->Pp, ekf->P);

    mat_dump(ekf->P, N, N, "%+10.4f"); exit(0);
}
