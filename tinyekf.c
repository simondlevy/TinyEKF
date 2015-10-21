#include<stdlib.h>

#include "linalg.h"
#include "tinyekf.h"

void ekf_init(ekf_t * ekf, int n, int m)
{
    vec_init(&ekf->X, n);

    ekf->P = newmat(n, n);
    ekf->Q = newmat(n, n);
    ekf->R = newmat(m, m);
    ekf->G = newmat(n, m);

    ekf->H   = newmat(m, n);
    ekf->fy  = newmat(n, n);

    ekf->Xp  = newvec(n);
    ekf->gXp = newvec(m);

    ekf->Pp   = newmat(n, n);

    ekf->Ht   = newmat(n, m);
    ekf->fyt  = newmat(n, n);

    ekf->eye_n_n = newmat(n, n);
    eye(ekf->eye_n_n, 1);

    ekf->tmp_m    = newvec(m);
    ekf->tmp_n    = newvec(n);
    ekf->tmp_m_m  = newmat(m, m);
    ekf->tmp2_m_m  = newmat(m, m);
    ekf->tmp_m_n  = newmat(m, n);
    ekf->tmp_n_m  = newmat(n, m);
    ekf->tmp2_n_m = newmat(n, m);
    ekf->tmp_n_n  = newmat(n, n);
}

void ekf_delete(ekf_t * ekf)
{
    deletemat(ekf->P);
    deletemat(ekf->Q);
    deletemat(ekf->R);
    deletemat(ekf->G);
    deletemat(ekf->H);

    deletemat(ekf->fy);
    deletemat(ekf->fyt);        

    vec_free(ekf->X);

    deletevec(ekf->Xp);
    deletevec(ekf->gXp);
    deletemat(ekf->Pp);
    deletemat(ekf->Ht);

    deletemat(ekf->eye_n_n);

    deletemat(ekf->tmp_n_n);
    deletemat(ekf->tmp_m_n);

    deletemat(ekf->tmp_n_m);        
    deletemat(ekf->tmp2_n_m);
    deletemat(ekf->tmp_m_m);
    deletevec(ekf->tmp_n);
}

void ekf_setP(ekf_t * ekf, int i, int j, double value)
{
    ekf->P->data[i][j] = value;
}

void ekf_setQ(ekf_t * ekf, int i, int j, double value)
{
    ekf->Q->data[i][j] = value;
}

void ekf_setR(ekf_t * ekf, int i, int j, double value)
{
    ekf->R->data[i][j] = value;
}

void ekf_setX(ekf_t * ekf, int i, double value)
{
    ekf->X.data[i] = value;
}

static void ekf_pre_update(
        ekf_t * ekf, 
        void (*f)(double *, double *, double **), 
        void (*g)(double *, double *, double **))
{
    // 1, 2
    zeros(ekf->fy);
    f(ekf->X.data, ekf->Xp->data, ekf->fy->data);

    // 3
    zeros(ekf->H);
    g(ekf->Xp->data, ekf->gXp->data, ekf->H->data);     
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
    ekf->tmp_m->data = Z;
    sub(ekf->tmp_m, ekf->gXp);
    mulvec(ekf->G, ekf->tmp_m, &ekf->X);

    // 7
    mulmat(ekf->G, ekf->H, ekf->tmp_n_n);
    negate(ekf->tmp_n_n);
    add(ekf->tmp_n_n, ekf->eye_n_n);
    mulmat(ekf->tmp_n_n, ekf->Pp, ekf->P);
}
