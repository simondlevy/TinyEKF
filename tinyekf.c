#include<stdlib.h>

#include "linalg.h"
#include "tinyekf.h"

ekf_t * ekf_init(int m, int n)
{
    ekf_t * ekf = (ekf_t *)malloc(sizeof(ekf_t));

    ekf->X = newvec(n);

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

    return ekf;
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

    deletevec(ekf->X);
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

    free(ekf);
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
    ekf->X->data[i] = value;
}

void ekf_update(
        ekf_t * ekf, 
        double * Z, 
        void (*f)(double *, double *, double **), 
        void (*g)(double *, double *, double **))
{        
    // 1, 2
    zeros(ekf->fy);
    f(ekf->X->data, ekf->Xp->data, ekf->fy->data);

    // 3
    zeros(ekf->H);
    g(ekf->Xp->data, ekf->gXp->data, ekf->H->data);     

    // 4
    mul(ekf->fy, ekf->P, ekf->tmp_n_n);
    transpose(ekf->fy, ekf->fyt);
    mul(ekf->tmp_n_n, ekf->fyt, ekf->Pp);
    add(ekf->Pp, ekf->Q);

    // 5
    transpose(ekf->H, ekf->Ht);
    mul(ekf->Pp, ekf->Ht, ekf->tmp_n_m);
    mul(ekf->H, ekf->Pp, ekf->tmp_m_n);
    mul(ekf->tmp_m_n, ekf->Ht, ekf->tmp2_m_m);
    add(ekf->tmp2_m_m, ekf->R);
    invert(ekf->tmp2_m_m, ekf->tmp_m_m);
    mul(ekf->tmp_n_m, ekf->tmp_m_m, ekf->G);

    // 6
    ekf->tmp_m->data = Z;
    sub(ekf->tmp_m, ekf->gXp);
    mul(ekf->G, ekf->tmp_m, ekf->X);

    // 7
    mul(ekf->G, ekf->H, ekf->tmp_n_n);
    negate(ekf->tmp_n_n);
    add(ekf->tmp_n_n, ekf->eye_n_n);
    mul(ekf->tmp_n_n, ekf->Pp, ekf->P);

    dumpmat(ekf->P, "%+10.4f"); exit(0);
}
