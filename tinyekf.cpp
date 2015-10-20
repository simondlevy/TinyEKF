#include<stdlib.h>

#include "tinyekf.hpp"
#include "linalg.hpp"

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

} contents_t; 

TinyEKF::TinyEKF(int n, int m)
{
    contents_t * contents = new contents_t;

    this->contents = contents;

    contents->X = newvec(n);

    contents->P = newmat(n, n);
    contents->Q = newmat(n, n);
    contents->R = newmat(m, m);
    contents->G = newmat(n, m);

    contents->H   = newmat(m, n);
    contents->fy  = newmat(n, n);

    contents->Xp  = newvec(n);
    contents->gXp = newvec(m);

    contents->Pp   = newmat(n, n);

    contents->Ht   = newmat(n, m);
    contents->fyt  = newmat(n, n);

    contents->eye_n_n = newmat(n, n);
    eye(contents->eye_n_n, 1);

    contents->tmp_m    = newvec(m);
    contents->tmp_n    = newvec(n);
    contents->tmp_m_m  = newmat(m, m);
    contents->tmp2_m_m  = newmat(m, m);
    contents->tmp_m_n  = newmat(m, n);
    contents->tmp_n_m  = newmat(n, m);
    contents->tmp2_n_m = newmat(n, m);
    contents->tmp_n_n  = newmat(n, n);
}

TinyEKF::~TinyEKF()
{
    contents_t * contents = (contents_t *)this->contents;

    deletemat(contents->P);
    deletemat(contents->Q);
    deletemat(contents->R);
    deletemat(contents->G);
    deletemat(contents->H);

    deletemat(contents->fy);
    deletemat(contents->fyt);        

    deletevec(contents->X);
    deletevec(contents->Xp);
    deletevec(contents->gXp);
    deletemat(contents->Pp);
    deletemat(contents->Ht);

    deletemat(contents->eye_n_n);

    deletemat(contents->tmp_n_n);
    deletemat(contents->tmp_m_n);

    deletemat(contents->tmp_n_m);        
    deletemat(contents->tmp2_n_m);
    deletemat(contents->tmp_m_m);
    deletevec(contents->tmp_n);
}

void TinyEKF::setP(int i, int j, double value)
{
    ((contents_t *)this->contents)->P->data[i][j] = value;
}

void TinyEKF::setQ(int i, int j, double value)
{
    ((contents_t *)this->contents)->Q->data[i][j] = value;
}

void TinyEKF::setR(int i, int j, double value)
{
    ((contents_t *)this->contents)->R->data[i][j] = value;
}

void TinyEKF::setX(int i, double value)
{
    ((contents_t *)this->contents)->X->data[i] = value;
}

void TinyEKF::update(double * Z)
{        
    contents_t * contents = (contents_t *)this->contents;

    // 1, 2
    zeros(contents->fy);
    this->f(contents->X->data, contents->Xp->data, contents->fy->data);

    // 3
    zeros(contents->H);
    this->g(contents->Xp->data, contents->gXp->data, contents->H->data);     

    // 4
    mul(contents->fy, contents->P, contents->tmp_n_n);
    transpose(contents->fy, contents->fyt);
    mul(contents->tmp_n_n, contents->fyt, contents->Pp);
    add(contents->Pp, contents->Q);

    // 5
    transpose(contents->H, contents->Ht);
    mul(contents->Pp, contents->Ht, contents->tmp_n_m);
    mul(contents->H, contents->Pp, contents->tmp_m_n);
    mul(contents->tmp_m_n, contents->Ht, contents->tmp2_m_m);
    add(contents->tmp2_m_m, contents->R);
    invert(contents->tmp2_m_m, contents->tmp_m_m);
    mul(contents->tmp_n_m, contents->tmp_m_m, contents->G);

    // 6
    contents->tmp_m->data = Z;
    sub(contents->tmp_m, contents->gXp);
    mul(contents->G, contents->tmp_m, contents->X);

    // 7
    mul(contents->G, contents->H, contents->tmp_n_n);
    negate(contents->tmp_n_n);
    add(contents->tmp_n_n, contents->eye_n_n);
    mul(contents->tmp_n_n, contents->Pp, contents->P);

    dump(contents->P, "%+10.4f"); exit(0);
}
