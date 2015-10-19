/*
 * TinyEKF - Extended Kalman Filter in a small amount of code
 *
 * Based on Extended_KF.m by Chong You
 * https://sites.google.com/site/chongyou1987/
 *
 * Syntax:
 * [Xo,Po] = Extended_KF(f,g,Q,R,Z,Xi,Pi)
 *
 * State Equation:
 * X(n+1) = f(X(n)) + w(n)
 * where the state X has the dimension N-by-1
 * Observation Equation:
 * Z(n) = g(X(n)) + v(n)
 * where the observation y has the dimension M-by-1
 * w ~ N(0,Q) is gaussian noise with covariance Q
 * v ~ N(0,R) is gaussian noise with covariance R
 * Input:
 * f: function for state transition, it takes a state variable Xn and
 * returns 1) f(Xn) and 2) Jacobian of f at Xn. As a fake example:
 * function [Val, Jacob] = f(X)
 * Val = sin(X) + 3;
 * Jacob = cos(X);
 * end
 * g: function for measurement, it takes the state variable Xn and
 * returns 1) g(Xn) and 2) Jacobian of g at Xn.
 * Q: process noise covariance matrix, N-by-N
 * R: measurement noise covariance matrix, M-by-M
 *
 * Z: current measurement, M-by-1
 *
 * Xi: "a priori" state estimate, N-by-1
 * Pi: "a priori" estimated state covariance, N-by-N
 * Output:
 * Xo: "a posteriori" state estimate, N-by-1
 * Po: "a posteriori" estimated state covariance, N-by-N
 *
 * Algorithm for Extended Kalman Filter:
 * Linearize input functions f and g to get fy(state transition matrix)
 * and H(observation matrix) for an ordinary Kalman Filter:
 * State Equation:
 * X(n+1) = fy * X(n) + w(n)
 * Observation Equation:
 * Z(n) = H * X(n) + v(n)
 *
 * 1. Xp = f(Xi)                     : One step projection, also provides
 * linearization point
 *
 * 2.
 * d f    |
 * fy = -----------|                 : Linearize state equation, fy is the
 * d X    |X=Xp                        Jacobian of the process model
 *
 *
 * 3.
 * d g    |
 * H  = -----------|                 : Linearize observation equation, H is
 * d X    |X=Xp                        the Jacobian of the measurement model
 *
 *
 * 4. Pp = fy * Pi * fy' + Q         : Covariance of Xp
 *
 * 5. K = Pp * H' * inv(H * Pp * H' + R): Kalman Gain
 *
 * 6. Xo = Xp + K * (Z - g(Xp))      : Output state
 *
 * 7. Po = [I - K * H] * Pp          : Covariance of Xo
 */

#include "linalg.hpp"

class TinyEKF {

protected:

    double *  X;    // state
    double ** P;    // covariance of prediction
    double ** Q;    // covariance of process noise
    double ** R;    // covariance of measurement noise
    
    virtual void f(double * Xp, double ** fy) = 0;
    
    virtual void g(double * Xp, double * gXp, double ** H) = 0;    

    TinyEKF(int n, int m)
    {
        this->n = n;
        this->m = m;
        
        this->P = newmat(n, n);
        this->Q = newmat(n, n);
        this->R = newmat(m, m);
        this->G = newmat(n, m);
        
        this->H   = newmat(m, n);
        this->fy  = newmat(n, n);

        this->X   = new double [n];
        this->Xp  = new double [n];
        this->gXp = new double[m];

        this->Pp   = newmat(n, n);
        
        this->Ht   = newmat(n, m);
        this->fyt  = newmat(n, n);

        this->tmp_n    = new double [n];
        this->tmp_m_m  = newmat(m, m);
        this->tmp_m_n  = newmat(m, n);
        this->tmp_n_m  = newmat(n, m);
        this->tmp2_n_m = newmat(n, m);
        this->tmp_n_n  = newmat(n, n);
    }
    
   ~TinyEKF()
    {
        deletemat(this->P, this->n);
        deletemat(this->Q, this->n);
        deletemat(this->R, this->m);
        deletemat(this->G, this->n);
        
        deletemat(this->H, this->m);
        deletemat(this->fy, n);
        
        delete this->X;
        delete this->Xp;
        delete this->gXp;
        
        deletemat(this->tmp_n_n, this->n);
        deletemat(this->fyt, this->n);        
        deletemat(this->Pp,  this->n);
        deletemat(this->Ht, this->n);
        deletemat(this->tmp_m_n,  this->m);

        deletemat(this->tmp_n_m, this->n);        
        deletemat(this->tmp2_n_m,  this->m);
        deletemat(this->tmp_m_m,  this->m);
        delete this->tmp_n;
    }
  
private:
    
    int n;          // state values
    int m;          // measurement values
    
    double ** G;    // Kalman gain; a.k.a. K
    
    double *  Xp;   // output of state-transition function
    double ** fy;   // Jacobean of process model
    double ** H;    // Jacobean of measurement model
    double *  gXp;
    
    // temporary storage
    double ** Ht;
    double ** tmp_n_m;
    double ** tmp_n_n;
    double ** fyt;
    double ** Pp;
    double ** tmp_m_n;

    double  * tmp_n;
    double ** tmp2_n_m;
    double ** tmp_m_m;

public:
   
    void update(double * Z)
    {        
        // 1, 2
        this->f(this->Xp, this->fy);           
        
        // 3
        this->g(this->Xp, this->gXp, this->H);     
        
        // 4
        matmul(this->fy, this->P, this->tmp_n_n, this->n, this->n, this->n);
        transpose(this->fy, this->fyt, this->n, this->n);
        matmul(this->tmp_n_n, this->fyt, this->Pp, this->n, this->n, this->n);
        add(this->Pp, this->Q, this->n, this->n);

        // 5
        transpose(this->H, this->Ht, this->m, this->n);
        matmul(this->Pp, this->Ht, this->tmp_n_m, this->n, this->m, this->n);
        matmul(this->H, this->Pp, this->tmp_m_n, this->m, this->n, this->n);
        matmul(this->tmp_m_n, this->Ht, this->tmp2_n_m, this->m, this->m, this->n);
        add(this->tmp2_n_m, this->R, this->m, this->m);
        invert(this->tmp2_n_m, this->tmp_m_m, this->tmp_n, this->m);
        matmul(this->tmp_n_m, this->tmp_m_m, this->G, this->n, this->m, this->m);

        dump(this->G, this->n, this->m); exit(0);

        // 6
        //Xo = Xp + K * (Z - gXp);
    }
};
