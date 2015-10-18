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
 * d X    |X=Xp               Jacobian of the process model
 *
 *
 * 3.
 * d g    |
 * H  = -----------|                 : Linearize observation equation, H is
 * d X    |X=Xp               the Jacobian of the measurement model
 *
 *
 * 4. fy_P = fy * Pi * fy' + Q         : Covariance of Xp
 *
 * 5. K = fy_P * H' * inv(H * fy_P * H' + R): Kalman Gain
 *
 * 6. Xo = Xp + K * (Z - g(Xp))      : Output state
 *
 * 7. Po = [I - K * H] * fy_P          : Covariance of Xo
 */


class TinyEKF {
    
public:
    
    TinyEKF(int n, int m)
    {
        this->n = n;
        this->m = m;
        
        this->P = newmat(n, n);
        this->Q = newmat(n, n);
        this->R = newmat(m, m);
        
        this->H   = newmat(m, n);
        this->fy  = newmat(n, n);

        this->X   = new double [n];
        this->Xp  = new double [n];
        this->gXp = new double[m];
        
        this->Ht   = newmat(n, m);
        this->Pp_Ht = newmat(n, m);
        this->fy_P  = newmat(n, n);
        this->fyt  = newmat(n, n);
        this->Pp   = newmat(n, n);
        this->H_Pp = newmat(m, n);
        this->H_Pp_Ht = newmat(m, m);
    }
    
    ~TinyEKF()
    {
        deletemat(this->P, this->n);
        deletemat(this->Q, this->n);
        deletemat(this->R, this->m);
        
        deletemat(this->H, this->m);
        deletemat(this->fy, n);
        
        delete this->X;
        delete this->Xp;
        delete this->gXp;
        
        deletemat(this->fy_P, this->n);
        deletemat(this->fyt, this->n);        
        deletemat(this->Pp,  this->n);
        deletemat(this->Ht, this->n);
        deletemat(this->Pp_Ht, this->n);        
        deletemat(this->H_Pp,  this->m);
        deletemat(this->H_Pp_Ht,  this->m);
    }
    
    void update(double * X)
    {        
        // 1, 2
        this->f(this->X, this->Xp, this->fy);           
        
        // 3
        this->g(this->Xp, this->gXp, this->H);     
        
        // 4
        matmul(this->fy, this->P, this->fy_P, this->n, this->n);
        transpose(this->fy, this->fyt, this->n, this->n);
        matmul(this->fy_P, this->fyt, this->Pp, this->n, this->n);
        add(this->Pp, this->Q, this->n, this->n);
        
        // 5
        transpose(this->H, this->Ht, this->m, this->n);
        matmul(this->Pp, this->Ht, this->Pp_Ht, this->n, this->m);
        matmul(this->H, this->Pp, this->H_Pp, this->m, this->n);
        matmul(this->H_Pp, this->Ht, this->H_Pp_Ht, this->m, this->n);
        dump(this->H_Pp_Ht, this->m, this->m);
        exit(0);            
    }
    
    static double ** newmat(int m, int n)
    {
        double ** a = new double * [m];
        
        for (int i=0; i<m; ++i)
            a[i] = new double [n];
        
        return a;
    }
    
    static void eye(double ** a, int n, double s)
    {
        zeros(a, n, n);
        
        for (int k=0; k<n; ++k)
            a[k][k] = s;
    }

    static void dump(double ** a, int m, int n)
    {
        for (int i=0; i<m; ++i) {
            dump(a[i], n);
        }
    }
    
protected:
    
    virtual void f(double * X, double * Xp, double ** fy) = 0;
    
    virtual void g(double * Xp, double * gXp, double ** H) = 0;
    
    static void dump(double * x, int n)
    {
        for (int j=0; j<n; ++j)
            printf("%10.6f ", x[j]);
        printf("\n");
    }
    
    
    static double deletemat(double ** a, int m)
    {
        for (int i=0; i<m; ++i)
            delete a[i];
    }
    
    static void zeros(double * a, int n)
    {
        bzero(a, n*sizeof(double));
    }
    
    static void zeros(double ** a, int m, int n)
    {
        for (int i=0; i<m; ++i)
            zeros(a[i], n);
    }
    
    static void copy(double * dst, double * src, int n)
    {
        memcpy(dst, src, n*sizeof(double));
    }
    
    static double copy(double ** dst, double ** src, int m, int n)
    {
        for (int i=0; i<m; ++i)
            copy(dst[i], src[i], n);
    }
    
    // XXX see what can be made private
    
    int n;          // state values
    int m;          // measurement values
    
    double ** P;    // covariance of prediction
    double ** Q;    // covariance of process noise
    double ** R;    // covariance of measurement noise
    
    double *  X;    // state
    double *  Xp;   // output of state-transition function
    double ** fy;   // Jacobean of process model
    double ** H;    // Jacobean of measurement model
    double *  gXp;
    
    // temporary storage
    double ** Ht;
    double ** Pp_Ht;
    double ** fy_P;
    double ** fyt;
    double ** Pp;
    double ** H_Pp;
    double ** H_Pp_Ht;
    
private:
    
    static void matmul(double ** a, double **b, double **c, int rows, int cols)
    {
        for (int i=0; i<rows; ++i)
            for (int j=0; j<cols; ++j)
                c[i][j] = dotprod(a, b, i, j, cols);
    }
    
    static double dotprod(double ** a, double ** b, int row, int col, int rowscols) 
    {
        double d = 0;
        
        for (int l=0; l<rowscols; ++l)
            d += a[row][l] * b[l][col];
                
        return d;
    }
    
    static void transpose(double ** a, double ** at, int rows, int cols)
    {
        for (int i=0; i<rows; ++i)
            for (int j=0; j<cols; ++j) 
                at[j][i] = a[i][j];
    }
    
    // A = A + B
    static void add(double ** a, double ** b, int rows, int cols) 
    {        
        for (int i=0; i<rows; ++i)
            for (int j=0; j<cols; ++j)
                a[i][j] += b[i][j];
    }    
};
