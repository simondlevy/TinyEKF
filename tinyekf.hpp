/*
TinyEKF - Extended Kalman Filter in a small amount of code

 Based on Extended_KF.m by Chong You
 https://sites.google.com/site/chongyou1987/

 Syntax:
     [Xo,Po] = Extended_KF(f,g,Q,R,Z,Xi,Pi)

 State Equation:
     X(n+1) = f(X(n)) + w(n)
     where the state X has the dimension N-by-1
 Observation Equation:
     Z(n) = g(X(n)) + v(n)
     where the observation y has the dimension M-by-1
     w ~ N(0,Q) is gaussian noise with covariance Q
     v ~ N(0,R) is gaussian noise with covariance R     
 Input:
     f: function for state transition, it takes a state variable Xn and
       returns 1) f(Xn) and 2) Jacobian of f at Xn. As a fake example:
           function [Val, Jacob] = f(X)
           Val = sin(X) + 3;
           Jacob = cos(X);
           end
     g: function for measurement, it takes the state variable Xn and
       returns 1) g(Xn) and 2) Jacobian of g at Xn.
     Q: process noise covariance matrix, N-by-N
     R: measurement noise covariance matrix, M-by-M
     
     Z: current measurement, M-by-1
     
     Xi: "a priori" state estimate, N-by-1
     Pi: "a priori" estimated state covariance, N-by-N
 Output:
     Xo: "a posteriori" state estimate, N-by-1
     Po: "a posteriori" estimated state covariance, N-by-N

 Algorithm for Extended Kalman Filter:
 Linearize input functions f and g to get fy(state transition matrix)
 and H(observation matrix) for an ordinary Kalman Filter:
 State Equation:
     X(n+1) = fy * X(n) + w(n)
 Observation Equation:
     Z(n) = H * X(n) + v(n)

 1. Xp = f(Xi)                     : One step projection, also provides 
                                     linearization point
 
 2. 
          d f    |
 fy = -----------|                 : Linearize state equation, fy is the
          d X    |X=Xp               Jacobian of the process model
       
 
 3.
          d g    |
 H  = -----------|                 : Linearize observation equation, H is
          d X    |X=Xp               the Jacobian of the measurement model
             
       
 4. Pp = fy * Pi * fy' + Q         : Covariance of Xp
 
 5. K = Pp * H' * inv(H * Pp * H' + R): Kalman Gain
 
 6. Xo = Xp + K * (Z - g(Xp))      : Output state
 
 7. Po = [I - K * H] * Pp          : Covariance of Xo
*/

static void dump(double * a, int m, int n)
{
    for (int j=0; j<m; ++j) {
        for (int k=0; k<n; ++k)
            printf("%10.6e ", a[j*n+k]);
        printf("\n");
    }
}

class TinyEKF {

    public:

        TinyEKF(int n, int m)
        {
            this->Xp  = new double [n];
            this->fy  = new double [n*n];
            this->gXp = new double[m];
            this->H   = new double[m*n];
        }
        
        ~TinyEKF()
        {
            delete this->Xp;
        }

        void update(double * Q, double * R, double * Z, double * X, double * P)
        {
            this->f(X, this->Xp);                      // 1
             
            this->df(this->Xp, this->fy);              // 2
            
            this->g(this->Xp, this->gXp, this->H);     // 3

            dump(this->gXp, 4, 1);
            exit(0);
        }

    protected:

        virtual void f(double * x, double * fx) = 0;
        
        virtual void df(double * x, double * dfx) = 0;

        virtual void g(double * xp, double * gx, double * dgx) = 0;
        
    private:
        
        double * Xp;  // output of state-transition function
        double * fy;  // Jacobean of process model
        double * H;   // Jacobean of measurement model
        double * gXp; 
};
