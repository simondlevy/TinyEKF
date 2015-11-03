#include <stdlib.h>
#include <strings.h>

#define N 8
#define M 4

typedef double number_t;

class TinyEKF {

    private:

        int n;
        int m;

        number_t x[N];      // state

        number_t P[N][N];   // prediction error covariance
        number_t Q[N][N];   // process noise covariance
        number_t R[M][M];   // measurement error covariance

        number_t G[N][M];   // Kalman gain; a.k.a. K
        number_t F[N][N];   // Jacobian of process model
        number_t H[M][N];   // Jacobian of measurement model

        number_t Ht[N][M];  // transpose of measurement Jacobian
        number_t Ft[N][N];  // transpose of process Jacobian
        number_t Pp[N][N];  // P, post-prediction, pre-update

        number_t fx[N];     // f(x)
        number_t hx[N];     // h(x)

        // temporary storage
        number_t tmp1[N*N];
        number_t tmp2[M*N];
        number_t tmp3[M*M];
        number_t tmp4[M*M];
        number_t tmp5[M];

    protected:

        TinyEKF(int n, int m);

        ~TinyEKF();

        virtual void f(double x[N], double fx[N], double F[N][N]) = 0;

        virtual void h(double fx[N], double hx[N], double H[M][N]) = 0;    

    public:

        void setP(int i, int j, double value);

        void setQ(int i, int j, double value);

        void setR(int i, int j, double value);

        void setX(int i, double value);

        double getX(int i);

        void step(double * Z);
};
