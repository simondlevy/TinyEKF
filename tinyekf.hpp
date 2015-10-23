#include "tinyekf.h"

#include <stdlib.h>

class TinyEKF {

    private:

        ekf_t ekf;

        double fy[N][N];

    protected:

        virtual void f(double X[N], double Xp[N], double fy[N][N]) = 0;

        virtual void g(double Xp[N], double gXp[N], double H[M*N]) = 0;    

    public:

        void setP(int i, int j, double value)
        {
            ekf_setP(&this->ekf, i, j, value);
        }

        void setQ(int i, int j, double value)
        {
            ekf_setQ(&this->ekf, i, j, value);
        }

        void setR(int i, int j, double value)
        {
            ekf_setR(&this->ekf, i, j, value);
        }

        void setX(int i, double value)
        {
            ekf_setX(&this->ekf, i, value);
        }

        void update(double * Z)
        {        
            ekf_t ekf = this->ekf;

            // 1, 2
            this->f(ekf.X, ekf.Xp, this->fy); 
            memcpy(ekf.fy, this->fy, N*N*sizeof(double));

            // 3
            this->g(ekf.Xp, ekf.gXp, ekf.H);     

            // 4,5,6,7
            ekf_post_update(&ekf, Z);
        }
};
