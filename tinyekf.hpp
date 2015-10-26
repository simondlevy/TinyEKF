#include "tinyekf.h"

#include <stdlib.h>
#include <strings.h>

class TinyEKF {

    private:

        ekf_t ekf;

    protected:

        virtual void f(double X[N], double Xp[N], double F[N][N]) = 0;

        virtual void g(double Xp[N], double gXp[N], double H[M][N]) = 0;    

    public:

        void setP(int i, int j, double value)
        {
            this->ekf.P[i][j] = value;
        }

        void setQ(int i, int j, double value)
        {
            this->ekf.Q[i][j] = value;
        }

        void setR(int i, int j, double value)
        {
            this->ekf.R[i][j] = value;
        }

        void setX(int i, double value)
        {
            this->ekf.X[i] = value;
        }

        void update(double * Z)
        {        
            // 1, 2
            this->f(this->ekf.X, this->ekf.Xp, this->ekf.F); 

            // 3
            this->g(this->ekf.Xp, this->ekf.gXp, this->ekf.H);     
 
            // 4,5,6,7
            ekf_post_update(&this->ekf, Z);
        }
};
