#include "tinyekf.h"

#include <stdlib.h>
#include <strings.h>

class TinyEKF {

    private:

        ekf_t ekf;

    protected:

        virtual void f(double x[N], double F[N][N]) = 0;

        virtual void h(double x[N], double hx[N], double H[M][N]) = 0;    

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
            this->ekf.x[i] = value;
        }

        void step(double * Z)
        {        
            // Model
            this->f(this->ekf.x, this->ekf.F); 
            this->h(this->ekf.x, this->ekf.hx, this->ekf.H);     
 
            // Predict, update
            ekf_predict_and_update(&this->ekf, Z);
        }
};
