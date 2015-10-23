#include "tinyekf.h"

class TinyEKF {

    private:

        ekf_t ekf;

    protected:

        virtual void f(double * X, double * Xp, double * fy) = 0;

        virtual void g(double * Xp, double * gXp, double * H) = 0;    

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
            this->f(ekf.X, ekf.Xp, ekf.fy);

            // 3
            this->g(ekf.Xp, ekf.gXp, ekf.H);     

            // 4,5,6,7
            ekf_post_update(&ekf, Z);
        }
};
