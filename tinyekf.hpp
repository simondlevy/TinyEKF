#include<stdlib.h>
#include<stdio.h>

#include "linalg.h"
#include "tinyekf.h"

class TinyEKF {

    private:

        ekf_t * ekf;

    protected:

        virtual void f(double * X, double * Xp, double ** fy) = 0;

        virtual void g(double * Xp, double * gXp, double ** H) = 0;    


    public:

        TinyEKF(int n, int m)
        {
            this->ekf = ekf_init(n, m);
        }

        ~TinyEKF()
        {
            ekf_delete((ekf_t *)this->ekf);
        }

        void setP(int i, int j, double value)
        {
            ekf_setP((ekf_t *)this->ekf, i, j, value);
        }

        void setQ(int i, int j, double value)
        {
            ekf_setQ((ekf_t *)this->ekf, i, j, value);
        }

        void setR(int i, int j, double value)
        {
            ekf_setR((ekf_t *)this->ekf, i, j, value);
        }

        void setX(int i, double value)
        {
            ekf_setX((ekf_t *)this->ekf, i, value);
        }

        void update(double * Z)
        {        
            ekf_t * ekf = (ekf_t *)this->ekf;

            // 1, 2
            zeros(ekf->fy);
            this->f(ekf->X->data, ekf->Xp->data, ekf->fy->data);

            // 3
            zeros(ekf->H);
            this->g(ekf->Xp->data, ekf->gXp->data, ekf->H->data);     

            // 4,5,6,7
            ekf_post_update(ekf, Z);

            dumpmat(ekf->P, "%+10.4f"); exit(0);
        }

};
