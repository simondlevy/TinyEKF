typedef double number_t;

class TinyEKF {

    private:

        int n;          // number of states
        int m;          // number of measurements

        number_t * P;     // prediction error covariance
        number_t * Q;     // process noise covariance
        number_t * R;     // measurement error covariance

        number_t * G;     // Kalman gain; a.k.a. K
        number_t * F;     // Jacobian of process model
        number_t * H;     // Jacobian of measurement model

        number_t * Ht;    // transpose of measurement Jacobian
        number_t * Ft;    // transpose of process Jacobian
        number_t * Pp;    // P, post-prediction, pre-update

        number_t * fx;    // output of user defined f() state-transition function
        number_t * hx;    // output of user defined h() measurement function

        // temporary storage
        number_t * tmp1;
        number_t * tmp2;
        number_t * tmp3;
        number_t * tmp4;
        number_t * tmp5;

    protected:

        number_t * x;     // state

        TinyEKF(int n, int m);

        ~TinyEKF();

        virtual void model(number_t * fx, number_t * F, number_t * hx, number_t * H) = 0;

        void set(number_t * A, int i, int j, number_t value);

        void setP(int i, int j, number_t value);

        void setQ(int i, int j, number_t value);

        void setR(int i, int j, number_t value);

    public:

        number_t getX(int i);

        void step(number_t * z);
};
