class TinyEKF {

    public:

        TinyEKF(int m, int n);

        void update(double * Q, double * R, double * Z, double * X, double * P);

    protected:

        virtual void f(float * x, float * fx, float * dfx) = 0;

        virtual void g(float * x, float * gx, float * dgx) = 0;
};
