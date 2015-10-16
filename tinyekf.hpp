class TinyEKF {

    public:

        TinyEKF(int m, int n);

        void update(double * Q, double * R, double * Z, double * X, double * P);
};
