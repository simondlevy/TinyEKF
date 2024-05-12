/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

class TinyEKF {

    public:

        /**
         * Returns the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @return state value at index
         */
        float getX(int i) 
        { 
            return x[i]; 
        }

        /**
         * Sets the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @param value value to set
         */
        void setX(int i, float value) 
        { 
            x[i] = value; 
        }

        /**
          Performs one step of the prediction and update.
         * @param z observation vector, length <i>m</i>
         * @return true on success, false on failure caused by
         * non-positive-definite matrix.
         */
        bool step(float * z) 
        { 
            float fx[EKF_N] = {};
            float F[EKF_N][EKF_N] = {};
            float hx[EKF_M] = {};
            float H[EKF_M][EKF_N] = {};

            this->model(fx, F, hx, H); 

            return ekf_step(z, fx, F, hx, H) ? false : true;
        }

    private:

        // Cholesky-decomposition matrix-inversion code, adapated from
        // http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt

        static int choldc1(float * a, float * p, int n) {
            int i,j,k;
            float sum;

            for (i = 0; i < n; i++) {
                for (j = i; j < n; j++) {
                    sum = a[i*n+j];
                    for (k = i - 1; k >= 0; k--) {
                        sum -= a[i*n+k] * a[j*n+k];
                    }
                    if (i == j) {
                        if (sum <= 0) {
                            return 1; /* error */
                        }
                        p[i] = sqrt(sum);
                    }
                    else {
                        a[j*n+i] = sum / p[i];
                    }
                }
            }

            return 0; /* success */
        }

        static int choldcsl(float * A, float * a, float * p, int n) 
        {
            int i,j,k; float sum;
            for (i = 0; i < n; i++) 
                for (j = 0; j < n; j++) 
                    a[i*n+j] = A[i*n+j];
            if (choldc1(a, p, n)) return 1;
            for (i = 0; i < n; i++) {
                a[i*n+i] = 1 / p[i];
                for (j = i + 1; j < n; j++) {
                    sum = 0;
                    for (k = i; k < j; k++) {
                        sum -= a[j*n+k] * a[k*n+i];
                    }
                    a[j*n+i] = sum / p[j];
                }
            }

            return 0; /* success */
        }


        static int cholsl(float * A, float * a, float * p, int n) 
        {
            int i,j,k;
            if (choldcsl(A,a,p,n)) return 1;
            for (i = 0; i < n; i++) {
                for (j = i + 1; j < n; j++) {
                    a[i*n+j] = 0.0;
                }
            }
            for (i = 0; i < n; i++) {
                a[i*n+i] *= a[i*n+i];
                for (k = i + 1; k < n; k++) {
                    a[i*n+i] += a[k*n+i] * a[k*n+i];
                }
                for (j = i + 1; j < n; j++) {
                    for (k = j; k < n; k++) {
                        a[i*n+j] += a[k*n+i] * a[k*n+j];
                    }
                }
            }
            for (i = 0; i < n; i++) {
                for (j = 0; j < i; j++) {
                    a[i*n+j] = a[j*n+i];
                }
            }

            return 0; /* success */
        }

        /* C <- A * B */
        static void mulmat(
                float * a, float * b, float * c, 
                int arows, int acols, int bcols)
        {
            int i, j,l;

            for(i=0; i<arows; ++i)
                for(j=0; j<bcols; ++j) {
                    c[i*bcols+j] = 0;
                    for(l=0; l<acols; ++l)
                        c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
                }
        }

        static void mulvec(float * a, float * x, float * y, int m, int n)
        {
            int i, j;

            for(i=0; i<m; ++i) {
                y[i] = 0;
                for(j=0; j<n; ++j)
                    y[i] += x[j] * a[i*n+j];
            }
        }

        static void transpose(float * a, float * at, int m, int n)
        {
            int i,j;

            for(i=0; i<m; ++i)
                for(j=0; j<n; ++j) {
                    at[j*m+i] = a[i*n+j];
                }
        }

        /* A <- A + B */
        static void accum(float * a, float * b, int m, int n)
        {        
            int i,j;

            for(i=0; i<m; ++i)
                for(j=0; j<n; ++j)
                    a[i*n+j] += b[i*n+j];
        }

        /* C <- A + B */
        static void add(float * a, float * b, float * c, int n)
        {
            int j;

            for(j=0; j<n; ++j)
                c[j] = a[j] + b[j];
        }


        /* C <- A - B */
        static void sub(float * a, float * b, float * c, int n)
        {
            int j;

            for(j=0; j<n; ++j)
                c[j] = a[j] - b[j];
        }

        static void negate(float * a, int m, int n)
        {        
            int i, j;

            for(i=0; i<m; ++i)
                for(j=0; j<n; ++j)
                    a[i*n+j] = -a[i*n+j];
        }

        static void mat_addeye(float * a, int n)
        {
            int i;
            for (i=0; i<n; ++i)
                a[i*n+i] += 1;
        }

        int ekf_step(
                float * z, 
                float fx[EKF_N], 
                float F[EKF_N][EKF_N],
                float hx[EKF_N],
                float H[EKF_M][EKF_N])
        {        
            float _F[EKF_N*EKF_N] = {};
            memcpy(_F, F, EKF_N*EKF_N*sizeof(float));

            float _H[EKF_M*EKF_N] = {};
            memcpy(_H, H, EKF_M*EKF_N*sizeof(float));

            float tmp0[EKF_N*EKF_N] = {};
            float tmp1[EKF_N*EKF_M] = {};
            float tmp2[EKF_M*EKF_N] = {};
            float tmp3[EKF_M*EKF_M] = {};
            float tmp4[EKF_M*EKF_M] = {};
            float tmp5[EKF_M] = {}; 

            float G[EKF_N*EKF_M] = {};  
            float Ht[EKF_N*EKF_M] = {}; 
            float Ft[EKF_N*EKF_N] = {};
            float Pp[EKF_N*EKF_N]; 

            /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
            mulmat(_F, P, tmp0, EKF_N,EKF_N, EKF_N);
            transpose(_F, Ft, EKF_N, EKF_N);
            mulmat(tmp0, Ft, Pp, EKF_N,EKF_N, EKF_N);
            accum(Pp, Q, EKF_N, EKF_N);

            /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
            transpose(_H, Ht, EKF_M, EKF_N);
            mulmat(Pp, Ht, tmp1, EKF_N,EKF_N, EKF_M);
            mulmat(_H, Pp, tmp2, EKF_M,EKF_N, EKF_N);
            mulmat(tmp2, Ht, tmp3, EKF_M,EKF_N, EKF_M);
            accum(tmp3, R, EKF_M, EKF_M);
            if (cholsl(tmp3, tmp4, tmp5, EKF_M)) return 1;
            mulmat(tmp1, tmp4, G, EKF_N,EKF_M, EKF_M);

            /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
            sub(z, hx, tmp5, EKF_M);
            mulvec(G, tmp5, tmp2, EKF_N, EKF_M);
            add(fx, tmp2, x, EKF_N);

            /* P_k = (I - G_k H_k) P_k */
            mulmat(G, _H, tmp0, EKF_N,EKF_M, EKF_N);
            negate(tmp0, EKF_N, EKF_N);
            mat_addeye(tmp0, EKF_N);
            mulmat(tmp0, Pp, P, EKF_N,EKF_N, EKF_N);

            /* success */
            return 0;
        }

    protected:

        /**
         * The current state.
         */
        float x[EKF_N];

        float P[EKF_N * EKF_N];

        float Q[EKF_N * EKF_N];  

        float R[EKF_M * EKF_M];  

        /**
         * Initializes a TinyEKF object.
         */
        TinyEKF() { 
        }

        /**
         * Deallocates memory for a TinyEKF object.
         */
        ~TinyEKF() { }

        /**
         * Implement this function for your EKF model.
         * @param fx gets output of state-transition function <i>f(x<sub>0 ..
         * n-1</sub>)</i> @param F gets <i>n &times; n</i> Jacobian of
         * <i>f(x)</i>
         * @param hx gets output of observation function <i>h(x<sub>0 ..
         * n-1</sub>)</i> @param H gets <i>m &times; n</i> Jacobian of
         * <i>h(x)</i>
         */
        virtual void model(
                float fx[EKF_N], 
                float F[EKF_N][EKF_N], 
                float hx[EKF_M], 
                float H[EKF_M][EKF_N]) = 0;

        /**
         * Sets the specified value of the prediction error covariance.
         * <i>P<sub>i,j</sub> = value</i> @param i row index
         * @param j column index
         * @param value value to set
         */
        void setP(int i, int j, float value) 
        { 
            P[i * EKF_N + j] = value; 
        }

        /** Sets the specified value of the process noise covariance.
         * <i>Q<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setQ(int i, int j, float value) 
        { 
            Q[i * EKF_N + j] = value; 
        }

        /**
         * Sets the specified value of the observation noise covariance.
         * <i>R<sub>i,j</sub> = value</i> @param i row index
         * @param j column index
         * @param value value to set
         */
        void setR(int i, int j, float value) 
        { 
            R[i * EKF_M + j] = value; 
        }

};
