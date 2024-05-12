class TinyEKF {

    private:

        float x[EKF_N];    // state vector 

        float P[EKF_N][EKF_N];  // prediction error covariance 
        float Q[EKF_N][EKF_N];  // process noise covariance 
        float R[EKF_M][EKF_M];  // measurement error covariance 

        float G[EKF_N][EKF_M];  // Kalman gain; a.k.a. K 

        float F[EKF_N][EKF_N];  // Jacobian of process model 
        float H[EKF_M][EKF_N];  // Jacobian of measurement model 

        float fx[EKF_N];   // output of user defined f() state-transition function 
        float hx[EKF_M];   // output of user defined h() measurement function 

        // Cholesky-decomposition matrix-inversion code, adapated from
        // http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt

        /*
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
                            return 1; // error 
                        }
                        p[i] = sqrt(sum);
                    }
                    else {
                        a[j*n+i] = sum / p[i];
                    }
                }
            }

            return 0; // success 
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

            return 0; // success 
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

            return 0; // success 
        }
        */

        // Matrix * Matrix

        static void multiply(
                const float a[EKF_N][EKF_N], 
                const float b[EKF_N][EKF_M], 
                float c[EKF_N][EKF_M]) 
        {
            for (uint8_t i=0; i<EKF_N; i++) {

                for (uint8_t j=0; j<EKF_M; j++) {

                    c[i][j] = 0;

                    for (uint8_t k=0; k<EKF_N; k++) {

                        c[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        }

         static void multiply(
                const float a[EKF_N][EKF_N], 
                const float b[EKF_N][EKF_N], 
                float c[EKF_N][EKF_N]) 
        {
            for (uint8_t i=0; i<EKF_N; i++) {

                for (uint8_t j=0; j<EKF_N; j++) {

                    c[i][j] = 0;

                    for (uint8_t k=0; k<EKF_N; k++) {

                        c[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        }

        static void transpose(
                const float a[EKF_N][EKF_N], float at[EKF_N][EKF_N]) 
        {
            for (uint8_t i=0; i<EKF_N; ++i) {
                for (uint8_t j=0; j<EKF_N; ++j) {
                    auto tmp = a[i][j];
                    at[j][i] = a[i][j];
                }
            }
        }

        static void transpose(
                const float a[EKF_M][EKF_N], float at[EKF_N][EKF_M]) 
        {
            for (uint8_t i=0; i<EKF_M; ++i) {
                for (uint8_t j=0; j<EKF_N; ++j) {
                    at[j][i] = a[i][j];
                }
            }
        }

        // A <- A + B 
        static void accum(float a[EKF_N][EKF_N], float  b[EKF_N][EKF_N])
        {        
            for (uint8_t i=0; i<EKF_N; ++i) {
                for (uint8_t j=0; j<EKF_N; ++j) {
                    a[i][j] += b[i][j];
                }
            }
        }
        static void accum(float a[EKF_M][EKF_M], float  b[EKF_N][EKF_M])
        {        
            for (uint8_t i=0; i<EKF_M; ++i) {
                for (uint8_t j=0; j<EKF_M; ++j) {
                    a[i][j] += b[i][j];
                }
            }
        }

    public:

        void step(float v[EKF_N], float z[EKF_M])
        {        
            // Predict -------------------------------------------------------

            // P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}
            float FP[EKF_N][EKF_N] = {};
            multiply(F, P, FP);
            float Ft[EKF_N][EKF_N] = {}; 
            transpose(F, Ft);
            float Pp[EKF_N][EKF_N] = {};
            multiply(FP, Ft, Pp);
            accum(Pp, Q);

            // Update --------------------------------------------------------

            //float tmp0[EKF_N][EKF_N];
            //float tmp1[EKF_N][EKF_M];
            //float tmp2[EKF_M][EKF_N];
            //float tmp3[EKF_M][EKF_M];
            //float tmp4[EKF_M][EKF_M];
            //float tmp5[EKF_M]; 

            // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
            float Ht[EKF_N][EKF_M] = {};
            transpose(H, Ht);
            float PpHt[EKF_N][EKF_M] = {};
            multiply(Pp, Ht, PpHt);
            float HPp[EKF_M][EKF_N] = {};
            multiply(H, Pp, HPp);
            float HPpHt[EKF_M][EKF_M];
            multiply(HPp, Ht, HPpHt);
            accum(HPpHt, R);
            /*
            cholsl(HPpHt, tmp4, tmp5, m);
            mulmat(PpHt, tmp4, G, n, m, m);

            // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
            sub(z, hx, tmp5, m);
            mulvec(G, tmp5, tmp2, n, m);
            add(fx, tmp2, x, n);

            // P_k = (I - G_k H_k) P_k
            mulmat(G, H, tmp0, n, m, n);
            negate(tmp0, n, n);
            mat_addeye(tmp0, n);
            mulmat(tmp0, Pp, P, n, n, n);
            */
        }

        /**
         * Returns the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @return state value at index
         */
        float getX(int i) 
        { 
            return this->x[i]; 
        }

        /**
         * Sets the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @param value value to set
         */
        void setX(int i, float value) 
        { 
            this->x[i] = value; 
        }

};
