/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

class TinyEKF {

    public:

        void initialize(const float Pdiag[EKF_N])
        {
        }

        /**
         * Returns the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @return state value at index
         */
        float get(const int i) 
        { 
            return _x[i]; 
        }

         /**
          Performs one step of the prediction and update.
         * @param z observation vector, length <i>m</i>
         * @return true on success, false on failure caused by
         * non-positive-definite matrix.
         */
        bool step(
                const float fx[EKF_N],
                const float F[EKF_N*EKF_N],
                const float Q[EKF_N*EKF_N],
                const float hx[EKF_M],
                const float H[EKF_M*EKF_N],
                const float R[EKF_M*EKF_M],
                const float z[EKF_M]) 
        { 
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

            // P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}
            mulmat(F, _P, tmp0, EKF_N,EKF_N, EKF_N);
            transpose(F, Ft, EKF_N, EKF_N);
            mulmat(tmp0, Ft, Pp, EKF_N,EKF_N, EKF_N);
            accum(Pp, Q, EKF_N, EKF_N);

            // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
            transpose(H, Ht, EKF_M, EKF_N);
            mulmat(Pp, Ht, tmp1, EKF_N,EKF_N, EKF_M);
            mulmat(H, Pp, tmp2, EKF_M,EKF_N, EKF_N);
            mulmat(tmp2, Ht, tmp3, EKF_M,EKF_N, EKF_M);
            accum(tmp3, R, EKF_M, EKF_M);
            if (cholsl(tmp3, tmp4, tmp5, EKF_M)) return false;
            mulmat(tmp1, tmp4, G, EKF_N,EKF_M, EKF_M);

            // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
            sub(z, hx, tmp5, EKF_M);
            mulvec(G, tmp5, tmp2, EKF_N, EKF_M);
            add(fx, tmp2, _x, EKF_N);

            // P_k = (I - G_k H_k) P_k
            mulmat(G, H, tmp0, EKF_N,EKF_M, EKF_N);
            negate(tmp0, EKF_N, EKF_N);
            mat_addeye(tmp0, EKF_N);
            mulmat(tmp0, Pp, _P, EKF_N,EKF_N, EKF_N);

            // success
            return true;
        }

    private:

        // Cholesky-decomposition matrix-inversion code, adapated from
        // http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt

        static int choldc1(float * a, float * p, const int n) {

            float sum = 0;

            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    sum = a[i*n+j];
                    for (int k = i - 1; k >= 0; k--) {
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

        static int choldcsl(const float * A, float * a, float * p, const int n) 
        {
            float sum = 0;

            for (int i = 0; i < n; i++) 
                for (int j = 0; j < n; j++) 
                    a[i*n+j] = A[i*n+j];
            if (choldc1(a, p, n)) return 1;
            for (int i = 0; i < n; i++) {
                a[i*n+i] = 1 / p[i];
                for (int j = i + 1; j < n; j++) {
                    sum = 0;
                    for (int k = i; k < j; k++) {
                        sum -= a[j*n+k] * a[k*n+i];
                    }
                    a[j*n+i] = sum / p[j];
                }
            }

            return 0; /* success */
        }


        static int cholsl(const float * A, float * a, float * p, const int n) 
        {
            if (choldcsl(A,a,p,n)) return 1;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    a[i*n+j] = 0.0;
                }
            }
            for (int i = 0; i < n; i++) {
                a[i*n+i] *= a[i*n+i];
                for (int k = i + 1; k < n; k++) {
                    a[i*n+i] += a[k*n+i] * a[k*n+i];
                }
                for (int j = i + 1; j < n; j++) {
                    for (int k = j; k < n; k++) {
                        a[i*n+j] += a[k*n+i] * a[k*n+j];
                    }
                }
            }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    a[i*n+j] = a[j*n+i];
                }
            }

            return 0; /* success */
        }

        // C <- A * B
        static void mulmat(
                const float * a, const float * b, float * c, 
                const int arows, const int acols, const int bcols)
        {
            for (int i=0; i<arows; ++i)
                for (int j=0; j<bcols; ++j) {
                    c[i*bcols+j] = 0;
                    for (int l=0; l<acols; ++l)
                        c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
                }
        }

        static void mulvec(
                const float * a, const float * x, float * y, 
                const int m, const int n)
        {
            for (int i=0; i<m; ++i) {
                y[i] = 0;
                for (int j=0; j<n; ++j)
                    y[i] += x[j] * a[i*n+j];
            }
        }

        static void transpose(
                const float * a, float * at, const int m, const int n)
        {
            for (int i=0; i<m; ++i)
                for (int j=0; j<n; ++j) {
                    at[j*m+i] = a[i*n+j];
                }
        }

        // A <- A + B
        static void accum(float * a, const float * b, const int m, const int n)
        {        
            for (int i=0; i<m; ++i)
                for (int j=0; j<n; ++j)
                    a[i*n+j] += b[i*n+j];
        }

        // C <- A + B
        static void add(const float * a, const float * b, float * c, const int n)
        {
            for (int j=0; j<n; ++j)
                c[j] = a[j] + b[j];
        }

        // C <- A - B
        static void sub(const float * a, const float * b, float * c, const int n)
        {
            for (int j=0; j<n; ++j)
                c[j] = a[j] - b[j];
        }

        static void negate(float * a, const int m, const int n)
        {        
            for (int i=0; i<m; ++i)
                for (int j=0; j<n; ++j)
                    a[i*n+j] = -a[i*n+j];
        }

        static void mat_addeye(float * a, const int n)
        {
            for (int i=0; i<n; ++i)
                a[i*n+i] += 1;
        }

        // State
        float _x[EKF_N];

        // Covariance matrix
        float _P[EKF_N * EKF_N];
};
