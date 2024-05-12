/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

typedef struct {

    int n;          
    int m;          

    float x[EKF_N];    

    float P[EKF_N][EKF_N];  
    float Q[EKF_N][EKF_N];  
    float R[EKF_M][EKF_M];  

    float G[EKF_N][EKF_M];  

    float F[EKF_N][EKF_N];  
    float H[EKF_M][EKF_N];  

    float Ht[EKF_N][EKF_M]; 
    float Ft[EKF_N][EKF_N]; 
    float Pp[EKF_N][EKF_N]; 

    float fx[EKF_N];   
    float hx[EKF_M];   

    
    float tmp0[EKF_N][EKF_N];
    float tmp1[EKF_N][EKF_M];
    float tmp2[EKF_M][EKF_N];
    float tmp3[EKF_M][EKF_M];
    float tmp4[EKF_M][EKF_M];
    float tmp5[EKF_M]; 

} ekf1_t;        

class TinyEKF {

    public:

        /**
         * Returns the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @return state value at index
         */
        float getX(int i) 
        { 
            return this->ekf.x[i]; 
        }

        /**
         * Sets the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @param value value to set
         */
        void setX(int i, float value) 
        { 
            this->ekf.x[i] = value; 
        }

        /**
          Performs one step of the prediction and update.
         * @param z observation vector, length <i>m</i>
         * @return true on success, false on failure caused by
         * non-positive-definite matrix.
         */
        bool step(float * z) 
        { 
            this->model(this->ekf.fx, this->ekf.F, this->ekf.hx, this->ekf.H); 

            return ekf_step(&this->ekf, z) ? false : true;
        }
    private:

        ekf1_t ekf;

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

        static void zeros(float * a, int m, int n)
        {
            int j;
            for (j=0; j<m*n; ++j)
                a[j] = 0;
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

        typedef struct {

            float * x;    /* state vector */

            float * P;  /* prediction error covariance */
            float * Q;  /* process noise covariance */
            float * R;  /* measurement error covariance */

            float * G;  /* Kalman gain; a.k.a. K */

            float * F;  /* Jacobian of process model */
            float * H;  /* Jacobian of measurement model */

            float * Ht; /* transpose of measurement Jacobian */
            float * Ft; /* transpose of process Jacobian */
            float * Pp; /* P, post-prediction, pre-update */

            float * fx;  /* output of user defined f() state-transition function */
            float * hx;  /* output of user defined h() measurement function */

            /* temporary storage */
            float * tmp0;
            float * tmp1;
            float * tmp2;
            float * tmp3;
            float * tmp4;
            float * tmp5; 

        } ekf2_t;

        static void unpack(void * v, ekf2_t * ekf, int n, int m)
        {
            /* skip over n, m in data structure */
            char * cptr = (char *)v;
            cptr += 2*sizeof(int);

            float * dptr = (float *)cptr;
            ekf->x = dptr;
            dptr += n;
            ekf->P = dptr;
            dptr += n*n;
            ekf->Q = dptr;
            dptr += n*n;
            ekf->R = dptr;
            dptr += m*m;
            ekf->G = dptr;
            dptr += n*m;
            ekf->F = dptr;
            dptr += n*n;
            ekf->H = dptr;
            dptr += m*n;
            ekf->Ht = dptr;
            dptr += n*m;
            ekf->Ft = dptr;
            dptr += n*n;
            ekf->Pp = dptr;
            dptr += n*n;
            ekf->fx = dptr;
            dptr += n;
            ekf->hx = dptr;
            dptr += m;
            ekf->tmp0 = dptr;
            dptr += n*n;
            ekf->tmp1 = dptr;
            dptr += n*m;
            ekf->tmp2 = dptr;
            dptr += m*n;
            ekf->tmp3 = dptr;
            dptr += m*m;
            ekf->tmp4 = dptr;
            dptr += m*m;
            ekf->tmp5 = dptr;
        }

        void ekf_init(void * v, int n, int m)
        {
            /* retrieve n, m and set them in incoming data structure */
            int * ptr = (int *)v;
            *ptr = n;
            ptr++;
            *ptr = m;

            /* unpack rest of incoming structure for initlization */
            ekf2_t ekf;
            unpack(v, &ekf, n, m);

            /* zero-out matrices */
            zeros(ekf.P, n, n);
            zeros(ekf.Q, n, n);
            zeros(ekf.R, m, m);
            zeros(ekf.G, n, m);
            zeros(ekf.F, n, n);
            zeros(ekf.H, m, n);
        }

        int ekf_step(void * v, float * z)
        {        
            /* unpack incoming structure */

            int * ptr = (int *)v;
            int n = *ptr;
            ptr++;
            int m = *ptr;

            ekf2_t ekf2;
            unpack(v, &ekf2, n, m); 

            /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
            mulmat(ekf2.F, ekf2.P, ekf2.tmp0, n, n, n);
            transpose(ekf2.F, ekf2.Ft, n, n);
            mulmat(ekf2.tmp0, ekf2.Ft, ekf2.Pp, n, n, n);
            accum(ekf2.Pp, ekf2.Q, n, n);

            /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
            transpose(ekf2.H, ekf2.Ht, m, n);
            mulmat(ekf2.Pp, ekf2.Ht, ekf2.tmp1, n, n, m);
            mulmat(ekf2.H, ekf2.Pp, ekf2.tmp2, m, n, n);
            mulmat(ekf2.tmp2, ekf2.Ht, ekf2.tmp3, m, n, m);
            accum(ekf2.tmp3, ekf2.R, m, m);
            if (cholsl(ekf2.tmp3, ekf2.tmp4, ekf2.tmp5, m)) return 1;
            mulmat(ekf2.tmp1, ekf2.tmp4, ekf2.G, n, m, m);

            /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
            sub(z, ekf2.hx, ekf2.tmp5, m);
            mulvec(ekf2.G, ekf2.tmp5, ekf2.tmp2, n, m);
            add(ekf2.fx, ekf2.tmp2, ekf2.x, n);

            /* P_k = (I - G_k H_k) P_k */
            mulmat(ekf2.G, ekf2.H, ekf2.tmp0, n, m, n);
            negate(ekf2.tmp0, n, n);
            mat_addeye(ekf2.tmp0, n);
            mulmat(ekf2.tmp0, ekf2.Pp, ekf2.P, n, n, n);

            /* success */
            return 0;
        }
    protected:

        /**
         * The current state.
         */
        float * x;

        /**
         * Initializes a TinyEKF object.
         */
        TinyEKF() { 
            ekf_init(&this->ekf, EKF_N, EKF_M); 
            this->x = this->ekf.x; 
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
            this->ekf.P[i][j] = value; 
        }

        /** Sets the specified value of the process noise covariance.
         * <i>Q<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setQ(int i, int j, float value) 
        { 
            this->ekf.Q[i][j] = value; 
        }

        /**
         * Sets the specified value of the observation noise covariance.
         * <i>R<sub>i,j</sub> = value</i> @param i row index
         * @param j column index
         * @param value value to set
         */
        void setR(int i, int j, float value) 
        { 
            this->ekf.R[i][j] = value; 
        }

};
