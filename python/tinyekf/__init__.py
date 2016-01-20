'''
    TinyEKF in Python

    Copyright (C) 2016 Simon D. Levy

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as 
    published by the Free Software Foundation, either version 3 of the 
    License, or (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
'''

import numpy as np

class EKF(object):

    def __init__(self, n, m):
        '''
        Creates an EKF object with n states and m observables.
        '''

        self.x = np.zeros((n))

        self.fx = np.zeros((n))
        self.hx = np.zeros((n))

        self.P = np.zeros((n, n))
        self.Q = np.zeros((n, n))
        self.R = np.zeros((m, m))
        self.G = np.zeros((n, m))
        self.F = np.zeros((n, n))
        self.H = np.zeros((m, n))

    def setP(self, i, j, value):
        '''
         Sets the value of the prediction error covariance P[i,j].
        '''
        self.P[i][j] = value

    def setQ(self, i, j, value):
        '''
         Sets the value of the process noise covariance Q[i,j].
        '''
        self.Q[i][j] = value

    def setR(self, i, j, value):
        '''
         Sets the value of the observation noise covariance R[i,j].
        '''
        self.R[i][j] = value

    def getX(self, i):
        '''
        Returns the state element at index i.
        '''
        return self.x[i]

    def setX(self, i, value):
        '''
        Sets the state element at index i.
        '''
        self.x[i] = value

    def step(self, z):
        '''
        Performs one step of the prediction and update based on observations in tuple z.
        Calls subclass model() method with following output arguments:
         fx gets output of state-transition function $f(x_{0 .. n-1})$
         F gets $n \times n$ Jacobian of $f(x)$
         hx gets output of observation function $h(x_{0 .. n-1})$
         H gets $m \times n$ Jacobian of $h(x)$
        '''
        self.model(self.x, self.fx, self.F, self.hx, self.H)
        
        # P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}
        #mulmat(ekf.F, ekf.P, ekf.tmp1, n, n, n);
        #transpose(ekf.F, ekf.Ft, n, n);
        #mulmat(ekf.tmp1, ekf.Ft, ekf.Pp, n, n, n);
        #accum(ekf.Pp, ekf.Q, n, n);

        # G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
        #transpose(ekf.H, ekf.Ht, m, n);
        #mulmat(ekf.Pp, ekf.Ht, ekf.tmp1, n, n, m);
        #mulmat(ekf.H, ekf.Pp, ekf.tmp2, m, n, n);
        #mulmat(ekf.tmp2, ekf.Ht, ekf.tmp3, m, n, m);
        #accum(ekf.tmp3, ekf.R, m, m);
        #if (cholsl(ekf.tmp3, ekf.tmp4, ekf.tmp5, m)) return 1;
        #mulmat(ekf.tmp1, ekf.tmp4, ekf.G, n, m, m);

        # \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k
        #sub(z, ekf.hx, ekf.tmp5, m);
        #mulvec(ekf.G, ekf.tmp5, ekf.tmp2, n, m);
        #add(ekf.fx, ekf.tmp2, ekf.x, n);

        # P_k = (I - G_k H_k) P_k
        #mulmat(ekf.G, ekf.H, ekf.tmp1, n, m, n);
        #negate(ekf.tmp1, n, n);
        #mat_addeye(ekf.tmp1, n);
        #mulmat(ekf.tmp1, ekf.Pp, ekf.P, n, n, n);

 
