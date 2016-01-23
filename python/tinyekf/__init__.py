'''
    TinyEKF in Python

    Copyright (C) 2016 Simon D. Levy

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as 
    published by the Free Software Foundation, either version 3 of the 
    License, or (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
'''

import numpy as np

class Matrix(object):

    def __init__(self, r=0, c=0):

        self.data = np.zeros((r,c)) if r>0 and c>0 else None

    def __str__(self):

        return str(self.data)

    def __mul__(self, other):

        new = Matrix()
        new.data = np.dot(self.data, other.data)
        return new

    def __add__(self, other):

        new = Matrix()
        new.data = self.data + other.data
        return new

    def __sub__(self, other):

        new = Matrix()
        new.data = self.data - other.data
        return new

    def __setitem__(self, key, value):

        self.data[key] = value

    def __getitem__(self, key):

        return self.data[key]

    def transpose(self):

        new = Matrix()
        new.data = self.data.T
        return new

    def invert(self):

        new = Matrix()
        try:
            new.data = np.linalg.inv(self.data)
        except Exception as e:
            print(self.data)
            print(e)
            exit(0)
        return new

    @staticmethod
    def eye(n):

        I = Matrix(n,n)

        for k in range(n):
            I[k,k] = 1

        return I

class Vector(Matrix):

    def __init__(self, n):

        Matrix.__init__(self, n, 1)

    @staticmethod
    def fromTuple(t):

        v = Vector(len(t))

        for k in range(len(t)):
            v[k] = t[k]

        return v



class EKF(object):

    def __init__(self, n, m):
        '''
        Creates an EKF object with n states and m observables.
        '''

        self.x = Vector(n)    # state vector 

        self.P = Matrix(n,n)  # prediction error covariance 
        self.Q = Matrix(n,n)  # process noise covariance 
        self.R = Matrix(m,m)  # measurement error covariance 

        self.G = Matrix(n,m)  # Kalman gain a.k.a. K 

        self.F = Matrix(n,n)  # Jacobian of process model 
        self.H = Matrix(m,n)  # Jacobian of measurement model 

        self.Ht = Matrix(n,m) # transpose of measurement Jacobian 
        self.Ft = Matrix(n,n) # transpose of process Jacobian 
        self.Pp = Matrix(n,n) # P, post-prediction, pre-update 

        self.fx = Vector(n)   # output of user defined f() state-transition function 
        self.hx = Vector(m)   # output of user defined h() measurement function 

        self.I = Matrix.eye(n)

    def setP(self, i, j, value):
        '''
         Sets the value of the prediction error covariance P[i,j].
        '''
        self.P[i,j] = value

    def setQ(self, i, j, value):
        '''
         Sets the value of the process noise covariance Q[i,j].
        '''
        self.Q[i,j] = value

    def setR(self, i, j, value):
        '''
         Sets the value of the observation noise covariance R[i,j].
        '''
        self.R[i,j] = value

    def getX(self, i):
        '''
        Returns the state element at index i.
        '''
        return self.x[i][0]

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
        self.Pp = self.Pp + self.F * self.P * self.F.transpose() + self.Q

        # G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
        self.G = self.Pp * self.H.transpose() * (self.H * self.Pp * self.H.transpose() + self.R).invert()

        # \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
        self.x = self.x + self.G * (Vector.fromTuple(z) - self.hx)

        # P_k = (I - G_k H_k) P_k
        self.P = (self.I - self.G * self.H) * self.Pp
