'''
    Linear KalmanFilter in Python

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

class KF(object):

    def __init__(self, n, m, pval=0.1, qval=1e-4, rval=0.1):
        '''
        Creates a KF object with n states and m observables.
        '''

        # No previous prediction noise covariance
        self.P_pre = None

        # Current state is zero, with diagonal noise covariance matrix
        self.x = Vector(n)
        self.P_post = Matrix.eye(n) * pval

        # State transition function is identity
        self.F = Matrix.eye(n)

        self.Q = Matrix.eye(n) * qval
        self.H = Matrix.eye(m, n)
        self.R = Matrix.eye(m) * rval

    def step(self, z):
        '''
        Runs one step of the EKF on observations z, where z is a tuple of length M.
        '''

        # Predict ----------------------------------------------------

        # $\hat{x}_k = f(\hat{x}_{k-1})$
        self.x = Vector.fromData(self.f(self.x.data))

        # $P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}$
        self.P_pre = self.F * self.P_post * self.F.transpose() + self.Q

        self.P_post = self.P_pre.copy()

        # Update -----------------------------------------------------

        # $G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}$
        G = self.P_pre * self.H.transpose() * (self.H * self.P_pre * self.H.transpose() + self.R).invert()

        # $\hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))$
        #self.x = self.x + G * (Vector.fromTuple(z) - Vector.fromData(self.h(self.x.data)))
        self.x = self.x + G * (Vector.fromTuple(z) - Vector.fromData(self.h(self.x.data)))

        # $P_k = (I - G_k H_k) P_k$
        self.P_post = self.P_pre - G * (self.H * self.P_pre)

        return self.x.asarray()

# Linear Algebra support =============================================

class Matrix(object):

    def __init__(self, r=0, c=0):

        self.data = np.zeros((r,c)) if r>0 and c>0 else None

    def __str__(self):

        return str(self.data) + " " + str(self.data.shape)

    def __mul__(self, other):

        new = Matrix()

        if type(other).__name__ in ['float', 'int']:
            new.data = np.copy(self.data)
            new.data *= other
        else:
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

    def asarray(self):

        return np.asarray(self.data[:,0])

    def copy(self):

        new = Matrix()
        new.data = np.copy(self.data)
        return new

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
    def eye(n, m=0):

        I = Matrix()

        if m == 0:
            m = n

        I.data = np.eye(n, m)

        return I

class Vector(Matrix):

    def __init__(self, n=0):

        self.data = np.zeros((n,1)) if n>0 else None

    @staticmethod
    def fromTuple(t):

        v = Vector(len(t))

        for k in range(len(t)):
            v[k] = t[k]

        return v


    @staticmethod
    def fromData(data):

        v = Vector()

        v.data = data

        return v





