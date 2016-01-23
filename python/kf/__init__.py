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
        self.x_Pre = None
        self.x_Post = Matrix(n, 1)

        self.F = Matrix.eye(n)

        self.Q = Matrix.eye(n)
        self.H = Matrix.eye(m, n)
        self.R = Matrix.eye(m) * rval

        self.P_Pre = None
        self.P_Post = Matrix(n,n)

        for j in range(n):
            self.Q[j,j] = qval
            self.P_Post[j,j] = pval

    def step(self, z):

        # Predict

        self.x_Pre = self.F * self.x_Post

        self.P_Pre = self.F * self.P_Post * self.F + self.Q

        self.x_Pre.copyTo(self.x_Post)
        self.P_Pre.copyTo(self.P_Post)

        # Update

        G = ((self.H * self.P_Pre * self.H.transpose() + self.R).invert() * (self.H * self.P_Pre)).transpose()

        self.x_Post = self.x_Pre + G * (Vector.fromTuple(z) - self.H * self.x_Pre)

        self.P_Post = self.P_Pre - G * (self.H * self.P_Pre)

    def getX(self, i):
        '''
        Returns the x element at index i.
        '''
        return self.x_Post[i][0]

# Linear Algebra support =============================================

class Matrix(object):

    def __init__(self, r=0, c=0):

        self.data = np.zeros((r,c)) if r>0 and c>0 else None

    def __str__(self):

        return str(self.data)

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

    def copyTo(self, other):

        other.data = np.copy(self.data)

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

    def __init__(self, n):

        Matrix.__init__(self, n, 1)

    @staticmethod
    def fromTuple(t):

        v = Vector(len(t))

        for k in range(len(t)):
            v[k] = t[k]

        return v




