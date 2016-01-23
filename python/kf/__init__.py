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

    def __init__(self, n, m):
        '''
        Creates a KF object with n states and m observables.
        '''
        self.statePre = Matrix(n, 1)
        self.statePost = Matrix(n, 1)
        self.F = Matrix.eye(n)

        self.Q = Matrix.eye(n)
        self.measurementMatrix = Matrix(m, n)
        self.R = Matrix.eye(m)

        self.errorCovPre = Matrix(n,n)
        self.errorCovPost = Matrix(n,n)
        self.G = Matrix(n, m)

        for j in range(n):
            self.Q[j,j] = 1e-4
            self.errorCovPost[j,j] = 0.1

        for j in range(m):
            self.R[j,j] = 1e-1
            self.measurementMatrix[j,j] = 1

    def step(self, z):

        # Predict

        self.statePre = self.F * self.statePost

        temp1 = self.F * self.errorCovPost

        self.errorCovPre = temp1 * self.F + self.Q

        self.statePre.copyTo(self.statePost)
        self.errorCovPre.copyTo(self.errorCovPost)

        # Update

        temp2 = self.measurementMatrix * self.errorCovPre

        temp3 = temp2 * self.measurementMatrix.transpose() + self.R

        temp4 = temp3.invert() * temp2

        G = temp4.transpose()

        temp5 = Vector.fromTuple(z) - self.measurementMatrix * self.statePre
        
        self.statePost = self.statePre + G * temp5

        self.errorCovPost = self.errorCovPre - G * temp2

    def getX(self, i):
        '''
        Returns the state element at index i.
        '''
        return self.statePost[i][0]

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




