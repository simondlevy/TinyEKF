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

        self.x = np.zeros((1, n))

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
        Performs one step of the prediction and update based on observations in tuple z
        '''
        None

