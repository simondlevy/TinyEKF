'''
    Extended Kalman Filter in Python

    Copyright (C) 2016 Simon D. Levy

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as 
    published by the Free Software Foundation, either version 3 of the 
    License, or (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARANTY without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
'''

import numpy as np
from abc import ABCMeta, abstractmethod

class EKF(object):
    '''
    A abstrat class for the Extended Kalman Filter, based on the tutorial in
    http://home.wlu.edu/~levys/kalman_tutorial.
    '''
    __metaclass__ = ABCMeta

    def __init__(self, n, m, pval=0.1, qval=1e-4, rval=0.1):
        '''
        Creates a KF object with n states, m observables, and specified values for 
        prediction noise covariance pval, process noise covariance qval, and 
        measurement noise covariance rval.
        '''

        # No previous prediction noise covariance
        self.P_pre = None

        # Current state is zero, with diagonal noise covariance matrix
        self.x = np.zeros((n,1))
        self.P_post = np.eye(n) * pval

        # Get state transition and measurement Jacobians from implementing class
        self.F = self.getF(self.x)
        self.H = self.getH(self.x)

        # Set up covariance matrices for process noise and measurement noise
        self.Q = np.eye(n) * qval
        self.R = np.eye(m) * rval
 
        # Identity matrix will be usefel later
        self.I = np.eye(n)

    def step(self, z):
        '''
        Runs one step of the EKF on observations z, where z is a tuple of length M.
        Returns a NumPy array representing the updated state.
        '''

        # Predict ----------------------------------------------------

        # $\hat{x}_k = f(\hat{x}_{k-1})$
        self.x = self.f(self.x)

        # $P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}$
        self.P_pre = self.F * self.P_post * self.F.T + self.Q

        # Update -----------------------------------------------------

        # $G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}$
        G = np.dot(self.P_pre * self.H.T, np.linalg.inv(self.H * self.P_pre * self.H.T + self.R))

        # $\hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))$
        self.x += np.dot(G, (np.array(z) - self.h(self.x).T).T)

        # $P_k = (I - G_k H_k) P_k$
        self.P_post = np.dot(self.I - np.dot(G, self.H), self.P_pre)

        # return self.x.asarray()
        return self.x

    @abstractmethod
    def f(self, x):
        '''
        Your implementing class should define this method for the state transition function f(x),
        returning a NumPy array of n elements.  Typically this is just the identity function np.copy(x).
        '''
        raise NotImplementedError()    

    @abstractmethod
    def getF(self, x):
        '''
        Your implementing class should define this method for returning the n x n Jacobian matrix F of the 
        state transition function as a NumPy array.  Typically this is just the identity matrix np.eye(n).
        '''
        raise NotImplementedError()    

    @abstractmethod
    def h(self, x):
        '''
        Your implementing class should define this method for the observation function h(x), returning
        a NumPy array of m elements. For example, your function might include a component that
        turns barometric pressure into altitude in meters.
        '''
        raise NotImplementedError()    

    @abstractmethod
    def getH(self, x):
        '''
        Your implementing class should define this method for returning the m x n Jacobian matirx H of the 
        observation function as a NumPy array.
        '''
        raise NotImplementedError()    

