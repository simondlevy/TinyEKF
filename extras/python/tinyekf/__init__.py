'''
Extended Kalman Filter in Python

Copyright (C) 2016 Simon D. Levy

MIT License
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
        self.x = np.zeros(n)
        self.P_post = np.eye(n) * pval

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
        self.x, F = self.f(self.x)

        # $P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}$
        self.P_pre = F * self.P_post * F.T + self.Q

        # Update -----------------------------------------------------

        h, H = self.h(self.x)

        # $G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}$
        G = np.dot(self.P_pre.dot(H.T), np.linalg.inv(H.dot(self.P_pre).dot(H.T) + self.R))

        # $\hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))$
        self.x += np.dot(G, (np.array(z) - h.T).T)

        # $P_k = (I - G_k H_k) P_k$
        self.P_post = np.dot(self.I - np.dot(G, H), self.P_pre)

        # return self.x.asarray()
        return self.x

    @abstractmethod
    def f(self, x):
        '''
        Your implementing class should define this method for the state-transition function f(x).
        Your state-transition fucntion should return a NumPy array of n elements representing the
        new state, and a nXn NumPy array of elements representing the the Jacobian of the function
        with respect to the new state.  Typically this is just the identity
        function np.copy(x), so the Jacobian is just np.eye(len(x)).  '''
        raise NotImplementedError()    

    @abstractmethod
    def h(self, x):
        '''
        Your implementing class should define this method for the observation function h(x), returning
        a NumPy array of m elements, and a NumPy array of m x n elements representing the Jacobian matrix
        H of the observation function with respect to the observation. For
        example, your function might include a component that turns barometric
        pressure into altitude in meters.
        '''
        raise NotImplementedError()    
