'''
Extended Kalman Filter in Python

Copyright (C) 2016 Simon D. Levy

MIT License
'''

import numpy as np


class EKF(object):
    '''
    A simple class for the Extended Kalman Filter, based on the tutorial in
    http://home.wlu.edu/~levys/kalman_tutorial.
    '''

    def __init__(self, n, m, P=1, Q=1, R=1):
        '''
        Creates a KF object with n states, m observables, and specified values
        for prediction noise covariance P, process noise covariance Q,
        and measurement noise covariance R.  You can also pass in your own
        covariance matrix for each of these parameters.
        '''

        # No previous prediction noise covariance
        self.P_pre = None

        # Current state is zero, with diagonal noise covariance matrix
        self.x = np.zeros(n)

        # Set up covariance matrices
        self.P_post = self._covar(n, P)
        self.Q = self._covar(n, Q)
        self.R = self._covar(m, R)

        # Identity matrix will be usefel later
        self.eye = np.eye(n)

    def step(self, z):
        '''
        Runs one step of the EKF 
        '''

        # Predict ----------------------------------------------------

        # $\hat{x}_k = f(\hat{x}_{k-1})$
        self.x, F = self.f(self.x)

        # $P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}$
        self.P_pre = np.dot(F, self.P_post).dot(F.T) + self.Q

        # Update -----------------------------------------------------

        h, H = self.h(self.x)

        # $G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}$
        G = np.dot(self.P_pre.dot(H.T),
                   np.linalg.inv(H.dot(self.P_pre).dot(H.T) + self.R))

        # $\hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))$
        self.x += np.dot(G, (np.array(z) - h.T).T)

        # $P_k = (I - G_k H_k) P_k$
        self.P_post = np.dot(self.eye - np.dot(G, H), self.P_pre)

    def get(self):

        return self.x

    def f(self, x):
        '''
        Your implementing class should override this method for the
        state-transition function f(x).  Your state-transition function should
        return a NumPy array of n elements representing the new state, and a
        nXn NumPy array of elements representing the the Jacobian of the
        function with respect to the new state.  By default this is just the
        identity function np.copy(x), so the Jacobian is just np.eye(len(x)).
        '''
        return np.copy(x), np.eye(len(x))

    def h(self, x):
        '''
        Your implementing class should override this method for the observation
        function h(x), returning a NumPy array of m elements, and a NumPy array
        of m x n elements representing the Jacobian matrix H of the observation
        function with respect to the observation. For example, your function
        might include a component that turns barometric pressure into altitude
        in meters.By default this function is just the identity function
        np.copy(x), so the Jacobian is just np.eye(len(x)).
        '''
        return np.copy(x), np.eye(len(x))

    def _covar(self, siz, arg):

        return arg * np.eye(siz) if np.shape(arg) == () else arg
