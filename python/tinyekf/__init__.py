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

    def __init__(self, P):
        '''
        Creates a KF object with an initial covariance matrix P
        '''

        n, _ = P.shape

        # No previous prediction noise covariance
        self.P_pre = None

        # Current state is zero, with diagonal noise covariance matrix
        self.x = np.zeros(n)

        # Set up covariance matrices
        self.P_post = P

        # Identity matrix will be usefel later
        self.eye = np.eye(n)

    def predict(self, fx, F, Q):
        '''
        Runs the prediction step
        '''

        # Predict ----------------------------------------------------

        # $P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}$
        self.P_pre = np.dot(F, self.P_post).dot(F.T) + Q

    def update(self, z, hx, H, R):
        '''
        Runs the udpate step
        '''

        # $G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}$
        G = np.dot(self.P_pre.dot(H.T),
                   np.linalg.inv(H.dot(self.P_pre).dot(H.T) + R))

        # $\hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))$
        self.x = self.x + np.dot(G, (np.array(z) - hx.T).T)

        # $P_k = (I - G_k H_k) P_k$
        self.P_post = np.dot(self.eye - np.dot(G, H), self.P_pre)

    def get(self):

        return self.x
