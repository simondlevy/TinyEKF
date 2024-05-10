'''
Extended Kalman Filter in Python

Copyright (C) 2016 Simon D. Levy

MIT License
'''

import numpy as np


class TinyEkf(object):
    '''
    A simple class for the Extended Kalman Filter, based on the tutorial in
    http://home.wlu.edu/~levys/kalman_tutorial.
    '''

    def __init__(self, diag):
        '''
        '''

        self.n = len(diag)

        self.p = np.eye(self.n)
 
        np.fill_diagonal(self.p, diag)

        # Identity matrix will be usefel later
        self.eye = np.eye(self.n)

    def predict(self, xold):
        '''
        '''

        xnew = xold

        F = np.zeros((self.n, self.n))

        return xnew, F

