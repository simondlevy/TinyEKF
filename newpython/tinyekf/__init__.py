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

    def __init__(self, diag,
            nowMsec=0,
            predictionIntervalMsec=0,
            lastProcessUpdateNoiseMsec=0,
            lastPredictionMsec=0,
            minCovariance=-np.inf,
            maxCovariance=+np.inf):
        '''
        '''

        self.n = len(diag)

        self.p = np.eye(self.n)
 
        np.fill_diagonal(self.p, diag)

        # Identity matrix will be usefel later
        self.eye = np.eye(self.n)

        self.predictionIntervalMsec = predictionIntervalMsec

        self.lastProcessNoiseUpdateMsec = nowMsec

        self.lastPredictionMsec = nowMsec

        self._min_covariance = minCovariance
        self._max_covariance = maxCovariance

        self.isUpdated = False
        self.nextPredictionMsec = 0

    def predict(self, xold, nowMsec=0):
        '''
        '''

        if nowMsec > self.nextPredictionMsec:
            self.nextPredictionMsec = nowMsec + self.predictionIntervalMsec

        if nowMsec >= self.nextPredictionMsec:

                self.isUpdated = True

                shouldAddProcessNoise = (
                        nowMsec - self.lastProcessNoiseUpdateMsec > 0)


        xnew = xold

        F = np.zeros((self.n, self.n))

        return xnew, F

