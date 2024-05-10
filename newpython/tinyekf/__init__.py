'''
Extended Kalman Filter in Python

Copyright (C) 2016 Simon D. Levy

MIT License
'''

import numpy as np
from abc import ABC, abstractmethod


class TinyEkf(ABC):
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

        self.x = np.zeros(self.n)

        self.p = np.eye(self.n)

        np.fill_diagonal(self.p, diag)

        # Identity matrix will be usefel later
        self.eye = np.eye(self.n)

        self.predictionIntervalMsec = predictionIntervalMsec

        self.lastProcessNoiseUpdateMsec = nowMsec

        self.lastPredictionMsec = nowMsec

        self.min_covariance = minCovariance
        self.max_covariance = maxCovariance

        self.isUpdated = False
        self.nextPredictionMsec = 0

    def predict(self, xold, nowMsec=0):
        '''
        '''

        if nowMsec > self.nextPredictionMsec:
            self.nextPredictionMsec = nowMsec + self.predictionIntervalMsec

        if nowMsec >= self.nextPredictionMsec:

            self.isUpdated = True

            shouldAddProcessNoise = nowMsec == 0 or (
                    nowMsec - self.lastProcessNoiseUpdateMsec > 0)

            dt = (nowMsec - self.lastPredictionMsec) / 1000

            xnew, F = self.get_prediction(xold, dt, shouldAddProcessNoise)

            self._multiplyCovariance(F)

            self._cleanupCovariance()

            if shouldAddProcessNoise:

                self.lastProcessNoiseUpdateMsec = nowMsec

                self.x = xnew

        xnew = xold

        return xnew, F

    def update(self, h, error, stdMeasureNoise):

        pass

    def _multiplyCovariance(self, a):

        self.p = np.dot(np.dot(a, self.p), a.transpose())

    def _cleanupCovariance(self):

        for i in range(self.n):

            for j in range(self.n):

                pval = (self.p[i][j] + self.p[j][i]) / 2

                self.p[i][j] = self.p[j][i] = (
                    self.max_covariance
                    if pval > self.max_covariance
                    else self.min_covariance
                    if i == j and pval < self.min_covariance
                    else pval)

    @abstractmethod
    def get_prediction(self, xold, shouldAddProcessNoise):

        pass
