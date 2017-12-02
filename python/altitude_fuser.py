#!/usr/bin/env python3
'''
altitude_fuser.py - Sonar / Barometer fusion example using TinyEKF.  

We model a single state variable, altitude above sea level (ASL) in centimeters.
This is obtained by fusing the barometer pressure in Pascals and sonar above-ground level
(ASL) in centimters.

Also requires RealtimePlotter: https://github.com/simondlevy/RealtimePlotter

Copyright (C) 2016 Simon D. Levy

MIT License
'''

# for plotting
BARO_RANGE    = 20
SONAR_RANGE   = 200
BARO_BASELINE = 97420

import numpy as np
from tinyekf import EKF
from realtime_plot import RealtimePlotter
from time import sleep
import threading
from math import sin, pi

# ground-truth AGL to sonar measurement, empirically determined:
# see http://diydrones.com/profiles/blogs/altitude-hold-with-mb1242-sonar
def sonarfun( agl):

    return 0.933 * agl - 2.894

# Convert ASL cm to Pascals: see http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
def asl2baro( asl):

    return 101325 * pow((1 - 2.25577e-7 * asl), 5.25588)

# Convert Pascals to cm ASL
def baro2asl( pa):

    return (1.0 - pow(pa/ 101325.0, 0.190295)) * 4433000.0


class ASL_EKF(EKF):
    '''
    An abstract class for fusing baro and sonar.  
    '''

    def __init__(self):

        # One state (ASL), two measurements (baro, sonar), with larger-than-usual
        # measurement covariance noise to help with sonar blips.
        EKF.__init__(self, 1, 2, rval=.5)

    def f(self, x):

        # State-transition function is identity
        return np.copy(x), np.eye(1)


    def h(self, x):

        # State value is ASL
        asl = x[0]

        # Convert ASL cm to sonar AGL cm by subtracting off ASL baseline from baro
        s = sonarfun(asl - baro2asl(BARO_BASELINE))

        # Convert ASL cm to Pascals: see http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
        b = asl2baro(asl)

        h = np.array([b, s])

        # First derivative of nonlinear baro-measurement function
        # Used http://www.wolframalpha.com
        dpdx = -0.120131 * pow((1 - 2.2577e-7 * x[0]), 4.25588)

        # Sonar response is linear, so derivative is constant
        dsdx = 0.933

        H = np.array([[dpdx], [dsdx]])

        return h, H


class ASL_Plotter(RealtimePlotter):
    '''
    An abstract class plotting Above Sea Level altitude.
    Implementing class should define getSensors(self), returning baro and sonar readings.
    '''

    def __init__(self):

        self.ekf = ASL_EKF()

        self.baro = 0
        self.sonar = 0

        baromin = BARO_BASELINE - BARO_RANGE
        baromax = BARO_BASELINE + BARO_RANGE

        max_asl_cm      = int(baro2asl(baromin))
        min_asl_cm      = int(baro2asl(baromax))

        RealtimePlotter.__init__(self, [(min_asl_cm,max_asl_cm), (baromin,baromax), (0,SONAR_RANGE)], 
                window_name='Altitude Sensor Fusion',
                yticks = [
                    range(min_asl_cm, max_asl_cm, 50),  # Fused
                    range(baromin,baromax,int((baromax-baromin)/10.)),    # Baro
                    range(0, SONAR_RANGE, 20)                             # Sonar
                    ],
                styles = ['r','b', 'g'], 
                ylabels=['Fused ASL (cm)', 'Baro (Pa)', 'Sonar ASL (cm)'])

        self.xcurr = 0
        self.fused = 0

    def update(self):

        while True:

            self.baro, self.sonar = self.getSensors()

            # Run the EKF on the current baro and sonar measurements, getting back
            # an updated state estimate made by fusing them.
            # Fused state comes back as an array, so grab first element
            self.fused = self.ekf.step((self.baro, self.sonar))[0]

            self.xcurr += 1
            sleep(.001)

    def getValues(self):

        return self.fused, self.baro, self.sonar

# Simulation ===============================================================================


class _Sim_ASLPlotter(ASL_Plotter):

    def __init__(self):

        ASL_Plotter.__init__(self)
        self.count = 0

    def getSensors(self):

        LOOPSIZE = 5000

        # Model up-and-down motion with a sine wave
        self.count = (self.count + 1) % LOOPSIZE
        sine = sin(self.count/float(LOOPSIZE) * 2 * pi)

        baro  = BARO_BASELINE + sine * BARO_RANGE

        # Add noise to simulated sonar at random intervals
        sonar = sonarfun(50*(1-sine)) + (50 if np.random.rand()>0.9 else 0)

        return baro, sonar

if __name__ == '__main__':

    plotter = _Sim_ASLPlotter()

    thread = threading.Thread(target=plotter.update)
    thread.daemon = True

    thread.start()
    plotter.start()
