#!/usr/bin/env python3
'''
altitude_fuser.py - Sonar / Barometer fusion example using TinyEKF.  

We model a single state variable, altitude above sea level (ASL) in centimeters.
This is obtained by fusing the barometer pressure in Pascals and sonar above-ground level
(ASL) in centimters.

Also requires RealtimePlotter: https://github.com/simondlevy/RealtimePlotter

Copyright (C) 2016 Simon D. Levy

This code is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this code. If not, see <http://www.gnu.org/licenses/>.
'''

# for plotting
BASELINE_ASL_CM = 33000
RANGE_CM = 300

import numpy as np
from tinyekf import EKF
from realtime_plot import RealtimePlotter
from time import sleep
import threading
from math import sin, pi

class ASL_EKF(EKF):

    def __init__(self):

        # One state (ASL), two measurements (baro, sonar), with larger-than-usual
        # measurement covariance noise to help with sonar blips.
        EKF.__init__(self, 1, 2, rval=.5)

    def f(self, x):

        # State-transition function is identity
        return np.copy(x)

    def getF(self, x):

        # So state-transition Jacobian is identity matrix
        return np.eye(1)

    def h(self, x):

        # State value is ASL
        asl = x[0]

        # Convert ASL cm to sonar AGL cm by subtracting off ASL baseline
        s = self.sonarfun(asl - self.getBaselineASL())

        # Convert ASL cm to Pascals: see http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
        b = self.barofun(asl)

        return np.array([b, s])

    def getH(self, x):

        # First derivative of nonlinear baro-measurement function
        # Used http://www.wolframalpha.com
        dpdx = -0.120131 * pow((1 - 2.2577e-7 * x[0]), 4.25588)

        # Sonar response is linear, so derivative is constant
        dsdx = 0.933

        return np.array([[dpdx], [dsdx]])

    # ground-truth AGL to sonar measurement, empirically determined:
    # see http://diydrones.com/profiles/blogs/altitude-hold-with-mb1242-sonar
    def sonarfun(self, agl):

        return 0.933 * agl - 2.894

    # Convert ASL cm to Pascals: see http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
    def barofun(self, asl):

        return 101325 * pow((1 - 2.25577e-7 * asl), 5.25588)


class ASL_Plotter(RealtimePlotter):

    def __init__(self, ekf):

        self.ekf = ekf

        baromin = int(self.ekf.barofun(BASELINE_ASL_CM + RANGE_CM))
        baromax = int(self.ekf.barofun(BASELINE_ASL_CM - RANGE_CM))

        RealtimePlotter.__init__(self, [(BASELINE_ASL_CM,BASELINE_ASL_CM+RANGE_CM), (baromin,baromax), (0,RANGE_CM)], 
                window_name='Altitude Sensor Fusion',
                yticks = [
                    range(int(BASELINE_ASL_CM-RANGE_CM/2), BASELINE_ASL_CM+RANGE_CM, 50),  # Fused
                    range(baromin,baromax,int((baromax-baromin)/10.)),    # Baro
                    range(-20, RANGE_CM, 20)                             # Sonar
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

            plotter.xcurr += 1
            sleep(.001)

    def getValues(self):

        return self.fused, self.baro, self.sonar

# Simulation ===============================================================================

class _Sim_ASL_EKF(ASL_EKF):

    def __init__(self):

        ASL_EKF.__init__(self)

    def getBaselineASL(self):

        return BASELINE_ASL_CM

class _Sim_ASLPlotter(ASL_Plotter):

    def __init__(self):

        ASL_Plotter.__init__(self, _Sim_ASL_EKF())
        self.count = 0

    def getSensors(self):

        LOOPSIZE = 5000

        # Model up-and-down motion with a sine wave
        self.count = (self.count + 1) % LOOPSIZE
        climb = 50 * (sin(self.count/float(LOOPSIZE) * 2 * pi) + 1)

        baro  = self.ekf.barofun(BASELINE_ASL_CM+climb) 

        # Add noise to simulated sonar at random intervals
        climb += 50 if np.random.rand()>0.9 else 0
        sonar = self.ekf.sonarfun(climb)

        return baro, sonar

if __name__ == '__main__':

    plotter = _Sim_ASLPlotter()

    thread = threading.Thread(target=plotter.update)
    thread.daemon = True

    thread.start()
    plotter.start()
