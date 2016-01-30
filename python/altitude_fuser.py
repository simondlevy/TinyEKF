#!/usr/bin/env python3
'''
agl_fuser.py - Sonar / Barometer fusion example using TinyEKF.  

We model a single state variable, altitude above sea level (ASL) in centimeters.
This is obtained by fusing the barometer pressure in Pascals and sonar above-ground level
(AGL) in centimters.

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

# above seal-level altitude in centimeters
BASELINE_ASL_CM = 33000

RANGE_CM = 300
LOOPSIZE = 5000

import numpy as np
from tinyekf import EKF
from realtime_plot import RealtimePlotter
from time import sleep
import threading
from math import sin, pi


# ground-truth AGL to sonar measurement, empirically determined:
# see http://diydrones.com/profiles/blogs/altitude-hold-with-mb1242-sonar
def sonarfun(groundtruthAGL):

    return 0.933 * groundtruthAGL - 2.894

# cm to Pascals: see http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
def barofun(groundtruthASL):

    return 101325 * pow((1 - 2.25577e-7 * groundtruthASL), 5.25588)

class AGL_EKF(EKF):

    def __init__(self):

        EKF.__init__(self, 1, 2)

    def f(self, x):

        # State-transition function is identity
        return np.copy(x)

    def getF(self, x):

        # So state-transition Jacobian is identity matrix
        return np.eye(1)

    def h(self, x):

        # Mock-up sonar AGL by using altitude of floor as a baseline
        s = sonarfun(x[0] - BASELINE_ASL_CM)

        # Since state is ASL, we can apply the baro-response function directly to it
        b = barofun(x[0])

        return np.array([b, s])

    def getH(self, x):

        # First derivative of nonlinear baro-measurement function
        # Used http://www.wolframalpha.com
        dpdx = -0.120131 * pow((1 - 2.2577e-7 * x[0]), 4.25588)

        # Sonar response is linear, so derivative is constant
        dsdx = 0.933

        return np.array([[dpdx], [dsdx]])

class AGLPlotter(RealtimePlotter):

    def __init__(self):

        baromin = int(barofun(BASELINE_ASL_CM + RANGE_CM))
        baromax = int(barofun(BASELINE_ASL_CM - RANGE_CM))

        RealtimePlotter.__init__(self, [(BASELINE_ASL_CM,BASELINE_ASL_CM+RANGE_CM), (baromin,baromax), (0,RANGE_CM)], 
                window_name='Altitude Sensor Fusion',
                yticks = [
                    range(int(BASELINE_ASL_CM-RANGE_CM/2), BASELINE_ASL_CM+RANGE_CM, 50),  # Fused
                    range(baromin,baromax,int((baromax-baromin)/10.)),    # Baro
                    range(-20, RANGE_CM, 20)                             # Sonar
                    ],
                styles = ['r','b', 'g'], 
                ylabels=['Fused ASL (cm)', 'Baro (Pa)', 'Sonar AGL (cm)'])

        self.xcurr = 0
        self.fused = 0
        self.count = 0

        self.ekf = AGL_EKF()

    def update(self):

        while True:

            # Model up-and-down motion with a sine wave
            self.count = (self.count + 1) % LOOPSIZE
            climb = 50 * (sin(self.count/float(LOOPSIZE) * 2 * pi) + 1)

            # Simulate a noisy observation of baro and sonar
            self.baro  = barofun(BASELINE_ASL_CM+climb) 
            self.sonar = sonarfun(climb)
     
            # Run the EKF on the current baro and sonar measurements, getting back
            # an updated state estimate made by fusing them.
            # Fused state comes back as an array, so grab first element
            self.fused = self.ekf.step((self.baro, self.sonar))[0]

            plotter.xcurr += 1
            sleep(.001)

    def getValues(self):

        return self.fused, self.baro, self.sonar

if __name__ == '__main__':

    plotter = AGLPlotter()

    thread = threading.Thread(target=plotter.update)
    thread.daemon = True

    thread.start()
    plotter.start()
