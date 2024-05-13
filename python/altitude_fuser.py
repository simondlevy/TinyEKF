#!/usr/bin/env python3
'''
altitude_fuser.py - Sonar / Barometer fusion example using TinyEKF.

We model a single state variable, altitude above sea level (ASL) in
centimeters.  This is obtained by fusing the barometer pressure in Pascals and
sonar above-ground level (ASL) in centimters.

Copyright (C) 2016 Simon D. Levy

MIT License
'''

import numpy as np
import matplotlib.pyplot as plt
from tinyekf import EKF
from math import sin, pi

# for plotting
BARO_RANGE = 20
SONAR_RANGE = 200
BARO_BASELINE = 97420


# ground-truth AGL to sonar measurement, empirically determined:
# see http://diydrones.com/profiles/blogs/altitude-hold-with-mb1242-sonar
def sonarfun(agl):

    return 0.933 * agl - 2.894


# Convert ASL cm to Pascals: see
# http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
def asl2baro(asl):

    return 101325 * pow((1 - 2.25577e-7 * asl), 5.25588)


# Convert Pascals to cm ASL
def baro2asl(pa):

    return (1.0 - pow(pa / 101325.0, 0.190295)) * 4433000.0


if __name__ == '__main__':

    LOOPSIZE = 100

    count = 0

    N = 100

    # One state (ASL), two measurements (baro, sonar), with
    # larger-than-usual measurement covariance noise to help with sonar
    # blips.
    P = np.eye(1) * 1e-1
    Q = np.eye(1) * 1e-4
    R = np.eye(2) * 5e-1

    ekf = EKF(P)

    baro = np.zeros(N)
    sonar = np.zeros(N)
    fused = np.zeros(N)

    for k in range(N):

        # Model up-and-down motion with a sine wave
        sine = sin((k % LOOPSIZE)/LOOPSIZE * 2 * pi)

        baro[k] = BARO_BASELINE + sine * BARO_RANGE

        # Add noise to simulated sonar at random intervals
        sonar[k] = (sonarfun(50 * (1 - sine)) +
                    (50 if np.random.rand() > 0.9 else 0))

        # Observation is baro and sonar
        z = baro[k], sonar[k]

        # State-transition function f is identity
        fx = ekf.get()

        # So first derivative of f is identity matrix
        F = np.eye(1)

        # State value is ASL
        asl = ekf.get()[0]

        # Convert ASL cm to sonar AGL cm by subtracting off ASL baseline from
        # baro
        s = sonarfun(asl - baro2asl(BARO_BASELINE))

        # Convert ASL cm to Pascals: see
        # http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
        b = asl2baro(asl)

        hx = np.array([b, s])

        # First derivative of nonlinear baro-measurement function
        # Used http://www.wolframalpha.com
        dpdx = -0.120131 * pow((1 - 2.2577e-7 * asl), 4.25588)

        # Sonar response is linear, so derivative is constant
        dsdx = 0.933

        H = np.array([[dpdx], [dsdx]])

        ekf.predict(fx, F, Q)

        ekf.update(z, hx, H, R)

        fused[k] = ekf.get()[0]

    plt.subplot(3, 1, 1)
    plt.plot(fused, 'r')
    plt.ylabel('Fused ASL (cm)')

    plt.subplot(3, 1, 2)
    plt.plot(baro, 'b')
    plt.ylabel('Baro (Pa)')

    plt.subplot(3, 1, 3)
    plt.plot(sonar, 'g')
    plt.ylabel('Sonar ASL (cm)')

    plt.show()
