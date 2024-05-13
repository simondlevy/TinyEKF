#!/usr/bin/env python3

'''
kalman_mousetracker.py - OpenCV mouse-tracking demo using TinyEKF

Adapted from

   http://www.morethantechnical.com/2011/06/17/
     simple-kalman-filter-for-tracking-using-opencv-2-2-w-code/

Copyright (C) 2016 Simon D. Levy

MIT License
'''

import cv2
import numpy as np
from sys import exit

from tinyekf import EKF

# This delay will affect the Kalman update rate
DELAY_MSEC = 100

# Arbitrary display params
WINDOW_NAME = 'Kalman Mousetracker [ESC to quit]'
WINDOW_SIZE = 500


class MouseInfo(object):
    '''
    A class to store X,Y points
    '''

    def __init__(self):

        self.x, self.y = -1, -1

    def __str__(self):

        return '%4d %4d' % (self.x, self.y)


def mouseCallback(event, x, y, flags, mouse_info):
    '''
    Callback to update a MouseInfo object with new X,Y coordinates
    '''

    mouse_info.x = x
    mouse_info.y = y


def drawCross(img, center, r, g, b):
    '''
    Draws a cross a the specified X,Y coordinates with color RGB
    '''

    d = 5
    t = 2

    color = (r, g, b)

    ctrx = center[0]
    ctry = center[1]

    cv2.line(img, (ctrx - d, ctry - d), (ctrx + d, ctry + d), color, t,
             cv2.LINE_AA)
    cv2.line(img, (ctrx + d, ctry - d), (ctrx - d, ctry + d), color, t,
             cv2.LINE_AA)


def drawLines(img, points, r, g, b):
    '''
    Draws lines
    '''

    cv2.polylines(img, [np.int32(points)], isClosed=False, color=(r, g, b))


def newImage():
    '''
    Returns a new image
    '''

    return np.zeros((500, 500, 3), dtype=np.uint8)


if __name__ == '__main__':

    # Create a new image in a named window
    img = newImage()
    cv2.namedWindow(WINDOW_NAME)

    # Create an X,Y mouse info object and set the window's mouse callback to
    # modify it
    mouse_info = MouseInfo()
    cv2.setMouseCallback(WINDOW_NAME, mouseCallback, mouse_info)

    # Loop until mouse inside window
    while True:

        if mouse_info.x > 0 and mouse_info.y > 0:
            break

        cv2.imshow(WINDOW_NAME, img)
        if cv2.waitKey(1) == 27:
            exit(0)

    # These will get the trajectories for mouse location and Kalman estiamte
    measured_points = []
    kalman_points = []

    P = np.eye(2) * 1e-1

    Q = np.eye(2) * 1e-4

    R = np.eye(2) * 1e-1

    # Create a new Kalman filter for mouse tracking
    ekf = EKF(P)

    # Loop till user hits escape
    while True:

        # Serve up a fresh image
        img = newImage()

        # Grab current mouse position and add it to the trajectory
        measured = (mouse_info.x, mouse_info.y)
        measured_points.append(measured)

        # Update the Kalman filter with the mouse point, getting the estimate.

        # State-transition function f is identity
        fx = ekf.get()

        # So first derivative of f is identity matrix
        F = np.eye(2)

        # Observation is mouse coordinates
        z = mouse_info.x, mouse_info.y

        # Observiation function h is also identity
        hx = ekf.get()

        # So first deriviative of H is also identity matrix
        H = np.eye(2)

        ekf.predict(fx, F, Q)

        ekf.update(z, hx, H, R)

        estimate = ekf.get()

        # Add the estimate to the trajectory
        estimated = [int(c) for c in estimate]
        kalman_points.append(estimated)

        # Display the trajectories and current points
        drawLines(img, kalman_points,   0,   255, 0)
        drawCross(img, estimated,       255, 255, 255)
        drawLines(img, measured_points, 255, 255, 0)
        drawCross(img, measured, 0,   0,   255)

        # Delay for specified interval, quitting on ESC
        cv2.imshow(WINDOW_NAME, img)
        if cv2.waitKey(DELAY_MSEC) & 0xFF == 27:
            break
