#!/usr/bin/env python
'''
Plots output of Arduino EKF sketch.

Copyright (C) 2015 Simon D. Levy

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
'''

ARDUINO_PORT = '/dev/ttyACM0'
ARDUINO_BAUD = 9600

from serial import Serial
from realtime_plot import RealtimePlotter
import numpy as np

class EKF_Plotter(RealtimePlotter):

    def __init__(self):

        self.port = Serial(ARDUINO_PORT, ARDUINO_BAUD)

        RealtimePlotter.__init__(self, [(-1,+1), (-1,+1)], 
                window_name='EKF demo',
                yticks = [(-1,0,+1),(-1,0,+1)],
                styles = ['r--', 'b-'], 
                ylabels=['Slow', 'Fast'])

        self.xcurr,self.ycurr = 0,0
        self.msg = ''

    def getValues(self):

        c = self.port.read(1)

        if c == '\n':
            try:
                self.xcurr = float(self.msg)
            except:
                None
            self.msg = ''
        else:
            self.msg += c

        return self.xcurr, self.ycurr

def _update(plotter):

    from time import sleep

    while True:

        sleep(.002)

if __name__ == '__main__':

    import threading

    plotter = EKF_Plotter()

    thread = threading.Thread(target=_update, args = (plotter,))
    thread.daemon = True
    thread.start()

    plotter.start()
 
