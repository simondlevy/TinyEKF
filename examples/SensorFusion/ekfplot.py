#!/usr/bin/env python
'''
Real-time plot demo using sine waves.

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

from realtime_plot import RealtimePlotter
import numpy as np

# Simple example with threading

class _SinePlotter(RealtimePlotter):

    def __init__(self):

        RealtimePlotter.__init__(self, [(-1,+1), (-1,+1)], 
                window_name='Sinewave demo',
                yticks = [(-1,0,+1),(-1,0,+1)],
                styles = ['r--', 'b-'], 
                ylabels=['Slow', 'Fast'])

        self.xcurr = 0

    def getValues(self):

        return self._getRow(1), self._getRow(2)

    def _getRow(self, row):

        size = len(self.x)
        
        return np.sin(row*2*np.pi*(float(self.xcurr)%size)/size)


def _update(plotter):

    from time import sleep

    while True:

        plotter.xcurr += 1
        sleep(.002)

if __name__ == '__main__':

    import threading

    plotter = _SinePlotter()

    thread = threading.Thread(target=_update, args = (plotter,))
    thread.daemon = True
    thread.start()

    plotter.start()
 
