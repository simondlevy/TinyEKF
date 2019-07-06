#!/usr/bin/env python3
'''
Plots output of Arduino EKF sketch.

Copyright (C) 2015 Simon D. Levy

MIT License
'''

ARDUINO_PORT = '/dev/ttyACM0'
ARDUINO_BAUD = 9600

from serial import Serial
from realtime_plot import RealtimePlotter
import numpy as np
from time import sleep

class EKF_Plotter(RealtimePlotter):

    def __init__(self):

        self.port = Serial(ARDUINO_PORT, ARDUINO_BAUD)

        p_lo = 970
        p_hi = 990
        p_range = p_lo, p_hi

        t_lo = 20
        t_hi = 40 
        t_range = t_lo, t_hi

        t_ticks = tuple(range(t_lo,t_hi,2))
        p_ticks = tuple(range(p_lo,p_hi,2))

        RealtimePlotter.__init__(self, [p_range, t_range, t_range, p_range, t_range],
                window_name='EKF demo',
                yticks = [p_ticks, t_ticks, t_ticks, p_ticks, t_ticks],
                styles = ['r--', 'b-', 'g-', 'k-', 'k-'], 
                ylabels=['Baro Press (mb)','Baro Temp(C)', 'LM35 Temp(C)', 'Baro EKF', 'Temp EKF'],
                interval_msec=1)

        self.pbaro, self.tbaro, self.tlm35, self.p_ekf, self.t_ekf = 0,0,0,0,0
        self.msg = ''

    def getValues(self):

       return self.pbaro, self.tbaro, self.tlm35, self.p_ekf, self.t_ekf

def _update(plotter):

    while True:

        try:
            c = plotter.port.read(1).decode('utf-8')

            if c == '\n':
                try:
                    plotter.pbaro, plotter.tbaro, plotter.tlm35, plotter.p_ekf, plotter.t_ekf = \
                        map(lambda s:float(s), plotter.msg.split())
                except:
                    None
                plotter.msg = ''
            else:
                plotter.msg += c
        except:
            pass

        # yield to other thread
        sleep(0)

if __name__ == '__main__':

    import threading

    plotter = EKF_Plotter()

    thread = threading.Thread(target=_update, args = (plotter,))
    thread.daemon = True
    thread.start()

    plotter.start()
 
