#!/usr/bin/env python3
'''
TinyEKF mouse-tracking example.  

Inspired by

  http://opencvexamples.blogspot.com/2014/01/kalman-filter-implementation-tracking.html

Copyright (C) 2016 Simon D. Levy

'''

import tkinter as tk
root = tk.Tk()

def motion(event):
    x, y = event.x, event.y
    print('{}, {}'.format(x, y))

root.bind('<Motion>', motion)
root.mainloop()
