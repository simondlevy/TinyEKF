#!/usr/bin/env python3

'''
    TinyEKF mouse-tracking example.  

    Inspired by

      http://opencvexamples.blogspot.com/2014/01/kalman-filter-implementation-tracking.html

    Copyright (C) 2016 Simon D. Levy

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as 
    published by the Free Software Foundation, either version 3 of the 
    License, or (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
'''

DISPLAY_SIZE       = 600
DISPLAY_BORDER     = 4
CANVAS_MARGIN      = 20
DISPLAY_BACKGROUND = 'black'
MOUSE_COLOR        = 'green'
LINE_WIDTH         = 3

import tkinter as tk
from tinyekf import EKF

class Tracker(tk.Frame):

    def __init__(self):

        tk.Frame.__init__(self, borderwidth = DISPLAY_BORDER, relief = 'sunken')
        self.master.geometry(str(DISPLAY_SIZE)+ "x" + str(DISPLAY_SIZE))
        self.master.title('Mouse Tracker')
        self.grid()
        self.master.rowconfigure(0, weight = 1)
        self.master.columnconfigure(0, weight = 1)
        self.grid(sticky = tk.W+tk.E+tk.N+tk.S)

        canvas_width = DISPLAY_SIZE-2*DISPLAY_BORDER
        canvas_height = DISPLAY_SIZE-2*DISPLAY_BORDER
        
        self.canvas =  tk.Canvas(self, width = canvas_width, height = canvas_height, background = DISPLAY_BACKGROUND)
        self.canvas.grid(row = 0, column = 0)

        self.bind('<Key>', self.handle_key)
        self.canvas.bind('<Motion>', self.handle_motion)

        self.reset_lines()

        self.focus_set()

    def reset_lines(self):

            self.lines = []

            self.x = -1
            self.y = -1

    def out_of_bounds(self, coord, dimname):

            return coord < DISPLAY_BORDER or coord > int(self.canvas[dimname])-DISPLAY_BORDER
 
    def handle_motion(self, event):

        x, y = event.x, event.y

        if self.x != -1:
            if self.out_of_bounds(x, 'width') or self.out_of_bounds(y, 'height'):
                [self.canvas.delete(line) for line in self.lines]
                self.reset_lines()
            else:
                self.lines.append(self.canvas.create_line(self.x, self.y, x, y, fill=MOUSE_COLOR, width=LINE_WIDTH))

        self.x = x
        self.y = y

                 
    def handle_key(self, event):

        # Make sure the frame is receiving input!
        self.focus_force()
        if event.keysym == 'Escape':
            exit(0)

if __name__ == '__main__':
    
    game = Tracker()
    
    game.mainloop()




