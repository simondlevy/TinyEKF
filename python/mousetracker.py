#!/usr/bin/env python3

DISPLAY_SIZE       = 600
DISPLAY_BORDER     = 4
CANVAS_MARGIN      = 20
DISPLAY_BACKGROUND = 'black'
MOUSE_COLOR        = 'green'

import tkinter as tk

class RoleGame(tk.Frame):

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

        self.bind('<Key>', self._handle_key)
        self.canvas.bind('<Motion>', self._handle_motion)

        self._reset_lines()

        self.focus_set()

    def _reset_lines(self):

            self.lines = []

            self.x = -1
            self.y = -1

    def _handle_motion(self, event):

        x, y = event.x, event.y

        if self.x != -1:
            if x < 10 or x > int(self.canvas['width'])-10 or y < 10 or y > int(self.canvas['height'])-10:
                [self.canvas.delete(line) for line in self.lines]
                self._reset_lines()
            else:
                self.lines.append(self.canvas.create_line(self.x, self.y, x, y, fill=MOUSE_COLOR))

        self.x = x
        self.y = y

                 
    def _handle_key(self, event):

        # Make sure the frame is receiving input!
        self.focus_force()
        if event.keysym == 'Escape':
            exit(0)

if __name__ == '__main__':
    
    game = RoleGame()
    
    game.mainloop()




