#!/usr/bin/env python3
'''
Python distutils setup file for installing TinyEKF.

Copyright (C) 2016 Simon D. Levy

MIT License
'''

from distutils.core import setup

setup (name = 'TinyEKF',
       packages = ['tinyekf'],
       requires = ['numpy'],
       version = '0.1',
       description = 'A simple Extended Kalman Filter in Python',
       author_email='simon.d.levy@gmail.com',
       url='https://github.com/simondlevy/TinyEKF',
       license='LGPL',
       platforms='Linux; Windows; OS X')

