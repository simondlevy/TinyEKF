# TinyEKF: Lightweight C/C++ Extended Kalman Filter with Arduino Sensor Fusion Example

TinyEKF is a simple C/C++ implementation of the Extended Kalman filter that is general enough to use on different 
projects.  In order to make it practical for running on Arduino, STM32, and other microcontrollers, it uses static 
(compile-time) memory allocation (no "new" or "malloc").  The examples folder includes both a "pure C" example 
from the literature, as well as an Arduino example of sensor fusion.  

Arduino users can simply install or drag the whole folder into their Arduino libraries folder. This folder
contains a little sensor fusion example using a [BMP180 barometer](https://www.sparkfun.com/products/11824) and 
[LM35 temperature sensor](http://www.robotshop.com/en/dfrobot-lm35-linear-temperature-sensor.html).
I have run this example on an Arduino Uno and a Teensy 3.2.

To learn about the Extended Kalman Filter and why it is so useful, try 
this [interactive tutorial](http://home.wlu.edu/~levys/kalman_tutorial/).
