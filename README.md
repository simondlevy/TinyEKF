# TinyEKF: Lightweight C/C++ Extended Kalman Filter with Arduino Sensor Fusion Example

TinyEKF is simple C/++ implementation of the Extended Kalman filter that is general enough to use on different projects.  In order to make it practical for running on Arduino, STM32, and other microcontrollers, it uses static (compile-time) memory allocation.  The examples folder includes both a "pure C" example from the literature, as well as an Arduino example of sensor fusion.

Arduino users can simply install or drag the whole folder into their Arduino libraries folder, then launch the example by going to File/Sketchbook/libraries/TinyEKF/SensorFusion.  As of 13 November 2015 the example is working with some faked-up sensor data, but will shortly support an actual sensor suite (GPS / magnetometer).

To learn about the Extended Kalman Filter and why it is so useful, try this [interactive tutorial](http://home.wlu.edu/~levys/kalman_tutorial/).
