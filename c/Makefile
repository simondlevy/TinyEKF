#
# Makefile for TinyEKF GPS example
#
# Copyright (C) 2015 Simon D. Levy
#
# MIT License

all: gps_ekf

run: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.c tinyekf.c tinyekf.h
	gcc -Wall -I. -DEKF_N=8 -DEKF_M=4 -o gps_ekf gps_ekf.c tinyekf.c -lm

edit:
	vim gps_ekf.c

clean:
	rm -f gps_ekf *.o *~ ekf.csv
