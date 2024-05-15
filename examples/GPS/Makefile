#
# Makefile for TinyEKF GPS example
#
# Copyright (C) 2015 Simon D. Levy
#
# MIT License

SRC = ../../src

all: gps_ekf

run: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.c $(SRC)/tinyekf.c $(SRC)/tinyekf.h ekf_config.h
	gcc -Wall -I. -I$(SRC) -o gps_ekf gps_ekf.c $(SRC)/tinyekf.c -lm

edit:
	vim gps_ekf.c

clean:
	rm -f gps_ekf *.o *~ ekf.csv
