#
# Makefile for TinyEKF GPS example
#
# Copyright (C) 2015 Simon D. Levy
#
# MIT License

SRC = ../../src

all: gps

run: gps
	./gps

gps: gps.c $(SRC)/tinyekf.h
	gcc -Wall -I. -I$(SRC) -o gps gps.c -lm

edit:
	vim gps.c

clean:
	rm -f gps *.o *~ ekf.csv
