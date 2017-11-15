#
# Makefile for TinyEKF GPS example
#
# Copyright (C) 2015 Simon D. Levy
#
# MIT License

all: gps_ekf

run: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.c ../tiny_ekf.c ../tiny_ekf.h
	gcc -Wall -I. -I.. -o gps_ekf gps_ekf.c ../tiny_ekf.c -lm

edit:
	vim gps_ekf.c

clean:
	rm -f gps_ekf *.o *~ ekf.csv

commit:
	git commit -a --allow-empty-message -m ''
	git push
