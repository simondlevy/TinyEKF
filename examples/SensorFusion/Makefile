# Makefile for building, uploading, and maintaining Arduino SensorFusion example
# for those who want to avoid the Arduino IDE.
#
# Copyright (C) 2015 Simon D. Levy
# 
# MIT License

PORT = ACM0

#all: build SensorFusion
all: build

edit:
	vim SensorFusion.ino

build: SensorFusion.ino
	arduino --verify SensorFusion.ino

upload: SensorFusion.ino
	arduino --upload SensorFusion.ino

probe:
	ls /dev/ttyA*

listen:
	miniterm.py /dev/tty$(PORT) 9600

plot:
	./ekfplot.py

commit:
	git commit -a --allow-empty-message -m ''
	git push

clean:
	rm -f *~
