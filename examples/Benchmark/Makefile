# Makefile for building, uploading, and maintaining Arduino TinyEKF Benchmark example
# for those who want to avoid the Arduino IDE.
#
# Copyright (C) 2015 Simon D. Levy
#
# MIT License

PORT = ACM0

#all: build Benchmark
all: build

edit:
	vim Benchmark.ino

build: Benchmark.ino
	arduino --verify Benchmark.ino

upload: Benchmark.ino
	arduino --upload Benchmark.ino

probe:
	ls /dev/ttyA*

listen:
	miniterm.py /dev/tty$(PORT) 9600

commit:
	git commit -a --allow-empty-message -m ''
	git push

clean:
	rm -f *~
