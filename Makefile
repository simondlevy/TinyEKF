# Makefile for maintaining TinyEKF
#
# Copyright (C) 2015 Simon D. Levy
#
# MIT License

doc:
	doxygen

clean:
	rm -rf Documentation

commit:
	git commit -a --allow-empty-message -m ''
	git push
