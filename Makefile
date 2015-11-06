all: gps_ekf

gps: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.cpp TinyEKF.h TinyEKF.cpp
	g++ -Wall -o gps_ekf gps_ekf.cpp TinyEKF.cpp

clean:
	rm -f gps_ekf *.o *~

commit:
	git commit -a --allow-empty-message -m ''
	git push
