all: gps_ekf

run: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.cpp ../TinyEKF.h ../TinyEKF.cpp
	g++ -Wall -I.. -o gps_ekf gps_ekf.cpp ../TinyEKF.cpp

edit:
	vim gps_ekf.cpp

clean:
	rm -f gps_ekf *.o *~ ekf.csv
