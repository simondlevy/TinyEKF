all: gps_ekf

gps: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.cpp
	g++ -o gps_ekf gps_ekf.cpp

commit:
	git commit -a
	git push
