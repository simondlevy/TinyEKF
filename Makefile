all: gps_ekf

gps: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.cpp tinyekf.hpp tinyekf.cpp
	g++ -Wall -o gps_ekf gps_ekf.cpp tinyekf.cpp

clean:
	rm -f gps_ekf *.o *~

commit:
	git commit -a --allow-empty-message -m ''
	git push
