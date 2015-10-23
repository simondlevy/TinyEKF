all: gps_ekf

gps: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.cpp tinyekf.hpp tinyekf.h tinyekf.c
	g++ -Wall -o gps_ekf gps_ekf.cpp tinyekf.c

clean:
	rm -f gps_ekf *.o *~

commit:
	git commit -a --allow-empty-message -m ''
	git push
