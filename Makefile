all: gps_ekf

gps: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.cpp tinyekf.hpp tinyekf.h linalg.h tinyekf.c linalg.c
	g++ -Wall -o gps_ekf gps_ekf.cpp tinyekf.c linalg.c

clean:
	rm -f gps_ekf *.o *~

commit:
	git commit -a --allow-empty-message -m ''
	git push
