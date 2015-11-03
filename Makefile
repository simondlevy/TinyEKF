all: gps_ekf fusion

gps: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.cpp tinyekf.hpp tinyekf.cpp
	g++ -Wall -o gps_ekf gps_ekf.cpp tinyekf.cpp

fuse: fusion
	./fusion

fusion: fusion.cpp tinyekf.hpp tinyekf.cpp
	g++ -Wall -o fusion fusion.cpp tinyekf.cpp

clean:
	rm -f gps_ekf *.o *~

commit:
	git commit -a --allow-empty-message -m ''
	git push
