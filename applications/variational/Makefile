
variational_2D.o: variational_2D.cc constants.o
	g++ -g -O0 -c variational_2D.cc constants.o -o variational_2D.o 

variational_3D.o: variational_3D.cc constants.o
	g++ -g -O0 -c variational_3D.cc constants.o -o variational_3D.o 

constants.o: constants.cc
	g++ -g -O0 -c constants.cc 

configuration.o: configuration.hh configuration.cc
	g++ -g -O0 -c configuration.cc 

quaternion.o: quaternion.hh quaternion.cc
	g++ -g -O0 -c quaternion.cc 

# Tests

test_configuration: configuration.o
	g++ -g -O0 test_configuration.cc configuration.o -o test_configuration

test_variational_2D: test_variational_2D.cc variational_2D.o constants.o configuration.o
	g++ -g -O0 test_variational_2D.cc variational_2D.o configuration.o constants.o -o test_variational_2D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

test_variational_3D: test_variational_3D.cc variational_3D.o constants.o configuration.o
	g++ -g -O0 test_variational_3D.cc variational_3D.o configuration.o constants.o -o test_variational_3D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

2Dplot: test_variational_2D
	./test_variational_2D -n_var_points 25 -n_iter 25 | tail -29 | head -27 |sed -e "s/#####//" > x.dat 
	gnuplot -p -e "plot 'x.dat'" 

view_docs: share/doc/variational.md
	typora share/doc/variational.md

clean:
	rm -f test_rebalance_1D test_rebalance_2D test_rebalance_3D test_configuration test_variational_2D *.o a.out
