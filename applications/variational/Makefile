# This was the original, now refactoring to 2D etc
#variational: variational.cc rebalance.o
#	g++ -g -O0 variational.cc rebalance.o -o variational -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

variational_2D.o: variational_2D.cc constants.o
	g++ -g -O0 -c variational_2D.cc constants.o -o variational_2D.o 

variational_3D.o: variational_3D.cc constants.o
	g++ -g -O0 -c variational_3D.cc constants.o -o variational_3D.o 

constants.o: constants.cc
	g++ -g -O0 -c constants.cc 

configuration.o: configuration.hh configuration.cc
	g++ -g -O0 -c configuration.cc 

# Tests

test_configuration: configuration.o
	g++ -g -O0 test_configuration.cc configuration.o -o test_configuration

test_variational_2D: test_variational_2D.cc variational_2D.o constants.o configuration.o
	g++ -g -O0 test_variational_2D.cc variational_2D.o configuration.o constants.o -o test_variational_2D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

test_variational_3D: test_variational_3D.cc variational_3D.o constants.o configuration.o
	g++ -g -O0 test_variational_3D.cc variational_3D.o configuration.o constants.o -o test_variational_3D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

## These tests just check the rebalance routine
#
#test_rebalance_1D: test_rebalance_1D.cc rebalance.o
#	g++ test_rebalance_1D.cc rebalance.o -g -O0 -o test_rebalance_1D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/
#
#test_rebalance_2D: test_rebalance_2D.cc rebalance.o
#	g++ test_rebalance_2D.cc rebalance.o -g -O0 -o test_rebalance_2D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/
#
#test_rebalance_3D: test_rebalance_3D.cc rebalance.o
#	g++ test_rebalance_3D.cc rebalance.o -g -O0 -o test_rebalance_3D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/
#
#test_rebalance: test_rebalance_1D test_rebalance_2D test_rebalance_3D
#	./test_rebalance_1D
#	./test_rebalance_2D
#	./test_rebalance_3D

# Plot tests

#noplot: variational
#	./variational -n_var_points 25 -n_iter 25000
#
#plot: variational
#	./variational -n_var_points 25 -n_iter 25000|tail -27 > x.dat 
#	gnuplot -p -e "plot 'x.dat'" 

2Dplot: test_variational_2D
#	./variational_2D -n_var_points 25 -n_iter 25000|tail -25 > x.dat 
	./test_variational_2D -n_var_points 25 -n_iter 25 | tail -29 | head -27 |sed -e "s/#####//" > x.dat 
	gnuplot -p -e "plot 'x.dat'" 

view_docs: variational.md
	typora variational.md

clean:
	rm -f test_rebalance_1D test_rebalance_2D test_rebalance_3D test_configuration test_variational_2D *.o a.out
