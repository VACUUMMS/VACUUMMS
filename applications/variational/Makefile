variational: variational.cc rebalance.o
	g++ -g -O0 variational.cc rebalance.o -o variational -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

rebalance.o: rebalance.cc
	g++ -g -O0 -c rebalance.cc -lm -I/home/frankwillmore/VACUUMMS/include 

# Tests

test_variational: test_variational.cc rebalance.o
	g++ -g -O0 tester.cc rebalance.o -o tester -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

test_rebalance_1D: test_rebalance_1D.cc rebalance.o
	g++ test_rebalance_1D.cc rebalance.o -g -O0 -o test_rebalance_1D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

test_rebalance_2D: test_rebalance_2D.cc rebalance.o
	g++ test_rebalance_2D.cc rebalance.o -g -O0 -o test_rebalance_2D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

test_rebalance_3D: test_rebalance_3D.cc rebalance.o
	g++ test_rebalance_3D.cc rebalance.o -g -O0 -o test_rebalance_3D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

test_rebalance: test_rebalance_1D test_rebalance_2D test_rebalance_3D
	./test_rebalance_1D
	./test_rebalance_2D
	./test_rebalance_3D

# Plot tests

#noplot: variational
#	./variational -n_var_points 25 -n_iter 25000
#
#plot: variational
#	./variational -n_var_points 25 -n_iter 25000|tail -27 > x.dat 
#	gnuplot -p -e "plot 'x.dat'" 

clean:
	rm -f test_rebalance_1D test_rebalance_2D test_rebalance_3D rebalance.o variational
