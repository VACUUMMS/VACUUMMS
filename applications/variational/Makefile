variational: variational.cc rebalance.o
	g++ -g -O0 variational.cc rebalance.o -o variational -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

tester: tester.cc rebalance.o
	g++ -g -O0 tester.cc rebalance.o -o tester -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

runtester: tester
	./tester

rebalance.o: rebalance.cc
	g++ -g -O0 -c rebalance.cc 

noplot: variational
	./variational -n_var_points 25 -n_iter 25000

plot: variational
	./variational -n_var_points 25 -n_iter 25000|tail -27 > x.dat 
	gnuplot -p -e "plot 'x.dat'" 

