variational: variational.cc
	g++ -g -O0 variational.cc -o variational -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

plot: variational
	./variational -n_var_points 25 -n_iter 25000|tail -27 > x.dat 
	gnuplot -p -e "plot 'x.dat'" 

