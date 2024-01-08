
variational_2D.o: variational_2D.cc constants.o
	g++ -g -O0 -c variational_2D.cc constants.o -o variational_2D.o 

variational_3D.o: variational_3D.cc constants.o quaternion.o
	g++ -g -O0 -c variational_3D.cc constants.o -o variational_3D.o 

constants.o: constants.cc
	g++ -g -O0 -c constants.cc 

configuration.o: configuration.hh configuration.cc
	g++ -g -O0 -c configuration.cc 

quaternion.o: quaternion.hh quaternion.cc
	g++ -g -O0 -c quaternion.cc 

# Tests

test_variational_FCC: test_variational_FCC.cc variational_3D.o constants.o configuration.o quaternion.o
	g++ -g -O0 test_variational_FCC.cc -o test_variational_FCC variational_3D.o configuration.o constants.o quaternion.o -o test_variational_FCC -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/
test_configuration: configuration.o
	g++ -g -O0 test_configuration.cc configuration.o -o test_configuration

#test_variational_2D: test_variational_2D.cc variational_2D.o constants.o configuration.o
#	g++ -g -O0 test_variational_2D.cc variational_2D.o configuration.o constants.o -o test_variational_2D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

#test_variational_3D: test_variational_3D.cc variational_3D.o constants.o configuration.o
#	g++ -g -O0 test_variational_3D.cc variational_3D.o configuration.o constants.o -o test_variational_3D -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

2Dplot: test_variational_2D
	./test_variational_2D -n_var_points 25 -n_iter 25 | tail -29 | head -27 |sed -e "s/#####//" > x.dat 
	gnuplot -p -e "plot 'x.dat'" 

view_docs: share/doc/variational.md
	typora share/doc/variational.md

clean:
	rm -f test_rebalance_1D test_rebalance_2D test_rebalance_3D test_configuration test_variational_2D *.o a.out

vis:
	./test_variational_FCC -alpha 0.05 -n_iter 2500 -n_var_points 25 |grep "####" |tail -27 | head -26 |sed -e "s/#####//" | awk '{print $$1"\t"$$2"\t"$$3"\t0.05\t0"}' > fcc.vis

pov: vis
	vis2pov < fcc.vis  -camera 0.25 0.25 1.5 > fcc_top.pov 
	sed -ie "s/look_at <0,0,0>/look_at <0.25, 0.25, 0>/" fcc_top.pov
	povray fcc_top.pov 
	vis2pov < fcc.vis  -camera 1.31 -0.81 0.5 > fcc_right.pov 
	sed -ie "s/look_at <0,0,0>/look_at <0.25, 0.25, 0.5>/" fcc_right.pov
	povray fcc_right.pov 

vis-clean:
	rm -f fcc.vis fcc.pov fcc_top.pov fcc_right.pov
