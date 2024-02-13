
variational_2D.o: variational_2D.cc constants.o
	g++ -g -O0 -c variational_2D.cc constants.o -o variational_2D.o 

variational_3D.o: variational_3D.cc constants.o quaternion.o
	g++ -g -O0 -c variational_3D.cc constants.o -o variational_3D.o 

constants.o: constants.cc
	g++ -g -O0 -c constants.cc -I/home/frankwillmore/VACUUMMS/include

configuration.o: configuration.hh configuration.cc
	g++ -g -O0 -c configuration.cc -I/home/frankwillmore/VACUUMMS/include

quaternion.o: quaternion.hh quaternion.cc
	g++ -g -O0 -c quaternion.cc -I/home/frankwillmore/VACUUMMS/include

# Tests

test_variational_FCC: test_variational_FCC.cc variational_3D.o constants.o configuration.o quaternion.o
	g++ -g -O0 test_variational_FCC.cc -o test_variational_FCC variational_3D.o configuration.o constants.o quaternion.o -o test_variational_FCC -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

test_variational_lj: test_variational_lj.cc variational_3D.o constants.o configuration.o quaternion.o
	g++ -g -O0 test_variational_lj.cc -o test_variational_lj variational_3D.o configuration.o constants.o quaternion.o -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

test_variational_random: test_variational_random.cc variational_3D.o configuration.o constants.o quaternion.o
	g++ -g -O0 test_variational_random.cc -o test_variational_random variational_3D.o configuration.o constants.o quaternion.o -o test_variational_random -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/
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

fcc: fcc.gfg test_variational_FCC
	./test_variational_FCC -n_iter 12500 -n_var_points 125 |grep "#####" |tail -127 | sed -e "s/#####//" > fcc_var.out
	awk '{print $$1"\t"$$2"\t"$$3"\t0.01\t1.0"}' < fcc_var.out > fcc_var.gfg
	gfg2pov -no_header -color Green < fcc_var.gfg > fcc_var.pov 
	gfg2pov -no_header < fcc.gfg -transmit 0.6 -phong 1 > fcc_atoms.pov 
	gfg2voro -povray -box 2 2 2 < fcc.gfg > fcc_voro.pov
	sed -ie "s/,r/, 0.0025 texture{ pigment {color Yellow } }/" fcc_voro.pov
	povray -W3840 -H2160 fcc.pov

lj: lj.gfg test_variational_lj
	./test_variational_lj -n_iter 12500 -n_var_points 125 |grep "#####" |tail -127 | sed -e "s/#####//" > lj_var.out
	awk '{print $$1"\t"$$2"\t"$$3"\t0.01\t1.0"}' < lj_var.out > lj_var.gfg
	gfg2pov -no_header -color Green < lj_var.gfg > lj_var.pov 
	gfg2pov -no_header < lj.gfg -transmit 0.7 -phong 1 > lj_atoms.pov 
	gfg2voro -povray -box 2 2 2 < lj.gfg > lj_voro.pov
	sed -ie "s/,r/, 0.0025 texture{ pigment {color Yellow } }/" lj_voro.pov
	povray -W3840 -H2160 lj.pov

ljx.gfg:
	ljx -N 32 -box 3.174802103936399 3.174802103936399 3.174802103936399 -end_mcs 10000 |tail -34 | head -32 | awk '{print $3"\t"$4"\t"$5"\t1.0\t1.0"}' > ljx.gfg 

vis-clean:
	rm -f fcc.pove fcc.png fcc_var.pov fcc_atoms.pov fcc_voro.pov fcc_var.out facc_var.gfg 

veryclean: vis-clean clean

lj_all: lj.gfg lj_all.cc variational_3D.o constants.o configuration.o quaternion.o
	g++ -g -O0 lj_all.cc -o lj_all variational_3D.o configuration.o constants.o quaternion.o -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/

random.gfg:
	python -c "for i in range(10):from random import random; print(str(random()) + '\t' + str(random()) + '\t' + str(random()) +'\t0.25\t1.0')" > random.gfg


