# Generate a simulated LJ fluid configuration at density of 1 and temperature of 1
ljx.gfg:
	ljx -N 32 -box 3.174802103936399 3.174802103936399 3.174802103936399 -end_mcs 10000 |tail -34 | head -32 | awk '{print $3"\t"$4"\t"$5"\t1.0\t1.0"}' > ljx.gfg 

lj_atoms.pov: ljx.gfg
	gfg2pov -no_header < ljx.gfg -transmit 0.7 -phong 1 > lj_atoms.pov 

lj_voro.pov: ljx.gfg
	gfg2voro -povray -box 3.174802103936399 3.174802103936399 3.174802103936399 < ljx.gfg > lj_voro.pov
	sed -ie "s/,r/, 0.005000 texture{ pigment {color White } }/" lj_voro.pov

#lj_all: ljx.gfg lj_all.cc ../../variational_3D.o ../../constants.o ../../configuration.o ../../quaternion.o
#	g++ lj_all.cc -o lj_all ../../variational_3D.o ../../configuration.o ../../constants.o ../../quaternion.o -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/ -I../..

lj_all: lj_all.cc
	g++ lj_all.cc -o lj_all -lvacuumms_variational -lftw_general -lm -I${VACUUMMS_HOME}/include -L${VACUUMMS_HOME}/lib

lj_all.out: lj_all ljx.gfg
	./lj_all -n_iter 12500 -n_var_points 125 -box 3.174802103936399 3.174802103936399 3.174802103936399 > lj_all.out

# intermediate in .cav form
lj_all.pov: lj_all.out
	grep POV < lj_all.out | awk '{print $$2"\t"$$3"\t"$$4"\t0.01"}' | cav2pov -no_header > lj_all.pov 

lj_cylinder.pov: lj_all.out
	grep CYLINDER < lj_all.out | sed -e "s/###CYLINDER\s//" > lj_cylinder.pov

lj.png: lj_atoms.pov lj_voro.pov lj_all.pov lj_cylinder.pov
	povray -W3840 -H2160 lj.pov

inside-mp4:
	rm -rf pov
	mkdir pov
	ffmpeg -f image2 -framerate 30 -i pov/lj-%03d.png -vcodec libx264 -crf 22 -pix_fmt yuv420p lj_all_inside_HD.mp4

clean:
	rm -f lj_all.pov lj_atoms.pov lj_voro.pov lj_cylinder.pov lj.png *.pove

veryclean: clean
	rm -f lj_all lj_all.out ljx.gfg

