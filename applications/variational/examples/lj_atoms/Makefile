# Generate a simulated LJ fluid configuration at density of 1 and temperature of 1
ljx.gfg:
	ljx -N 32 -box 3.174802103936399 3.174802103936399 3.174802103936399 -end_mcs 10000 |tail -34 | head -32 | awk '{print $3"\t"$4"\t"$5"\t1.0\t1.0"}' > ljx.gfg 

#lj_var.gfg: lj.gfg lj.cav
#	../test_variational_lj -n_iter 12500 -n_var_points 125 -box 3.174802103936399 3.174802103936399 3.174802103936399 -start $$(head -1 lj.cav) -end $$(tail -1 lj.cav) | grep "#####" |tail -127 | sed -e "s/#####//" | awk '{print $$1"\t"$$2"\t"$$3"\t0.01\t1.0"}' > lj_var.gfg

#lj_var.pov: lj_var.gfg
#	gfg2pov -no_header -color Green < lj_var.gfg > lj_var.pov 

lj_atoms.pov: lj.gfg
	gfg2pov -no_header < ljx.gfg -transmit 0.7 -phong 1 > lj_atoms.pov 

lj_voro.pov: lj.gfg
	gfg2voro -povray -box 2 2 2 < lj.gfg > lj_voro.pov
	sed -ie "s/,r/, 0.0025 texture{ pigment {color Yellow } }/" lj_voro.pov

lj_all: ljx.gfg lj_all.cc ../../variational_3D.o ../../constants.o ../../configuration.o ../../quaternion.o
	g++ -g -O0 lj_all.cc -o lj_all ../../variational_3D.o ../../configuration.o ../../constants.o ../../quaternion.o -lm -I/home/frankwillmore/VACUUMMS/include /home/frankwillmore/VACUUMMS/cmake/libraries/libftw_general.so -Wl,-rpath,/home/frankwillmore/VACUUMMS/cmake/libraries/ -I../..

lj_all.out: lj_all
	./lj_all > lj_all.out

# intermediate in .cav form
lj_all.pov: lj_all.out
	grep POV < lj_all.out | awk '{print $$2"\t"$$3"\t"$$4"\t0.01"}' > lj_all.int
#	cav2pov -camera 3 4 5 -look_at 1.5 1.5 1.5 < lj_all_var.int > lj_all_var.pov 
	cav2pov -no_header < lj_all.int > lj_all.pov 
	rm lj_all.int

# this is in repo
#lj.pov: lj_all_var.pov 

lj.png: lj.pov
	povray -W3840 -H2160 lj.pov

inside-mp4:
	ffmpeg -f image2 -framerate 30 -i pov/lj-%03d.png -vcodec libx264 -crf 22 -pix_fmt yuv420p lj_all_inside_HD.mp4

clean:
	rm -f lj_all.int lj.png lj_all.out lj_all.pov lj.png

#veryclean: clean
#	rm -f lj_all lj_all.out ljx.gfg

