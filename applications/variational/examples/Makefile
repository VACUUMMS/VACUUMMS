
####### LJ ####### ####### ####### ####### ####### #######

ljx.gfg:
	ljx -N 32 -box 3.174802103936399 3.174802103936399 3.174802103936399 -end_mcs 10000 |tail -34 | head -32 | awk '{print $3"\t"$4"\t"$5"\t1.0\t1.0"}' > ljx.gfg 

# remove two atoms and get the variational trace between them
# use records 8 and 9
lj.gfg:
	head -3 < ljx.gfg > lj.gfg
	tail -27 < ljx.gfg >> lj.gfg
#	head -7 < ljx.gfg > lj.gfg
#	tail -23 < ljx.gfg >> lj.gfg

lj.cav:
	head -5 < ljx.gfg | tail -2 | awk '{print $$1"\t"$$2"\t"$$3}' > lj.cav
#	head -9 < ljx.gfg | tail -2 | awk '{print $$1"\t"$$2"\t"$$3}' > lj.cav

lj_var.gfg: lj.gfg lj.cav
	../test_variational_lj -n_iter 12500 -n_var_points 125 -box 3.174802103936399 3.174802103936399 3.174802103936399 -start $$(head -1 lj.cav) -end $$(tail -1 lj.cav) | grep "#####" |tail -127 | sed -e "s/#####//" | awk '{print $$1"\t"$$2"\t"$$3"\t0.01\t1.0"}' > lj_var.gfg

lj_var.pov: lj_var.gfg
	gfg2pov -no_header -color Green < lj_var.gfg > lj_var.pov 

lj_atoms.pov: lj.gfg
	gfg2pov -no_header < lj.gfg -transmit 0.7 -phong 1 > lj_atoms.pov 

lj_voro.pov: lj.gfg
	gfg2voro -povray -box 2 2 2 < lj.gfg > lj_voro.pov
	sed -ie "s/,r/, 0.0025 texture{ pigment {color Yellow } }/" lj_voro.pov

lj.pov: lj_var.pov lj_atoms.pov lj_voro.pov 
	cp _lj.pov lj.pov
	sed -ie "s/__start__/$$(head -1 lj.cav | awk '{print "<"$$1", "$$2", "$$3">"}')/g" lj.pov
	sed -ie "s/__end__/$$(tail -1 lj.cav | awk '{print "<"$$1", "$$2", "$$3">"}')/g" lj.pov
	
lj.png: lj.pov
	povray -W3840 -H2160 lj.pov

####### end LJ ####### ####### ####### ####### ####### #######

vis-clean:
	rm -f lj.gfg lj.cav lj_var.gfg lj_atoms.pov lj_voro.pov lj.pov lj.png

clean:
	rm -f *.o a.out

veryclean: vis-clean clean


	
	
	

#lj: lj.gfg test_variational_lj
#    ./test_variational_lj -n_iter 12500 -n_var_points 125 |grep "#####" |tail -127 | sed -e "s/#####//" > lj_var.out
#    awk '{print $$1"\t"$$2"\t"$$3"\t0.01\t1.0"}' < lj_var.out > lj_var.gfg
#    gfg2pov -no_header -color Green < lj_var.gfg > lj_var.pov 
#    gfg2pov -no_header < lj.gfg -transmit 0.7 -phong 1 > lj_atoms.pov 
#    gfg2voro -povray -box 2 2 2 < lj.gfg > lj_voro.pov
#    sed -ie "s/,r/, 0.0025 texture{ pigment {color Yellow } }/" lj_voro.pov
#    povray -W3840 -H2160 lj.pov

