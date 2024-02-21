
pmp: pmp.cc
	g++ pmp.cc -o pmp -lvacuumms_cpp -lvacuumms_variational -lftw_general -lm -I${VACUUMMS_ROOT}/include/include -L${VACUUMMS_ROOT}/lib -Wl,-rpath,${VACUUMMS_ROOT}/lib

pmp.out: pmp
	./pmp > pmp125.out 

PMP.cav: 
	pddx -box 20.862703481463092 20.862703481463092 20.862703481463092 -n_samples 1000000 -n_threads 16 < PMP.gfg > PMP.cav 
