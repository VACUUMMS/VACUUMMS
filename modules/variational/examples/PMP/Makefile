#################

# generate cavity-cavity pairs, all-to-all. The diagram is cluttered 
# but shows how the various paths tend to coalesce to the available passages

# External application that uses VACUUMMS API and library to run the variational calculation
cavity_pairs: cavity_pairs.cc
	g++ cavity_pairs.cc -o cavity_pairs -lvacuumms_cpp -lvacuumms_variational -lftw_general -lm -I${VACUUMMS_ROOT}/include/include -L${VACUUMMS_ROOT}/lib -Wl,-rpath,${VACUUMMS_ROOT}/lib

#################

# Voronoi var pair generator, map voro edges to verts to cavities found by insertion at verts, 
# and back to pairs based on voro edges

# Generate the gfg input file from the lammps data
PMP.gfg: PMP.lmps
	preprocess.sh    

# Extract the list of Voronoi vertices
PMP.vrt: PMP.gfg
	gfg2voro -box 20.862703481463092 20.862703481463092 20.862703481463092 -verts < PMP.gfg > PMP.vrt 

# Use the list of Voronoi vertices as starting points for test particle insertion, 
# and emit index of which vertex/insertion point finds which cavity
PMP.vcv: PMP.gfg PMP.vrt
	vddx -n_threads 16 -n_steps 10000 -box 20.862703481463092 20.862703481463092 20.862703481463092 \
         -vertices_file PMP.vrt \
         < PMP.gfg \
         > PMP.vcv

# Scrub the list of found cavities for duplicates. 
# Generate a list of indices that maps which vertex finds which cavity
PMP.vnq: PMP.vcv
	vuniq -box 20.862703481463092 20.862703481463092 20.862703481463092 < PMP.vcv -mapfile_name PMP.map > PMP.vnq

PMP.map: PMP.vnq

# Generate the list of voronoi edges.
PMP.edg: PMP.gfg
	gfg2voro -box 20.862703481463092 20.862703481463092 20.862703481463092 -edges < PMP.gfg > PMP.edg 

# Map the edges to their corresponding vertices, the vertices to the cavities they find. 
# Emit a list of the resulting pairs of cavity centers to be used as start/end points 
# for variational calculation of pore network. 
PMP.var: PMP.map PMP.vnq PMP.edg PMP.vrt
	edges2var -box 20.862703481463092 20.862703481463092 20.862703481463092 \
              -edges_file PMP.edg \
              -verts_file PMP.vrt \
              -map_file PMP.map \
              < PMP.vnq \
              > PMP.var

# External application that uses VACUUMMS API and library to run the variational calculation
voro_pair: voro_pair.cc
	g++ voro_pair.cc -o voro_pair -lvacuumms_cpp -lvacuumms_variational -lftw_general -lm \
            -I${VACUUMMS_ROOT}/include/include \
            -L${VACUUMMS_ROOT}/lib \
            -Wl,-rpath,${VACUUMMS_ROOT}/lib

# Run the variational calculation on each of the pairs
voro_pair.out: voro_pair PMP.var
	./voro_pair -box 20.862703481463092 20.862703481463092 20.862703481463092 \
                -delta_max 0.1 \
                -n_iter 225 \
                -n_var_points 45 \
                < PMP.var \
                | tee voro_pair.out

# Filter the output to generate a POVRAY scene file to display the generated variational paths
voro_pair.pov: voro_pair.out
	grep "###CYLINDER" < voro_pair.out | sed -e "s/###CYLINDER//" > voro_pair.pov

cavities.pov: PMP.vnq
	awk '{print $$2"\t"$$3"\t"$$4"\t"$$5}' < PMP.vnq | cav2pov -no_header -transmit 0.7 -phong 0.7 > cavities.pov
#       | awk '{print $$1"\t"$$2}' \
#	    | cav2pov -transmit 0.7 -phong 0.7 \
#        > cavities.pov

# Run POVRAY to generate the image, based on the includes in the base.pov file. 
base.png:
	povray -W3840 -H2160 base.pov   	 

twist: twist.c
	gcc twist.c -lm -o twist

inside-mp4: twist
	rm -rf pov
	mkdir pov
	./twist
	pwd
	./twist.sh
	ffmpeg -f image2 -framerate 30 -i pov/PMP-%03d.png -vcodec libx264 -crf 22 -pix_fmt yuv420p PMP_inside_HD.mp4

