# Scripts for manipulation of data and generation of images

## run_subsitutions.sh

Takes the original lammps file and atom data, renders it in the .gpg format used as input to vacuumms, and moves (crams) the atoms to set the simulation box origin to (0.0, 0.0, 0.0).

## run_cesa.sh

Run the Cavity Energetic Sizing Algorithm (CESA) to find cavities in the sample. Roughly a million samples are generated. 

## run_uniq.sh

Examine the cavities found by CESA and remove duplicates.

## run_cav2cluster.sh

Look for overlapping sets (clusters) of cavities from the set of unique cavities generated in the last step. Most clusters contain a single cavity.

## run_metrics.sh

Analyze the cavities and cavity clusters and extract statistics for size distribution, radius of gyration, end-to-end distance, surface area, and volume distribution.

## run_generate_sdl.sh

Generate POVRay Scene Definition Language (SDL) to describe the image of atoms and cavities, then render the images using POVRay. 

## run_gfg2fvi.sh

(requires NVidia GPU system) Generate discretized map of the free volume index FVI of the sample. 

## run_fvi2tiff.sh

Translate the discrete map of FVI to a TIFF file, to be used as input to Paraview, ML, etc. 
