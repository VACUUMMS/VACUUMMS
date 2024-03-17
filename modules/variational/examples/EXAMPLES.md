# VARIATIONAL MODULE EXAMPLES 

**The variational module enables calculation of paths of least energy between two points. Instead of a single test particle insertion, an entire curve of discrete points is inserted into the material of interest. The initial condition is a straight line spanning the start and end points. Successive iterations involve perturbing each of the points of the curve individually, then updating the entire set of points concurrently before performing subsequent iterations. Examples are provided for several models, including a simulated sample of PMP. **

## RANDOM

Simplest example. A single variational trajectiory is mapped and iterated in a small, randomly ordered configuration of particles. 

## LJ atoms

A small configuration (32 atoms) of a Lennard-Jones fluid is generated, and a map of variational trajectories is produced by removing pairs of atoms and using their coordinates as endpoints for the variational trajectories. A set of scripts is provided which will generate POVRAY SDL and render frames and video perspective of a full 360 degree rotation.

## PMP

See Makefile for usage and invocation. Tools used include: vddx, vuniq, edges2var.

### cavity_pairs.cc 

External application that uses VACUUMMS API and library to run the variational calculation. Generates cavity-cavity pairs, all-to-all. The generated diagram is cluttered but shows how the various paths tend to coalesce to the available passages

### preprocess.sh    

Script to preprocess data from lammps frame to gfg format for VACUUMMS

### vddx 

Use the list of Voronoi vertices as starting points for test particle insertion and emit index of which vertex/insertion point finds which cavity

### vuniq 

Scrub the list of found cavities for duplicates. Generate a list of indices that maps which vertex finds which cavity

### edges2var 

Map the edges to their corresponding vertices, the vertices to the cavities they find.  Emit a list of the resulting pairs of cavity centers to be used as start/end points for variational calculation of pore network. 

### voro_pair.cc

External application that uses VACUUMMS API and library to run the variational calculation Voronoi var pair generator, map voro edges to verts to cavities found by insertion at verts, and back to pairs based on voro edges

