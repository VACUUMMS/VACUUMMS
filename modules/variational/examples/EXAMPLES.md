# VARIATIONAL MODULE EXAMPLES 

**The variational module enables calculation of paths of least energy between two points. Instead of a single test particle insertion, an entire curve of discrete points is inserted into the material of interest. The initial condition is a straight line spanning the start and end points. Successive iterations involve perturbing each of the points of the curve individually, then updating the entire set of points concurrently before performing subsequent iterations. Examples are provided for several models, including a simulated sample of PMP. **

## RANDOM

Simplest example. A single variational trajectiory is mapped and iterated in a small, randomly ordered configuration of particles. 

## LJ atoms

A small configuration (32 atoms) of a Lennard-Jones fluid is generated, and a map of variational trajectories is produced by removing pairs of atoms and using their coordinates as endpoints for the variational trajectories. A script (twist.sh) and code (twist.c) are provided which will generate POVRAY SDL and render frames and video perspective of a full 360 degree rotation.

## PMP

See Makefile for usage and invocation. Tools used include: vddx, vuniq, edges2var. A similar script (twist.sh) and code (twist.c) are included, which allow a movie showing a view from the inside of the sample.
A summary of what is produced here is to take the voronoi vertices of a polymer configuration, use them as insertion points to find cavities, and to generate the connectivity graph for the cavities based on their adjacency in the voronoi graph. Not that not every Voronoi vertex will map to a *unique* cavity. When generating the final graph of variational curves, pairs of vertices are omitted when they either map to the same cavity (degenerate case) or when the pairs of vertices map to an already mapped pair of cavities (duplicate case). Also, pairs are discarded when they are clearly mapped across a periodic boundary, as determined by the distance between centers exceeding half of the box dimension, defined as sqrt(box_x * box_x + box_y * box_y + box_z * box_z).

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

