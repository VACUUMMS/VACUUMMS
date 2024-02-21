#include "colors.inc"
background {color Black}
//camera{location<3,4,5> look_at <1.5,1.5,1.5> right 1.0}
camera{location<23,24,25> look_at <1.5,1.5,1.5> right 1.0}
light_source{<125.000000,125.000000,125.000000> color White }

// all of the var points
//#include "lj_all.pov"

// all of the atoms
//#include "lj_atoms.pov"

// voronoi cells
//#include "lj_voro.pov"

// all variational paths built up of cylinders
//#include "lj_cylinder.pov"
#include "pmp125_cyl.pov"
