#include "colors.inc"
background {color Black}
camera{location<3.0,2.0,3.0> look_at <0,0,0> right 1.0}
light_source{<90.000000,120.000000,150.000000> color White }

// draw the straight line path
cylinder{<0.000000, 1.000000, 1.000000>, <1.500000, 1.500000, 1.000000>, 0.0050000 texture{ pigment {color Orange  }  } }

// start and end
sphere{<0.000000, 1.000000, 1.000000>, 0.250000 texture{ pigment {color Yellow  transmit 0.600000 }  finish {phong 0.800000}  } } 
sphere{<1.500000, 1.500000, 1.000000>, 0.250000 texture{ pigment {color Yellow  transmit 0.600000 }  finish {phong 0.800000}  } }  

#include "fcc_atoms.pov" // atoms
#include "fcc_voro.pov" // voronoi
#include "fcc_var.pov" // var points
