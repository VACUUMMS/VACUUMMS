#include "colors.inc"
background {color Black}

#include camera_light.pov
//camera{location<3.0,4.0,5.0> look_at <1,1,1> right 1.0}
//light_source{<90.000000,120.000000,150.000000> color White }

// draw the straight line path
// <2.554920, 1.927518, 1.946361>
// <0.622468, 0.276067, 0.880539>

cylinder{<2.554920, 1.927518, 1.946361>, <0.622468, 0.276067, 0.880539>, 0.0050000 texture{ pigment {color Orange  }  } }

// start and end
sphere{<2.554920, 1.927518, 1.946361>, 0.353553390593274 texture{ pigment {color Yellow  transmit 0.600000 }  finish {phong 0.800000}  } } 
sphere{<0.622468, 0.276067, 0.880539>, 0.353553390593274 texture{ pigment {color Yellow  transmit 0.600000 }  finish {phong 0.800000}  } }  

//#include "lj_atoms.pov" // atoms
//#include "lj_voro.pov" // voronoi
#include "lj_var.pov" // var points


