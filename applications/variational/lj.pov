#include "colors.inc"
background {color Black}
camera{location<3.0,2.0,3.0> look_at <0,0,0> right 1.0}
light_source{<90.000000,120.000000,150.000000> color White }

// draw the straight line path
// cylinder{<0.000000, 1.000000, 1.000000>, <1.500000, 1.500000, 1.000000>, 0.0050000 texture{ pigment {color Orange  }  } }
cylinder{<0.317301, 0.825600, 0.805090>, <1.691121, 1.815592, 1.190474>, 0.0050000 texture{ pigment {color Orange  }  } }
//  6 0.317301    0.825600    0.805090    1.0 1.0
// 24 1.691121    1.815592    1.190474    1.0 1.0

// start and end
sphere{<0.317301, 0.825600, 0.805090>, 0.250000 texture{ pigment {color Yellow  transmit 0.600000 }  finish {phong 0.800000}  } } 
sphere{<1.691121, 1.815592, 1.190474>, 0.250000 texture{ pigment {color Yellow  transmit 0.600000 }  finish {phong 0.800000}  } }  

//#include "lj_atoms.pov" // atoms
//#include "lj_voro.pov" // voronoi
#include "lj_var.pov" // var points
