#include "colors.inc"
background {color Black}

#include "camera_light"
//camera{location<3.0,4.0,5.0> look_at <1,1,1> right 1.0}
//light_source{<90.000000,120.000000,150.000000> color White }

cylinder{ __start__, __end__ 0.0050000 texture{ pigment {color Orange  }  } }

// start and end
sphere{ __start__ , 0.5 texture{ pigment {color Yellow  transmit 0.600000 }  finish {phong 0.800000}  } } 
sphere{ __end__ , 0.5 texture{ pigment {color Yellow  transmit 0.600000 }  finish {phong 0.800000}  } }  

//#include "lj_atoms.pov" // atoms
//#include "lj_voro.pov" // voronoi
#include "lj_var.pov" // var points

