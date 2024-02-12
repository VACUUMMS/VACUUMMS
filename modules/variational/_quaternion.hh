#include <stdio.h>
#include <math.h>

void rotate_vector(vacuumms_float v_x, vacuumms_float v_y, vacuumms_float v_z, vacuumms_float theta, 
            vacuumms_float u_x, vacuumms_float u_y, vacuumms_float u_z, 
            vacuumms_float *vout_x, vacuumms_float *vout_y, vacuumms_float *vout_z);

/* rotate a vector through an angle theta around a given axis 
   IN:  v_x, v_y, v_z - vector to be rotated
        u_x, u_y, u_z - unit axis of rotation
        theta - angle of rotation
   OUT: vout_x, vout_y, vout_z - rotated vector
*/
  
vacuumms_float extract_axis(vacuumms_float m_x, vacuumms_float m_y, vacuumms_float m_z, 
              vacuumms_float n_x, vacuumms_float n_y, vacuumms_float n_z, 
              vacuumms_float *u_x, vacuumms_float *u_y, vacuumms_float *u_z);

/* extract an axis perpendicular to an angle between two vectors
   IN:  m_x, m_y, m_z - first vector
        n_x, n_y, n_z - second vector
   OUT: u_x, u_y, u_z - unit axis of rotation
        return value - theta, angle of rotation (0 <= theta < PI/2)
*/
