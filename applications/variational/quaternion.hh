#include <stdio.h>
#include <math.h>

void rotate_vector(float v_x, float v_y, float v_z, float theta, 
            float u_x, float u_y, float u_z, 
            float *vout_x, float *vout_y, float *vout_z);

/* rotate a vector through an angle theta around a given axis 
   IN:  v_x, v_y, v_z - vector to be rotated
        u_x, u_y, u_z - unit axis of rotation
        theta - angle of rotation
   OUT: vout_x, vout_y, vout_z - rotated vector
*/
  
float extract_axis(float m_x, float m_y, float m_z, 
              float n_x, float n_y, float n_z, 
              float *u_x, float *u_y, float *u_z);

/* extract an axis perpendicular to an angle between two vectors
   IN:  m_x, m_y, m_z - first vector
        n_x, n_y, n_z - second vector
   OUT: u_x, u_y, u_z - unit axis of rotation
        return value - theta, angle of rotation (0 <= theta < PI/2)
*/
