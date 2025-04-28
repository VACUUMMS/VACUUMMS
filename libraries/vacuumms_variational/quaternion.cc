/* Quaternion rotation and axis extraction methods used in generating direction
 * for perturbation of variational points
 */

#include <math.h>
#include <vacuumms/types.h>


/* Rotate the vector (v_x, v_y, v_z) through an angle theta about the axis 
 * (u_x, u_y, u_z), returning the rotated vector as (vout_x, vout_y, vout_z).
 */

void rotate_vector(vacuumms_float v_x, vacuumms_float v_y, vacuumms_float v_z, vacuumms_float theta, 
                   vacuumms_float u_x, vacuumms_float u_y, vacuumms_float u_z, 
                   vacuumms_float *vout_x, vacuumms_float *vout_y, vacuumms_float *vout_z)
{
    vacuumms_float a = cos(theta/2);
    vacuumms_float b = u_x * sin(theta/2);
    vacuumms_float c = u_y * sin(theta/2);
    vacuumms_float d = u_z * sin(theta/2);

    *vout_x = (a*a + b*b - c*c - d*d) * v_x + 2*(b*c - a*d) * v_y + 2*(b*d + a*c) * v_z;
    *vout_y = 2*(b*c + a*d) * v_x + (a*a - b*b + c*c - d*d) * v_y + 2*(c*d - a*b) * v_z;
    *vout_z = 2*(b*d - a*c) * v_x + 2*(c*d + a*b) * v_y + (a*a - b*b - c*c + d*d) * v_z;
}


/* extract the vector (u_x, u_y, u_z) which is perpendicular to the input 
 * vectors (m_x, m_y, m_z) and (n_x, n_y, n_z) and return the angle theta 
 * between them.
 */

vacuumms_float extract_axis(vacuumms_float m_x, vacuumms_float m_y, vacuumms_float m_z, 
                            vacuumms_float n_x, vacuumms_float n_y, vacuumms_float n_z, 
                            vacuumms_float *u_x, vacuumms_float *u_y, vacuumms_float *u_z)
{
    // u = m X n / |m X n|:
    *u_x = (m_y * n_z) - (m_z * n_y);
    *u_y = (m_z * n_x) - (m_x * n_z);
    *u_z = (m_x * n_y) - (m_y * n_x);

    vacuumms_float magnitude = sqrt((*u_x)*(*u_x) + (*u_y)*(*u_y) + (*u_z)*(*u_z));

    *u_x /= magnitude;
    *u_y /= magnitude;
    *u_z /= magnitude;

    vacuumms_float dot = m_x * n_x + m_y * n_y + m_z * n_z;
    vacuumms_float theta = acos(dot);

    return theta;
}
