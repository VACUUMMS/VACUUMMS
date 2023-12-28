#include <math.h>

void rotate_vector(float v_x, float v_y, float v_z, float theta, 
                   float u_x, float u_y, float u_z, 
                   float *vout_x, float *vout_y, float *vout_z)
{
    float a = cos(theta/2);
    float b = u_x * sin(theta/2);
    float c = u_y * sin(theta/2);
    float d = u_z * sin(theta/2);

    *vout_x = (a*a + b*b - c*c - d*d) * v_x + 2*(b*c - a*d) * v_y + 2*(b*d + a*c) * v_z;
    *vout_y = 2*(b*c + a*d) * v_x + (a*a - b*b + c*c - d*d) * v_y + 2*(c*d - a*b) * v_z;
    *vout_z = 2*(b*d - a*c) * v_x + 2*(c*d + a*b) * v_y + (a*a - b*b - c*c + d*d) * v_z;
}

float extract_axis(float m_x, float m_y, float m_z, 
                   float n_x, float n_y, float n_z, 
                   float *u_x, float *u_y, float *u_z)
{
    // u = m X n / |m X n|:
    *u_x = (m_y * n_z) - (m_z * n_y);
    *u_y = (m_z * n_x) - (m_x * n_z);
    *u_z = (m_x * n_y) - (m_y * n_x);

    float magnitude = sqrt((*u_x)*(*u_x) + (*u_y)*(*u_y) + (*u_z)*(*u_z));

    *u_x /= magnitude;
    *u_y /= magnitude;
    *u_z /= magnitude;

    float dot = m_x * n_x + m_y * n_y + m_z * n_z;
//    printf("Dot product is %f\n", dot);
    float theta = acos(dot);
//    printf("Theta is %f\n", theta);

    return theta;
}
