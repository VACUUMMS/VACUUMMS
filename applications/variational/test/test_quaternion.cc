#include "quaternion.hh"

int main()
{
    float a = sqrt(3)/3;
    float b = sqrt(2);

    float m_x = 0, m_y = 0, m_z = 1;
    // float n_x = a, n_y = a, n_z = a;
    float n_x = 1, n_y = 0, n_z = 0;

    printf("rotating (%f, %f, %f) to (%f, %f, %f):\n", m_x, m_y, m_z, n_x, n_y, n_z);

    float u_x, u_y, u_z;
    float theta = extract_axis(m_x, m_y, m_z, n_x, n_y, n_z, &u_x, &u_y, &u_z);
    printf("extracted theta = %f, u = (%f, %f, %f)\n", theta, u_x, u_y, u_z);

    float vout_x, vout_y, vout_z;

    float v_x = 1, v_y = 0, v_z = 0;
    rotate_vector(v_x, v_y, v_z, theta, u_x, u_y, u_z, &vout_x, &vout_y, &vout_z);
    printf("rotated i vector to: (%f, %f, %f)\n", vout_x, vout_y, vout_z);
}
