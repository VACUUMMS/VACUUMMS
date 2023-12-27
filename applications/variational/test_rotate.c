#include <stdio.h>
#include <math.h>

void rotate(float v_x, float v_y, float v_z, float theta, 
            float u_x, float u_y, float u_z, 
            float *vout_x, float *vout_y, float *vout_z);

float extract(float m_x, float m_y, float m_z, 
              float v_x, float v_y, float v_z, 
              float *u_x, float *u_y, float *u_z);

int main()
{
    float a = sqrt(3)/3;
    float b = sqrt(2);

    float m_x = 0, m_y = 0, m_z = 1;
    // float n_x = a, n_y = a, n_z = a;
    float n_x = 1, n_y = 0, n_z = 0;

    printf("rotating (%f, %f, %f) to (%f, %f, %f):\n", m_x, m_y, m_z, n_x, n_y, n_z);

    float u_x, u_y, u_z;
    float theta = extract(m_x, m_y, m_z, n_x, n_y, n_z, &u_x, &u_y, &u_z);
    printf("extracted theta = %f, u = (%f, %f, %f)\n", theta, u_x, u_y, u_z);

    float vout_x, vout_y, vout_z;

    float v_x = 1, v_y = 0, v_z = 0;
    rotate(v_x, v_y, v_z, theta, u_x, u_y, u_z, &vout_x, &vout_y, &vout_z);
    printf("rotated i vector to: (%f, %f, %f)\n", vout_x, vout_y, vout_z);

    v_x = 0, v_y = 1, v_z = 0;
    rotate(v_x, v_y, v_z, theta, u_x, u_y, u_z, &vout_x, &vout_y, &vout_z);
    printf("rotated j vector: (%f, %f, %f)\n", vout_x, vout_y, vout_z);

    v_x = 0, v_y = 0, v_z = 1;
    rotate(v_x, v_y, v_z, theta, u_x, u_y, u_z, &vout_x, &vout_y, &vout_z);
    printf("rotated k vector: (%f, %f, %f)\n", vout_x, vout_y, vout_z);
}

void rotate(float v_x, float v_y, float v_z, float theta, 
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

float extract(float m_x, float m_y, float m_z, 
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
    printf("Dot product is %f\n", dot);
    float theta = acos(dot);
    printf("Theta is %f\n", theta);

    return theta;
}

