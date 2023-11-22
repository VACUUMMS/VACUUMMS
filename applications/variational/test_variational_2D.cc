#include <iostream>
#include <stdio.h>
#include <math.h>
#include "variational.hh"

extern "C" 
{
#include <ftw_param.h>
}; 

//float points_x[] = {-0.0149969,0.0149969};
//float points_y[] = {1.0,-1.0};
int n_points = 2;
int n_iter = 1;
int n_var_points = 25; 
int update = 1;
//float var_x[n_var_points], var_y[n_var_points];
float *var_x, *var_y;

float start_x = -2; 
float start_y = 0;
float end_x = 2;
float end_y = 0;

//void variational_2D(float start_x, float start_y, int n_iter, int update, int n_var_points, float(*energy)(float x, float y));

/* 
// Calculate energy at this point
float energy(float x, float y)
{
    float total = 0.0;
    for (int i=0; i<n_points; i++)
    {   
        float r_sq = (points_x[i] - x) * (points_x[i] - x) + (points_y[i] - y) * (points_y[i] - y); 
        float r_6 = r_sq * r_sq * r_sq;
        float r_12 = r_6 * r_6;
        total += (1.0/r_12 - 1.0/r_6);
    }   
    return total;    
}
*/

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getIntParam((char*)"-n_iter", &n_iter);
    getIntParam((char*)"-update", &update);
    getIntParam((char*)"-n_var_points", &n_var_points);

// this should instantiate, and then be destructed, right?
//    Variational2D v = Variational2D(start_x, start_y, end_x, end_y, n_iter, update, n_var_points, &energy);

    char filename[] = "x.gfg";
    Configuration c = Configuration(filename);
    Variational2D v = Variational2D(start_x, start_y, end_x, end_y, n_iter, update, n_var_points, c);
    v.printValues();

}


