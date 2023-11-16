#include <iostream>
#include <stdio.h>
#include <math.h>
#include "rebalance.hh"

extern "C"
{
#include <ftw_param.h>
}; 

float epsilon =  5.96e-08;
float sqrt_epsilon = 0.000244131112315;

int main(int argc, char** argv)
{
    int n_var_points = 9;
    getIntParam((char*)"-n_var_points", &n_var_points);

    float* var_x = new float[n_var_points];
    float* var_y = new float[n_var_points];

    float start_x = -1.0, start_y = 1.0;
    float end_x = 1.0, end_y = 1.0;

    // initialize set of points
    
    printf("%f %f\n", start_x, start_y);
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = start_x + (i + 1) * (end_x - start_x) / (n_var_points + 1);
//        var_y[i] = var_x[i] * var_x[i];
// make the simplest possible linear example y=x, should recover original points
// make parabolic example y=x*x
        var_y[i] = var_x[i] * var_x[i];
        printf("%f %f\n", var_x[i], var_y[i]);
    }
    printf("%f %f\n", end_x, end_y);

    // dump the points before rebalance
/*
    printf("%f %f\n", start_x, start_y);
    for (int i=0; i<n_var_points; i++)
    {
        printf("%f %f\n", var_x[i], var_y[i]);
    }
    printf("%f %f\n", end_x, end_y);
*/

    // call the update to re-space points:
    std::cout << "rebalancing " << std::endl;
    float shrinkage = rebalance_points(start_x, start_y, end_x, end_y, var_x, var_y, n_var_points);

    // dump the points after rebalance
    printf("%f %f\n", start_x, start_y);
    for (int i=0; i<n_var_points; i++) printf("%f %f\n", var_x[i], var_y[i]);
    printf("%f %f\n", end_x, end_y);

//    printf("shrinkage: %f\n", shrinkage);
}
