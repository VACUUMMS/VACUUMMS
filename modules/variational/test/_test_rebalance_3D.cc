#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vacuumms/variational/rebalance.hh>
#include <vacuumms/types.hh>

extern "C"
{
#include <ftw_param.h>
}; 

#define TAU 6.283185307179586

int main(int argc, char** argv)
{
    int n_var_points = 40;
    getIntParam((char*)"-n_var_points", &n_var_points);

    vacuumms_float* var_x = new vacuumms_float[n_var_points];
    vacuumms_float* var_y = new vacuumms_float[n_var_points];
    vacuumms_float* var_z = new vacuumms_float[n_var_points];

    // do a simple helix, 0 < x < TAU; y=cos(x); z=sin(x);
    vacuumms_float start_x = 0.0, start_y = 1.0, start_z = 0.0;
    vacuumms_float end_x = TAU, end_y = 1.0, end_z = 0.0;

    // initialize set of points
    for (int i=0; i<n_var_points; i++)
    {
        // evenly spaced in x
        var_x[i] = start_x + (i + 1) * (end_x - start_x) / (n_var_points + 1);
        var_y[i] = cos(var_x[i]);
        var_z[i] = sin(var_x[i]);
    }

    // dump the points before rebalance
    printf("%f %f %f\n", start_x, start_y, start_z);
    for (int i=0; i<n_var_points; i++) printf("%f %f %f\n", var_x[i], var_y[i], var_z[i]);
    printf("%f %f %f\n", end_x, end_y, end_z);

    // call the update to re-space points:
    std::cout << "rebalancing " << std::endl;
    vacuumms_float shrinkage = rebalance_points_3D(start_x, start_y, start_z, end_x, end_y, end_z, var_x, var_y, var_z, n_var_points);

    // dump the points after rebalance
    printf("%f %f %f\n", start_x, start_y, start_z);
    for (int i=0; i<n_var_points; i++) printf("%f %f %f\n", var_x[i], var_y[i], var_z[i]);
    printf("%f %f %f\n", end_x, end_y, end_z);

    printf("shrinkage: %f\n", shrinkage);
}
