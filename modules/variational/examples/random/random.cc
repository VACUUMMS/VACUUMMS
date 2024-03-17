#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vacuumms/variational/variational.hh>

#include <vacuumms/param.hh>

int n_iter = 1;
int n_var_points = 5; 
int update = 1;
float *var_x, *var_y;

// a single edge pulled from the voronoi graph
float start_x = 1 + 0.969667;
float start_y = 1 + 0.560653;
float start_z = 1 + 0.0434551;
float end_x = 1 + 0.290834;
float end_y = 1 + -0.21678;
float end_z = 1 + 0.529496;

float sigma = 1.0;
float epsilon = 1.0;

float alpha = 0.05;

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getFloatParam((char*)"-alpha", &alpha);
    getIntParam((char*)"-n_iter", &n_iter);
    getIntParam((char*)"-n_var_points", &n_var_points);

    char filename[] = "rando.gfg";
    Configuration c = Configuration(filename);

    printf("dumping configuration:\n");
    c.dumpContents();
    printf("done. \n\n");

    Variational3D v = Variational3D(start_x, start_y, start_z, end_x, end_y, end_z, sigma, epsilon, n_var_points, &c);
    v.setAlpha(alpha);
    printf("dumping original Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

    printf("iterating %d times:\n\n", n_iter);
    for (int i=0; i < n_iter; i++) 
    {
        printf("######################################################################################################\n");
        printf("iteration %d:\n\n", i);
        v.adaptiveIterateAndUpdate();
        v.printValues();
        printf("######################################################################################################\n");
    }

    printf("dumping final Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");
}

