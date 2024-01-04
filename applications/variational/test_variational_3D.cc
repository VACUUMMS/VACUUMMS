#include <iostream>
#include <stdio.h>
#include <math.h>
#include "variational.hh"

extern "C" 
{
#include <ftw_param.h>
}; 

//int n_points = 2;
int n_iter = 1;
int n_var_points = 25; 
int update = 1;
//float var_x[n_var_points], var_y[n_var_points];
float *var_x, *var_y;

float start_x = -2; 
float start_y = 0;
float start_z = -1;
float end_x = 2;
float end_y = 0;
float end_z = 1;

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getIntParam((char*)"-n_iter", &n_iter);
//    getIntParam((char*)"-update", &update);
    getIntParam((char*)"-n_var_points", &n_var_points);

    char filename[] = "x.gfg";
    Configuration c = Configuration(filename);

    printf("dumping configuration:\n");
    c.dumpContents();
    printf("done. \n\n");

    Variational3D v = Variational3D(start_x, start_y, start_z, end_x, end_y, end_z, n_var_points, &c);
    printf("dumping Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

    printf("iterating %d times:\n\n", n_iter);
    for (int i=0; i < n_iter; i++) 
    {
        v.iterate();
    }

    printf("dumping Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

    printf("skipping rebalancing:\n\n");
//    v.rebalancePoints3D();
    printf("dumping Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

}

