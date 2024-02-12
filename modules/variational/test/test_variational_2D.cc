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
float end_x = 2;
float end_y = 0;

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

    Variational2D v = Variational2D(start_x, start_y, end_x, end_y, n_var_points, &c);
    printf("dumping Variational2D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

    printf("iterating %d times:\n\n", n_iter);
    for (int i=0; i < n_iter; i++) 
    {
        v.iterate();
    }

    printf("dumping Variational2D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

    printf("rebalancing:\n\n");
    v.rebalancePoints2D();
    printf("dumping Variational2D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

}

