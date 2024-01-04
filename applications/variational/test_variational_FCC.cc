#include <iostream>
#include <stdio.h>
#include <math.h>
#include "variational.hh"

extern "C" 
{
#include <ftw_param.h>
}; 

int n_iter = 1;
int n_var_points = 5; 
int update = 1;
float *var_x, *var_y;

float start_x = 0; 
float start_y = 0;
float start_z = 0.5;
float end_x = 0.5;
float end_y = 0.5;
float end_z = 0.5;

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getIntParam((char*)"-n_iter", &n_iter);
    getIntParam((char*)"-n_var_points", &n_var_points);

    char filename[] = "fcc.gfg";
    Configuration c = Configuration(filename);

    printf("dumping configuration:\n");
    c.dumpContents();
    printf("done. \n\n");

    Variational3D v = Variational3D(start_x, start_y, start_z, end_x, end_y, end_z, n_var_points, &c);
    printf("dumping new Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

    printf("rebalancing:\n\n");
    v.rebalancePoints3D();
    printf("done. \n\n");


/*
    printf("iterating %d times:\n\n", n_iter);
    for (int i=0; i < n_iter; i++) 
    {
        printf("iterating:\n\n");
        v.iterate();
        v.printValues();
        printf("rebalancing:\n\n");
        v.rebalancePoints3D();
        v.printValues();
    }

    printf("dumping Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

//    printf("rebalancing:\n\n");
//    v.rebalancePoints3D();

//    printf("pseudo-rebalancing:\n\n");
//    v.pseudoRebalance();

    printf("after rebalancing:\n\n");
    printf("dumping Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");
*/
}

