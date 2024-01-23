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

float sigma = 0.5;
float epsilon = 1.0;

// delete atoms and calc path
float start_x = 0.0; 
float start_y = 1.0;
float start_z = 1.0;
float end_x = 1.5;
float end_y = 1.5;
float end_z = 1.0;

/* original edge, not in voronoi
float start_x = 0; 
float start_y = 0;
float start_z = 0.5;
float end_x = 0.5;
float end_y = 0.5;
float end_z = 0.5;
*/

/*
// this edge from voronoi
float start_x = 0.25; 
float start_y = 0.25;
float start_z = 0.25;
float end_x = 0.0;
float end_y = 0.5;
float end_z = 0.0;
*/

float alpha = 0.1;
float alpha_max = 1.0;

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getFloatParam((char*)"-alpha", &alpha);
    getIntParam((char*)"-n_iter", &n_iter);
    getIntParam((char*)"-n_var_points", &n_var_points);

    char filename[] = "fcc.gfg";
    Configuration c = Configuration(filename);
    c.setBoxDimensions(2,2,2);

    printf("dumping configuration:\n");
    c.dumpContents();
    printf("done. \n\n");

    Variational3D v = Variational3D(start_x, start_y, start_z, end_x, end_y, end_z, sigma, epsilon, n_var_points, &c);
    v.setAlpha(alpha);
    v.setAlphaMax(alpha_max);
    printf("dumping new Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

    printf("rebalancing:\n\n");
    float shrinkage=v.rebalancePoints3D();
    printf("Shrinkage: %f\n", shrinkage);
    printf("done. \n\n");

    printf("dumping new Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");


    printf("iterating %d times:\n\n", n_iter);
    for (int i=0; i < n_iter; i++) 
    {
        printf("###############################################################################################################\n");
//        printf("iteration %d:\n\n", i);
//        v.iterate();
        printf("adaptive iteration %d:\n\n", i);
        v.adaptiveIterateAndUpdate();
        v.printValues();
//        printf("###############################################################################################################\n");
//        printf("rebalancing:\n\n");
//        float shrinkage=v.rebalancePoints3D();
//        printf("Shrinkage: %f\n", shrinkage);
//        v.printValues();
        printf("###############################################################################################################\n");
    }

    printf("dumping Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

//    printf("rebalancing:\n\n");
//    v.rebalancePoints3D();

//    printf("pseudo-rebalancing:\n\n");
//    v.pseudoRebalance();

//    printf("after rebalancing:\n\n");
//    printf("dumping Variational3D object: %p\n", &v);

//    v.printValues();
//    printf("done. \n\n");

}

