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

//float sigma = 0.707106781186548;
float sigma = 1.0;
float epsilon = 1.0;

// delete atoms and calc path
// 2.554920	1.927518	1.946361	1.0	1.0
// 0.622468	0.276067	0.880539	1.0	1.0

// 0.212469	0.580506	1.796850	1.0	1.0
// 2.272498	0.502913	1.765567	1.0	1.0

double start_x = 0.212469;
double start_y = 0.580506;
double start_z = 1.796850;

double end_x = 2.272498;
double end_y = 0.502913;
double end_z = 1.765567;

double box_x = 3.174802103936399;
double box_y = 3.174802103936399;
double box_z = 3.174802103936399;

float alpha = 0.1;
float alpha_max = 1.0;

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getFloatParam((char*)"-alpha", &alpha);
    getFloatParam((char*)"-alpha_max", &alpha_max);
    getIntParam((char*)"-n_iter", &n_iter);
    getIntParam((char*)"-n_var_points", &n_var_points);
    getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
    getVectorParam((char*)"-start", &start_x, &start_y, &start_z);
    getVectorParam((char*)"-end", &end_x, &end_y, &end_z);

    char filename[] = "lj.gfg";
    Configuration c = Configuration(filename);
    c.setBoxDimensions(box_x, box_y, box_z);

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
        printf("adaptive iteration %d:\n\n", i);
        v.adaptiveIterateAndUpdate();
        v.printValues();
        printf("###############################################################################################################\n");
    }

    printf("dumping Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");
}

