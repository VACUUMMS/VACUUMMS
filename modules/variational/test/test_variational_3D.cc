#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vacuumms/variational/variational.hh>

extern "C" 
{
#include <ftw_param.h>
}; 

int n_iter = 1;
int n_var_points = 25; 
int update = 1;

vacuumms_float start_x = -2.0; 
vacuumms_float start_y = 0.0;
vacuumms_float start_z = -1.0;
vacuumms_float end_x = 2.0;
vacuumms_float end_y = 0.0;
vacuumms_float end_z = 1.0;
vacuumms_float sigma = 1.0;
vacuumms_float epsilon = 1.0;

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getIntParam((char*)"-n_iter", &n_iter);
    getIntParam((char*)"-n_var_points", &n_var_points);

    char filename[] = "x.gfg";
    Configuration c = Configuration(filename);

    printf("dumping configuration:\n");
    c.dumpContents();
    printf("done. \n\n");

    Variational3D v = Variational3D(start_x, start_y, start_z, end_x, end_y, end_z, sigma, epsilon, n_var_points, &c);

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

