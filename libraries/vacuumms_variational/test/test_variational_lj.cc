#include <iostream>
//#include <stdio.h>
//#include <math.h>

#include <vacuumms/variational/variational.hh>
#include <vacuumms/parameters.hh>

int n_iter = 1;
int n_var_points = 5; 

//float sigma = 0.707106781186548;
vacuumms_float sigma = 1.0;
vacuumms_float epsilon = 1.0;

// delete atoms and calc path
// 2.554920	1.927518	1.946361	1.0	1.0
// 0.622468	0.276067	0.880539	1.0	1.0

// 0.212469	0.580506	1.796850	1.0	1.0
// 2.272498	0.502913	1.765567	1.0	1.0

vacuumms_float start_x = 0.212469;
vacuumms_float start_y = 0.580506;
vacuumms_float start_z = 1.796850;

vacuumms_float end_x = 2.272498;
vacuumms_float end_y = 0.502913;
vacuumms_float end_z = 1.765567;

vacuumms_float box_x = 3.174802103936399;
vacuumms_float box_y = 3.174802103936399;
vacuumms_float box_z = 3.174802103936399;

vacuumms_float alpha = 0.1;
vacuumms_float alpha_max = 1.0;

int main(int argc, char** argv)
{
    Parameters p(argc, argv);
    p.getFloatParam((char*)"-alpha", &alpha);
    p.getFloatParam((char*)"-alpha_max", &alpha_max);
    p.getIntParam((char*)"-n_iter", &n_iter);
    p.getIntParam((char*)"-n_var_points", &n_var_points);
    p.getVectorParam((char*)"-start", &start_x, &start_y, &start_z);
    p.getVectorParam((char*)"-end", &end_x, &end_y, &end_z);

    // using the new containers
    // std::vector<vacuumms_float> box_dims = p.getVectorParam(std::string("-box"));
    std::vector<vacuumms_float> box_dims = {box_x, box_y, box_z};

    //char filename[] = "ljx.gfg";
    const char *filename = p.getStringParam((char*)"-filename");
    Configuration c = Configuration(filename);
    c.setBoxDimensions(box_dims);

    printf("dumping configuration:\n");
    c.dumpContents();
    printf("done. \n\n");

    Variational3D v = Variational3D(start_x, start_y, start_z, end_x, end_y, end_z, sigma, epsilon, n_var_points, &c);
    v.setAlpha(alpha);
    v.setAlphaMax(alpha_max);
    // printf("dumping new Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");

    printf("rebalancing:\n\n");
    float shrinkage=v.rebalancePoints3D();
    printf("Shrinkage: %f\n", shrinkage);
    printf("done. \n\n");

    // printf("dumping new Variational3D object: %p\n", &v);
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

    // printf("dumping Variational3D object: %p\n", &v);
    v.printValues();
    printf("done. \n\n");
}

