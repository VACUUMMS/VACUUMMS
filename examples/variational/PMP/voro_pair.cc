#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vacuumms/variational/variational.hh>
#include <vacuumms/cavity.hh>

extern "C" 
{
#include <ftw_param.h>
}; 

int n_iter = 125;
int n_var_points = 125; 
int mirror_depth = 1;

float sigma = 1.0;
float epsilon = 1.0;

double box_x = 0.0;
double box_y = 0.0;
double box_z = 0.0;

// Gradient descent scaling parameter alpha
float alpha = 0.1;
float alpha_max = 0.5;
float delta_max =0.01;

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getFloatParam((char*)"-alpha", &alpha);
    getFloatParam((char*)"-alpha_max", &alpha_max);
    getIntParam((char*)"-n_iter", &n_iter);
    getIntParam((char*)"-n_var_points", &n_var_points);
    getFloatParam((char*)"-sigma", &sigma);
    getFloatParam((char*)"-epsilon", &epsilon);
    getFloatParam((char*)"-delta_max", &delta_max);
    getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
    getIntParam((char*)"-mirror_depth", &mirror_depth);

    char configuration_filename[] = "PMP.gfg";
    Configuration c = Configuration(configuration_filename);
    c.setBoxDimensions(box_x, box_y, box_z);
    c.setMirrorDepth(mirror_depth);


int index=0;

    // read the var start/end from the var file and create variationals
    while (!feof(stdin))
    {
        vacuumms_float xA, yA, zA;
        vacuumms_float xB, yB, zB;

        // read each pair
        fscanf(stdin, "%f\t%f\t%f\t%f\t%f\t%f\n", &xA, &yA, &zA, &xB, &yB, &zB);

// is this an obtuse pair? From a periodic image?
vacuumms_float d_sq = (xB - xA) * (xB - xA)
                    + (yB - yA) * (yB - yA)             
                    + (zB - zA) * (zB - zA);
vacuumms_float box_sq = box_x * box_x + box_y * box_y + box_z * box_z;
// This excludes pairs with more than 0.25 of box diagonal distance separation. E.g. for 20x20x20 box, this is 8.66.
if (d_sq > 0.0625 * box_sq) 
{
    printf("discarding obtuse pair %f\t%f\t%f - %f\t%f\t%f\n", xA, yA, zA, xB, yB, zB);
    continue;
}

        // create, parameterize, and run the variational for the pair
        printf("creating %d\t%f\t%f\t%f - %f\t%f\t%f\n", index++, xA, yA, zA, xB, yB, zB);
        Variational3D v = Variational3D(xA, yA, zA, xB, yB, zB, sigma, epsilon, n_var_points, &c);

        // set numerical parameters
        v.setAlpha(alpha);
        v.setAlphaMax(alpha_max);
        v.setDeltaMax(delta_max);

        printf("iterating %d times:\n\n", n_iter);
        for (int iter=0; iter < n_iter; iter++) 
        {
            printf("adaptive iteration %d:\n\n", iter);
            v.adaptiveIterateAndUpdate();
        }

        printf("dumping Variational3D object: %p\n", &v);
        v.printValues();
        float* var_x = v.getX();
        float* var_y = v.getY();
        float* var_z = v.getZ();

        /* write var points to be used as  
        for (int point=0; point<n_var_points; point++) 
            printf("###POV_%d-%d\t%f\t%f\t%f\n", i, j, var_x[point], var_y[point], var_z[point]); 
        */

        // draw as cylinders for a continuous trajectory
        for (int point=0; point<n_var_points-1; point++) 
        {
            printf("###CYLINDER\tcylinder{<%f,%f,%f>, <%f,%f,%f>, 0.01 texture{ pigment {color Yellow } }}\n", 
                var_x[point], 
                var_y[point], 
                var_z[point], 
                var_x[point + 1], 
                var_y[point + 1], 
                var_z[point + 1]);
        }

        printf("done. \n\n");
    } // end of stdin
}
