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

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getFloatParam((char*)"-alpha", &alpha);
    getFloatParam((char*)"-alpha_max", &alpha_max);
    getIntParam((char*)"-n_iter", &n_iter);
    getIntParam((char*)"-n_var_points", &n_var_points);
    getFloatParam((char*)"-sigma", &sigma);
    getFloatParam((char*)"-epsilon", &epsilon);
    getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
    getIntParam((char*)"-mirror_depth", &mirror_depth);

    char configuration_filename[] = "PMP.gfg";
    Configuration c = Configuration(configuration_filename);
    c.setBoxDimensions(box_x, box_y, box_z);
    c.setMirrorDepth(mirror_depth);

    printf("dumping original configuration:\n");
//    c.dumpContents();
    printf("done. \n\n");

    // loop over all pairs of cavities 
    char cavity_filename[] = "PMP.unq";
    CavityConfiguration cavs = CavityConfiguration(cavity_filename);

vacuumms_float threshold = 2.5;
printf("using threshold %f\n", threshold);
    
    for (int i=0; i<cavs.getSize(); i++) 
    {
        if (cavs.recordAt(i).d < threshold)
        {
            printf("deleting %d: %f\t%f\t%f\t%f\n", i, cavs.recordAt(i).x, cavs.recordAt(i).y, cavs.recordAt(i).z, cavs.recordAt(i).d); 
            cavs.deleteRecordAt(i--);
        }
        //printf("%f\t%f\t%f\t%f\n", cavs.recordAt(i).x, cavs.recordAt(i).y, cavs.recordAt(i).z, cavs.recordAt(i).d); 
    }
printf("%d cavs.\n", cavs.getSize());
    
    int configuration_size = cavs.getSize();

    for (int i=0; i<configuration_size; i++)
        for (int j=i+1; j<configuration_size; j++)
        {
            // select the cavity pair
            Cavity start = cavs.recordAt(i);
            Cavity end = cavs.recordAt(j);

            // create, parameterize, and run the variational for the pair
            Variational3D v = Variational3D(start.x, start.y, start.z, end.x, end.y, end.z, sigma, epsilon, n_var_points, &c);
            v.setAlpha(alpha);
            v.setAlphaMax(alpha_max);

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

            for (int point=0; point<n_var_points; point++) printf("###POV_%d-%d\t%f\t%f\t%f\n", i, j, var_x[point], var_y[point], var_z[point]); 

            // draw as cylinders for a continuous trajectory
            for (int point=0; point<n_var_points-1; point++) 
                printf("###CYLINDER\tcylinder{<%f,%f,%f>, <%f,%f,%f>, 0.01 texture{ pigment {color Yellow } }}\n", var_x[point], var_y[point], var_z[point], var_x[point + 1], var_y[point + 1], var_z[point + 1]);

//cylinder{<1.0493,0.0344035,0.224887>,<1.23868,0.77789,-0.232467>, 0.0025 texture{ pigment {color Yellow } }}

            printf("done. \n\n");
        } // end loop over (i,j) pairs
}
