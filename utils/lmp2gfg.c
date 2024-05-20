// parser.c

#include <vacuumms/config_parser.h>
#include <vacuumms/param.h>
#include <vacuumms/types.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

int main(int argc, char *argv[]) 
{
  double box_x=10, box_y=10, box_z=10;
  char* pair_filename = "pair.dat";
  char* lmps_filename = "in.lmps";

  setCommandLineParameters(argc, argv);

  if (getFlagParam("-usage"))
  {
    printf("\n");
    printf("usage:     lmp2gfg        -box [10.0 10.0 10.0] \n");
    printf("                          -pair_file [pair.dat] \n");
    printf("                          -lmps_file [in.lmps]  \n");
    printf("\n");
    exit(0);
  }

  getVectorParam("-box", &box_x, &box_y, &box_z);
  getStringParam("-pair_file", &pair_filename);
  getStringParam("-lmps_file", &lmps_filename);

printf("pair file: %s\n", pair_filename);
printf("lmps file: %s\n", lmps_filename);

  FILE* pair_file = fopen(pair_filename, "r");
assert(pair_file);

  FILE* lmps_file = fopen(lmps_filename, "r");
assert(lmps_file);

char line[256];
int index;
vacuumms_float epsilon, sigma;

//  while(fscanf(pair_file, "%d %f %f", &index, &epsilon, &sigma) != EOF)
while(fgets(line, 256, lmps_file))
{
    float xlo, xhi;
    char* retval = strstr(line, "xlo");
    if (retval !=NULL)
    {
         sscanf(line, "%f %f xlo xhi", &xlo, &xhi);
printf("Got xlo, xhi = %f, %f\n", xlo, xhi);
    }

     // printf("Scanned: %d\t%f\t%f\n", index, epsilon, sigma);
    printf("skipping...>>>%s", line);
    //char* retval = strstr(line, "Atoms");
    retval = strstr(line, "Atoms");
    if (retval != NULL) 
    {
        // found it, so do the work
        int dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
        vacuumms_float f_dummy, x, y, z;
        while(fscanf(lmps_file, "%d %d %d %f %f %f %f %d %d %d", &dummy1, &dummy2, &dummy3, &f_dummy, &x, &y, &z, &dummy4, &dummy5, &dummy6) > 0)
            printf("got atom at (%f, %f, %f)\n", x, y, z);
       // this will return zero when scanf fails, dropping back to outer loop, with seek set to end of atom list.. i think... 
    }


    // now find pair stuff
    retval = strstr(line, "Pair");
    if (retval != NULL) 
    {
        // found it, so do the work
        int atom_type;
        vacuumms_float sigma, epsilon;
        while(fscanf(lmps_file, "%d %f %f", &atom_type, &epsilon, &sigma) > 0)
            printf("got pair (%d, %f, %f)\n", atom_type, epsilon, sigma);
       // this will return zero when scanf fails, dropping back to outer loop, with seek set to end of atom list.. i think... 
    }
}

/*
FILE* process = popen(command, "r");
while(fgets(line, 256, process))
{
double x0, y0, z0;
sscanf(line, "sphere{<%lf,%lf,%lf>,r}", &x0, &y0, &z0);
printf("%lf\t%lf\t%lf\n", x0, y0, z0);
}
fclose(process);
*/




/*
  fprintf(stderr,"reading configuration\n");
  fprintf(stderr, "calculating resolution = %d for %d potential\n", resolution, potential);
  fprintf(stderr, "using sigma = %f and epsilon = %f \n", sigma, epsilon);
*/

  return 0;

}
