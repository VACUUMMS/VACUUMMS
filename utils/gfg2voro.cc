/* gfg2voro.cc */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Hacking to interface to C linkage; since voro++ is C++, gfg2voro also needs to be C++
extern "C"
{
#include <ftw_param.h>
#include <ftw_pov.h>
}
#include "voro++.hh"

int main(int argc, char *argv[])
{

  // Parse command line options

  setCommandLineParameters(argc, argv);
  if (getFlagParam((char*)"-usage"))
  {
    printf("usage:       gfg2voro    [ -verts || -edges || -povray ]\n");
    printf("                         ( -povray dumps in povray format, otherwise just positions)\n");
    printf("                         -box [ n.nn n.nn n.nn ] \n");
    return 0;
  }
  double box_x=10, box_y=10, box_z=10;
  getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
  int verts = getFlagParam((char*)"-verts"); // dump vertices
  int edges = getFlagParam((char*)"-edges"); // dump edges
  int povray = getFlagParam((char*)"-povray"); // povray format


  // get centers from gfg and dump to temp file for input to voro++

  char temp_name[] = "/tmp/vacuumms-XXXXXX";
  if (mkstemp(temp_name) == -1) 
  {
    printf("vacuumms: could not create temp file %s\n", temp_name);
    exit(1);
  }
  FILE* tempfile=fopen(temp_name, "w");

  int particle_number=0;
  char line[256];
  char *xs, *ys, *zs, *ds, *es;
  double x, y, z, d, e;
  while (++particle_number) // loop over all lines
  {
    fgets(line, 256, stdin);
    if (feof(stdin)) break;
    
    xs = strtok(line, "\t");
    ys = strtok(NULL, "\t");
    zs = strtok(NULL, "\t");
    ds = strtok(NULL, "\t");
    es = strtok(NULL, "\n");

    x = strtod(xs, NULL);
    y = strtod(ys, NULL);
    z = strtod(zs, NULL);
    d = strtod(ds, NULL);
    e = strtod(es, NULL);
    
    // Grab the centers and write to the tempfile
    fprintf(tempfile, "%d %lf %lf %lf\n", particle_number, x, y, z);
  }
  
  // clean up
  fclose(stdin);
  fclose(tempfile);


  // Run voro++ on the temporary structure file
  
  // Set up constants for the container geometry
  const double x_min=0, x_max=box_x;
  const double y_min=0, y_max=box_y;
  const double z_min=0, z_max=box_z;
  
  // Set up the number of blocks that the container is divided into
  const int n_x=6,n_y=6,n_z=6;

  // create container and add particles into the container
  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z, true,true,true,8);
  con.import(temp_name);

  // Create a temp file for output from voro++
  char pov_temp_name[] = "/tmp/vacuumms-pov-XXXXXX";
  if (mkstemp(pov_temp_name) == -1) 
  {
    printf("vacuumms: could not create temp file %s\n", pov_temp_name);
    exit(1);
  }
  con.draw_cells_pov(pov_temp_name);


  // Format and dump the voro++ results back to stdout

  if (verts)
  {
    char command[1024];
    double x0, y0, z0;
    sprintf(command, "grep sphere < %s", pov_temp_name);

    FILE* process = popen(command, "r"); 
    while(fgets(line, 256, process))
    {
      sscanf(line, "sphere{<%lf,%lf,%lf>,r}", &x0, &y0, &z0);
      printf("%lf\t%lf\t%lf\n", x0, y0, z0);
    }
    fclose(process);
  }
  else if (edges)
  {
    char command[1024];
    // cylinder{<3,13,13>,<3,13,3>,r}

    sprintf(command, "grep cylinder < %s", pov_temp_name);
    FILE* process = popen(command, "r"); 
    double x0, y0, z0, x1, y1, z1;

    while(fgets(line, 256, process))
    {
        sscanf(line, "cylinder{<%lf,%lf,%lf>,<%lf,%lf,%lf>,r}", &x0, &y0, &z0, &x1, &y1, &z1);
        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x0, y0, z0, x1, y1, z1);
    }
    fclose(process);
  }
  else if (povray) // just dump the raw file to stdout
  {
    FILE* process = fopen(pov_temp_name, "r");
    while(fgets(line, 256, process)) printf("%s", line); 
    fclose(process);
  }
  

  // clean up

  remove(temp_name);
  remove(pov_temp_name);


  return 0;
}
