/* gfg2voro.cc */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Hacking to interface to C linkage; since voro++ is C++, this also needs to be C++
extern "C"
{
#include <ftw_param.h>
#include <ftw_pov.h>
}
#include "voro++.hh"

/*
#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=-5,x_max=5;
const double y_min=-5,y_max=5;
const double z_min=0,z_max=10;

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

int main() {

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
//      container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
//                      false,false,false,8);
        container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                        true,true,true,8);

        //Randomly add particles into the container
        con.import("pack_ten_cube");

        // Save the Voronoi network of all the particles to text files
        // in gnuplot and POV-Ray formats
        con.draw_cells_gnuplot("pack_ten_cube.gnu");
        con.draw_cells_pov("pack_ten_cube_v.pov");

        // Output the particles in POV-Ray format
        con.draw_particles_pov("pack_ten_cube_p.pov");
}

// end of example */

int main(int argc, char *argv[])
{
  char line[256];
  char *xs, *ys, *zs, *ds, *es;
  double x, y, z, d, e;
  const char *color="White";
  double box_x=10, box_y=10, box_z=10;
//  int verts = 0;
//  int edges = 0;
  char transmit_str[256] = "";
  double phong = 0.0;
  char phong_str[256] = "";

  setCommandLineParameters(argc, argv);
  if (getFlagParam((char*)"-usage"))
  {
    printf("usage:       gfg2voro    [ -verts || -edges || -povray ]\n");
    printf("                         ( -povray dumps in povray format, otherwise just positions)\n");
    printf("                         -box [ n.nn n.nn n.nn ] \n");
    return 0;
  }

  // Create a temp file for input to voro++
  char temp_name[] = "/tmp/vacuumms-XXXXXX";
  int fd = mkstemp(temp_name);
  printf("FTW: got temp_name = >>>%s<<<\n", temp_name);
  FILE* tempfile=fopen(temp_name, "w");

  // Parse command line
  getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
  int verts = getFlagParam((char*)"-verts"); // dump vertices
  int edges = getFlagParam((char*)"-edges"); // dump edges
  int povray = getFlagParam((char*)"-povray"); // povray format

  printf("FTW: got box dims = %lf x %lf x %lf \n", box_x, box_y, box_z);

  int particle_number=0;
  //while (1) // loop over all lines
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
    // printf("%lf,%lf,%lf\n", x, y, z);
    fprintf(tempfile, "%d %lf %lf %lf\n", particle_number, x, y, z);
  }
  fclose(stdin);
  fclose(tempfile);


  // Now do the voro++ stuff.
  
  // Set up constants for the container geometry
  const double x_min=0, x_max=box_x;
  const double y_min=0, y_max=box_y;
  const double z_min=0, z_max=box_z;
  
  // Set up the number of blocks that the container is divided into
  const int n_x=6,n_y=6,n_z=6;

  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z, true,true,true,8);
  //Randomly add particles into the container
  con.import(temp_name);

  // Create a temp file for output from voro++
  char pov_temp_name[] = "/tmp/vacuumms-pov-XXXXXX";
  int pov_fd = mkstemp(pov_temp_name);
  printf("Got pov_temp_name = %s\n", pov_temp_name);
  printf("Got pov_fd = %d\n", pov_fd);
  con.draw_cells_pov(pov_temp_name);


  // Dump the results back to stdout

  if (verts)
  {
      char command[1024];
      sprintf(command, "grep sphere < %s | sed -e \"s/sphere{<//g\" | sed -e \"s/>.*//g\" |tr \",\" \"\t\"", pov_temp_name);

      FILE* process = popen(command, "r"); 
      while(1)
      {
        fgets(line, 256, process);
        if (feof(process)) break;
        printf("%s", line);
      }
      fclose(process);
  }
  else if (edges)
  {
      char command[1024];
      // cylinder{<3,13,13>,<3,13,3>,r}

      sprintf(command, "grep cylinder < %s", pov_temp_name);
      FILE* process = popen(command, "r"); 
      double x0,y0,z0,x1,y1,z1;
// need to find a way to break loop
      while(1)
      {
//        fgets(line, 256, process);
        int retval = fscanf(process, "cylinder{<%lf,%lf,%lf>,<%lf,%lf,%lf>,r}", &x0, &y0, &z0, &x1, &y1, &z1);
printf("FTW retval is %d\n", retval);
if (retval == 0) break;
        if (feof(process)) break;
        printf("got <%lf,%lf,%lf> - <%lf,%lf,%lf>!\n", x0, y0, z0, x1, y1, z1);
      }
      fclose(process);
  } 
 
  
//  FILE *popen(const char *command, const char *mode); [Option End]

  // clean up
//  remove(temp_name);
//  remove(pov_temp_name);

  return 0;
}

