/* vddx_io_setup.c */

#include <ftw_std.h>
#include <ftw_rng.h>
#include <time.h>
#include <string.h>

#include <vacuumms/types.h>

extern vacuumms_float vertex_x[];
extern vacuumms_float vertex_y[];
extern vacuumms_float vertex_z[];

int loadVertices(char* vertices_filename)
{
  char line[256];
  char *xs, *ys, *zs;

  int number_of_vertices = 0;
  FILE* vertices_file = fopen(vertices_filename, "r");
  if (vertices_file == NULL) return -1;

  while (1)
  {
    fgets(line, 256, vertices_file);
    if (feof(vertices_file)) break;

    xs = strtok(line, "\t");
    ys = strtok(NULL, "\t");
    zs = strtok(NULL, "\n");

    vertex_x[number_of_vertices] = strtod(xs, NULL);
    vertex_y[number_of_vertices] = strtod(ys, NULL);
    vertex_z[number_of_vertices] = strtod(zs, NULL);
    number_of_vertices++;
  }
 
  fclose(vertices_file);

  return number_of_vertices;
}

