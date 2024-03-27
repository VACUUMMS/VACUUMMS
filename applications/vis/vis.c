/*************************************** vis.c ********************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <vacuumms/io_setup.h>
#include <vacuumms/graphics.h>
#include <vacuumms/science.h>
#include <vacuumms/std.h>

#ifndef MAX_NUMBER_MOLECULES
#define MAX_NUMBER_MOLECULES 16384
#endif

/* non-configurable global params */
double x[MAX_NUMBER_MOLECULES], y[MAX_NUMBER_MOLECULES], z[MAX_NUMBER_MOLECULES], d[MAX_NUMBER_MOLECULES];
int c[MAX_NUMBER_MOLECULES];

int wsize_x, wsize_y, wsize_z;

int number_of_molecules;

void parseCommandLineOptions(int argc, char *argv[]);

int main(int argc, char *argv[])
{
  parseCommandLineOptions(argc, argv);
printf("%d\n",   loadConfiguration());
  initializeDisplay();

  while (1)
  {
    drawGraphicalRepresentation();
    checkForWindowEvent();
    sleep(5);
  }

  return 0;

} /* end main */
