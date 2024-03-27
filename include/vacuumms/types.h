/* vacuumms/types.h */

/* Define as single or double precision.
 * This is new for 1.2.x, first in the variational module, and to be integrated
 * more widely later.
 */

#ifndef VACUUMMS_TYPES
#define VACUUMMS_TYPES

#define vacuumms_float float

#include <vacuumms/limits.h>


typedef struct vacuumms_vector
{
  double x;
  double y;
  double z;
} vacuumms_vector;

struct Atom
{
  float x;
  float y;
  float z;
  float sigma;
  float epsilon;
};

typedef struct Atom vacuumms_Atom;

struct vacuumms_Cavity
{
  float x;
  float y;
  float z;
  float diameter;
};

typedef struct vacuumms_Cavity vacuumms_Cavity;

/* Structure to hold a set of cavities */
typedef struct CAV65536
{
  vacuumms_Cavity cavity[65536];
  float box_x;
  float box_y;
  float box_z;
  int n_cavities;
} vacuumms_CAV65536;

struct vacuumms_Configuration
{
  vacuumms_Atom *atom;
  double box_x;
  double box_y;
  double box_z;
  int n_atoms;
};

typedef struct vacuumms_Configuration vacuumms_Configuration;

struct CommandLineOptions
{ 
  int verbose;
  float T;
  float verlet_cutoff_sq;
  float box_x;
  float box_y;
  float box_z;
  int n_threads;
  float r_i;
  float epsilon_i;
  int molecule;
  char *config_directory;
  char *config_name;
  int use_stdin;
  float gross_resolution;
  float fine_resolution;
};

typedef vacuumms_Configuration vacuumms_GFG;

/*  Adding the following definition because CUDA doesn't handle deep copies */
// struct GFG65536
// {
//   vacuumms_Atom atom[65536];
//   float box_x;
//   float box_y;
//   float box_z;
//   int n_atoms;
// };

/* Could it really be as easy as just making this number bigger? */
struct GFG65536
{
  vacuumms_Atom atom[VACUUMMS_MAX_NUMBER_OF_MOLECULES * 27];
  float box_x;
  float box_y;
  float box_z;
  int n_atoms;
};

typedef struct GFG65536 vacuumms_GFG65536;

struct EnergyArray256
{
  float energy[256][256][256];
};

typedef struct EnergyArray256 vacuumms_EnergyArray256;

struct EnergyArray512
{
  float energy[512][512][512];
};

typedef struct EnergyArray512 vacuumms_EnergyArray512;

struct EnergyArray1024
{
  float energy[1024][1024][1024];
};

typedef struct EnergyArray1024 vacuumms_EnergyArray1024;

struct FVI256
{
  float intensity[256][256][256];
};

typedef struct FVI256 vacuumms_FVI256;

struct EnergyArray
{
  float ***energy;
};

typedef struct EnergyArray vacuumms_EnergyArray;

struct FVI
{
  float ***intensity;
};

typedef struct FVI vacuumms_FVI;

#endif

