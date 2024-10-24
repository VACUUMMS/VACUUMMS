/* io_setup.c */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <vacuumms/std.h>
#include <vacuumms/rng.h>

/* prototypes from old io_setup.h */
void setInitialConditions();
void generateUniqueId();
void finalizeOutput();
void initializeOutput();
void generateOutput();
void loadConfiguration();
void readEnvironmentVariables();

extern char *simulation_unique_identifier;
extern double temperature;
extern int number_of_molecules;
extern int energy_report_frequency;
extern int configuration_threshold;
extern int configuration_frequency;
extern double box_x, box_y, box_z;
extern int monte_carlo_steps;
extern int end_mcs;
extern char* log_file_name;
extern char* output_file_name;
extern char* input_file_name;
extern double x[], y[], z[];
extern int x_laps[], y_laps[], z_laps[];

FILE *log_file;
FILE *output_file;
FILE *input_file;

extern char hostname[50];
char *log_path;
char *results_path;

time_t now;

extern int verbose;

void setInitialConditions()
{
  int i, j, k;
  int num_so_far=0;

  if (input_file_name != NULL) loadConfiguration();
  else 
  {
    monte_carlo_steps = 0;
    for (i=0; i<box_x; i++)
    for (j=0; j<box_y; j++)
    for (k=0; k<box_z; k++)
    {
      x[num_so_far] = i;   
      y[num_so_far] = j;   
      z[num_so_far] = k;   

      if (num_so_far++ >= number_of_molecules) return;
    }

    printf ("too many molecules...\n");
    exit(0);
  }
}

void loadConfiguration()
{
  FILE *datastream;
  char line[80];
  char *xs, *ys, *zs;

  number_of_molecules = 0;
  V printf("loading %s...\n", input_file_name);
  datastream = fopen(input_file_name, "r");

  while (1)
  {
    fgets(line, 80, datastream);
    if (feof(datastream)) break;

    xs = strtok(line, "\t");
    ys = strtok(NULL, "\t");
    zs = strtok(NULL, "\n");

    x[number_of_molecules] = strtod(xs, NULL);
    y[number_of_molecules] = strtod(ys, NULL);
    z[number_of_molecules++] = strtod(zs, NULL);
  }
 
  V printf("%d lines read.\n", number_of_molecules);
  fclose(datastream);
}

void generateUniqueId()
{
  int i;
  if (!strcmp(simulation_unique_identifier, "################"))
    for (i=0; i<16; i++) *(simulation_unique_identifier + i) = (char)(rnd() * 26 + 65);
}

void readEnvironmentVariables()
{
  gethostname(hostname, 50);
  log_path = getenv("LOG_PATH");
  log_file_name = strcat(log_path, "/");
  log_file_name = strcat(log_file_name, hostname);
  log_file_name = strcat(log_file_name, "-hs.log");
  results_path = getenv("RESULTS_PATH");
  output_file_name = strcat(results_path, "/");
  output_file_name = strcat(output_file_name, simulation_unique_identifier);
  output_file_name = strcat(output_file_name, "-hs.out");
printf("output_file_name = %s\n", output_file_name);
}

/* display headers on output */
void initializeOutput()
{
  now = time(NULL);

  if (verbose)
  {
    printf("#HT simulation %s started on %s:  %s", simulation_unique_identifier, hostname, ctime(&now));
    printf("#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n");
    printf("#H N=%d\n", number_of_molecules);
    printf("#H box dimensions:  \t%lf x %lf x %lf = %lf\n", box_x, box_y, box_z, box_x*box_y*box_z);
    printf("#H reduced density = %lf\n", (number_of_molecules / (box_x * box_y * box_z)));
    printf("#H\n");
    printf("#H configuration threshold:\t%d\n", configuration_threshold);
    printf("#H configuration frequency:\t%d\n", configuration_frequency);
    printf("#H run until mcs = %d\n", end_mcs);
    printf("#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n");
    printf("#H\n");
    printf("#H\n");
  }

  log_file = fopen(log_file_name, "a");
  fprintf(log_file, "simulation %s launched on %s:  %s", simulation_unique_identifier, hostname, ctime(&now));
  fclose(log_file);

  output_file = fopen(output_file_name, "w");
  fprintf(output_file, "#HT simulation %s started on %s:  %s", simulation_unique_identifier, hostname, ctime(&now));
  fprintf(output_file, "#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n");
  fprintf(output_file, "#H N=%d\n", number_of_molecules);
  fprintf(output_file, "#H box dimensions:  \t%lf x %lf x %lf = %lf\n", box_x, box_y, box_z, box_x*box_y*box_z);
  fprintf(output_file, "#H reduced density = %lf\n", (number_of_molecules / (box_x * box_y * box_z)));
  fprintf(output_file, "#H\n");
  fprintf(output_file, "#H configuration threshold:\t%d\n", configuration_threshold);
  fprintf(output_file, "#H configuration frequency:\t%d\n", configuration_frequency);
  fprintf(output_file, "#H run until mcs = %d\n", end_mcs);
  fprintf(output_file, "#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n");
  fprintf(output_file, "#H\n");
  fprintf(output_file, "#H\n");
  fflush(output_file);
}

void finalizeOutput()
{
  now = time(NULL);
  if (verbose) printf("#HT simulation %s finished on %s:  %s", simulation_unique_identifier, hostname, ctime(&now));
  fprintf(output_file, "#HT simulation %s finished on %s:  %s", simulation_unique_identifier, hostname, ctime(&now));
  fclose(output_file);

  log_file = fopen(log_file_name, "a");
  fprintf(log_file, "simulation %s finished on %s:  %s", simulation_unique_identifier, hostname, ctime(&now));
  fclose(log_file);
}

void generateOutput()
{
  int i;

  if ((monte_carlo_steps > configuration_threshold) && (monte_carlo_steps % configuration_frequency == 0))
  {
    now = time(NULL);
    log_file = fopen(log_file_name, "a");
    fprintf(log_file, "dumping configuration for %s on %s at %d steps:  %s", hostname, \
            simulation_unique_identifier, monte_carlo_steps, ctime(&now));
    fclose(log_file);
/*
for (i=0; i<number_of_molecules; i++)
{
printf("%d\t%d\t%d\t%d\n", i, x_laps[i], y_laps[i], z_laps[i]);
printf("%d\t%lf\t%lf\t%lf\n", i, x[i], y[i], z[i]);
}
*/

    if (verbose)
    {
      printf("#HC%06d\n", monte_carlo_steps);
      printf("#HC%06d dumping configuration at mcs=%d...\n", monte_carlo_steps, monte_carlo_steps);
      printf("#HC%06d\n", monte_carlo_steps);
      for (i=0; i<number_of_molecules; i++) printf("#C%06d\t%d\t%lf\t%lf\t%lf\n", monte_carlo_steps, i, 
                                                    box_x * x_laps[i] + x[i], 
                                                    box_y * y_laps[i] + y[i],
                                                    box_z * z_laps[i] + z[i]);
      printf("#HC%06d\n", monte_carlo_steps);
    }

    fprintf(output_file, "#HC%06d\n", monte_carlo_steps);
    fprintf(output_file, "#HC%06d dumping configuration at mcs=%d...\n", monte_carlo_steps, monte_carlo_steps);
    fprintf(output_file, "#HC%06d\n", monte_carlo_steps);
    for (i=0; i<number_of_molecules; i++) fprintf(output_file, "#C%06d\t%d\t%lf\t%lf\t%lf\n", monte_carlo_steps, i,
                                                    box_x * x_laps[i] + x[i], 
                                                    box_y * y_laps[i] + y[i],
                                                    box_z * z_laps[i] + z[i]);
    fprintf(output_file, "#HC%06d\n", monte_carlo_steps);
  }

  fflush(output_file);
}

