// vddx.c   This is shmem/pthread version of ddx starting from voronoi vertices 

#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <semaphore.h>
#include <unistd.h>
#include <math.h>

#include <vacuumms/types.h>
#include <ftw_param.h>
#include <ftw_prng.h>
#include <ftw_types.h>
#include <ftw_std.h>
#include <io_setup.h>

int loadVertices(char* vertices_file);

// setting defaults for ~640MB bss for 16M MOLECULES

#ifndef MAX_NUM_MOLECULES
#define MAX_NUM_MOLECULES 16777216
#endif

// setting defaults for 60MB of heap/thread for Verlet list of 1M

#ifndef MAX_CLOSE
#define MAX_CLOSE 1048576
#endif

typedef struct 
{
  int				thread_id;
  int   			close_molecules;
  int   			attempts;
  double 			test_x0, test_y0, test_z0;
  double 			test_x, test_y, test_z;
  double 			verlet_center_x, verlet_center_y, verlet_center_z;
  double 			diameter;
  double 			close_x[MAX_CLOSE], close_y[MAX_CLOSE], close_z[MAX_CLOSE];
  double 			close_sigma[MAX_CLOSE];
  double 			close_sigma6[MAX_CLOSE];
  double 			close_sigma12[MAX_CLOSE];
  double 			close_epsilon[MAX_CLOSE];
  double 			sq_distance_from_initial_pt;
  struct MersenneTwister 	rng;
} Trajectory;

double calculateRepulsion(Trajectory*);
double calculateEnergy(Trajectory*, double);
void generateTestPoint(Trajectory*);
double findEnergyMinimum(Trajectory*);
void makeVerletList(Trajectory*);
void expandTestParticle(Trajectory*);
void* ThreadMain(void *threadID);

pthread_t* 		threads;	// the threads
void**			passvals;	// values passed to each thread
sem_t			semaphore;	// semaphore to restrict number of threads running
sem_t			completion_semaphore;	// semaphore to count # of completed threads
pthread_mutex_t 	mutex;
int                     thread_idx; 

vacuumms_float vertex_x[MAX_NUM_MOLECULES];
vacuumms_float vertex_y[MAX_NUM_MOLECULES];
vacuumms_float vertex_z[MAX_NUM_MOLECULES];

// Accessible to all threads
double 	x[MAX_NUM_MOLECULES];
double 	y[MAX_NUM_MOLECULES];
double 	z[MAX_NUM_MOLECULES];
double 	sigma[MAX_NUM_MOLECULES];
double 	epsilon[MAX_NUM_MOLECULES];

// default values
double 	box_x=0.0, box_y=0.0, box_z=0.0;
double 	verlet_cutoff=100.0;
int 	n_steps = 1000000;
int 	include_center_energy = 0;
int     check_insertion_point = 0;
int 	show_steps = 0;
double  drift_threshold = INFINITY;
int 	number_of_molecules = 0;
int     n_threads = 1;

double 	min_diameter = 0.0;
double 	characteristic_length = 1.0;
double 	characteristic_energy = 1.0;
int 	seed = 1;

int main(int argc, char **argv)
{
  setCommandLineParameters(argc, argv);
  getIntParam("-n_threads", &n_threads);
  getIntParam("-seed", &seed);
  getVectorParam("-box", &box_x, &box_y, &box_z);
  getDoubleParam("-characteristic_length", &characteristic_length);
  getDoubleParam("-characteristic_energy", &characteristic_energy);
  getDoubleParam("-verlet_cutoff", &verlet_cutoff);
  getDoubleParam("-drift_threshold", &drift_threshold);
  getIntParam("-n_steps", &n_steps);
  check_insertion_point = getFlagParam("-check_insertion_point");

  char *vertices_file;
  getStringParam("-vertices_file", &vertices_file);

  include_center_energy = getFlagParam("-include_center_energy");
  show_steps = getFlagParam("-show_steps");
  getDoubleParam("-min_diameter", &min_diameter);

  if (getFlagParam("-usage"))
  {
    printf("\n");
    printf("vddx uses a set of Voronoi vertices as the initial condition\n");
    printf("for particle insertion, generating an indexed list of cavities.\n");
    printf("The index is based on order of vertices read from vertices file.\n");
    printf("Output shows index and drift (distance from insertion point) for\n");
    printf("each of the cavities, with options to show number of steps taken\n");
    printf("and center energy as additional columns.\n\n");
    printf("vddx usage:     -box [ 0.0 0.0 0.0 ]\n");
    printf("                -seed [ 1 ]\n");
    printf("                -characteristic_length [ 1.0 ] if in doubt, use largest sigma/diameter\n");
    printf("                -characteristic_energy [ 1.0 ] if in doubt, use largest energy/epsilon\n");
    printf("                -n_steps [ 1000000 ] maximum before giving up\n");
    printf("                -n_threads [ 1 ] \n");
    printf("                -verlet_cutoff [ 100.0 ]\n");
    printf("                -vertices_file [ filename ] \n");
    printf("                -include_center_energy \n");
    printf("                -show_steps (includes steps taken as final column)\n");
    printf("                -check_insertion_point \n");
    printf("                -min_diameter [ 0.0 ] 0.0 will give all \n");
    printf("\n");
    printf("                IN:  read from stdin < .gfg as: %%f     %%f     %%f     %%f     %%f\n");
    printf("                                                 x      y      z      sigma  epsilon\n");
    printf("                OUT: write to stdout > .vcv as: %%d     %%f     %%f     %%f     %%f     %%f\n");
    printf("                                                idx     x      y      z      d      drift\n");
    printf("\n");
    exit(0);
  }

  if (box_x * box_y * box_z <= 0.0) 
  {
    fprintf(stderr, "vacuumms/vddx: please supply suitable value of box size. (%lf x %lf x %lf)\n", box_x, box_y, box_z);
    exit(2);
  }

  loadConfiguration();

  // load verts
  int n_vertices = loadVertices(vertices_file);

  // make and verify all the threads and resources
  passvals = (void **)malloc(sizeof(void*) * n_vertices);
  assert(passvals);
  threads = (pthread_t*)malloc(sizeof(pthread_t) * n_vertices);
  assert(threads);

  // set stack size for threads
  size_t stacksize = (size_t)2048;

  pthread_attr_t thread_attr;
  pthread_attr_init(&thread_attr);
  pthread_attr_setstacksize(&thread_attr, stacksize);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED);

  // initialize the semaphores
  int status;
  status = sem_init(&semaphore, 0, n_threads);
  assert(status == 0);
  sem_init(&completion_semaphore, 0, 0);
  assert(status == 0);
  int complete=0;

  /* This is the loop where all threads are started, wait and run */
  for (thread_idx=0; thread_idx<n_vertices; thread_idx++) {
    /* Since the only value passed is an integer, we just cast to and pass as void* */
    passvals[thread_idx] = (void*)(long)thread_idx; 
    sem_wait(&semaphore); // thread waits to become eligible
    int rc;
    rc = pthread_create(&threads[thread_idx], &thread_attr, &ThreadMain, passvals[thread_idx]);
    assert(rc == 0);
  }

  //  spinlock to wait for completion
  while(complete < n_vertices) sem_getvalue(&completion_semaphore, &complete);
  
  free(threads);
  free(passvals);

  return 0;
} // end main()

void *ThreadMain(void* passval) { 
  Trajectory *p_traj = (Trajectory *)malloc(sizeof(Trajectory));
  assert(p_traj);
  /* passval is just an int wrapped as void*, need to cast down to int via long to match type size */
  p_traj->thread_id = (int)(long)passval; 
  MersenneInitialize(&(p_traj->rng), seed + p_traj->thread_id);

  /* guarantee insertion point not be located in a valley whose bottom is > 0 energy... 
   * otherwise cavity cannot be sized with real value */
  // while (calculateEnergy(p_traj, 0.0) > 0) generateTestPoint(p_traj);

  generateTestPoint(p_traj);

  if (check_insertion_point && calculateEnergy(p_traj, 0.0) > 0) 
  {
    pthread_mutex_lock(&mutex);
    fprintf(stderr, "Insertion failing energy criterium at vertex %d\t%lf\t%lf\t%lf\t%lf with zero diameter energy %f\n", p_traj->thread_id, p_traj->test_x, p_traj->test_y, p_traj->test_z, p_traj->diameter, calculateEnergy(p_traj, min_diameter) );
    fflush(stderr);
    pthread_mutex_unlock(&mutex);
    
    free(p_traj);
    sem_post(&semaphore);
    sem_post(&completion_semaphore);
    pthread_exit(NULL);
  } 
   
  // now find energy minimum point
  double drift = findEnergyMinimum(p_traj);
   
  p_traj->sq_distance_from_initial_pt 	= (p_traj->test_x-p_traj->test_x0)*(p_traj->test_x-p_traj->test_x0) 
				   	+ (p_traj->test_y-p_traj->test_y0)*(p_traj->test_y-p_traj->test_y0) 
				   	+ (p_traj->test_z-p_traj->test_z0)*(p_traj->test_z-p_traj->test_z0);

  if (drift < drift_threshold)
  {
    makeVerletList(p_traj);
    expandTestParticle(p_traj);
    if (p_traj->diameter > min_diameter) {
      // correct for box edges...
      while (p_traj->test_x >= box_x) 	p_traj->test_x -= box_x;
      while (p_traj->test_x < 0) 	p_traj->test_x += box_x;
      while (p_traj->test_y >= box_y) 	p_traj->test_y -= box_y;
      while (p_traj->test_y < 0) 	p_traj->test_y += box_y;
      while (p_traj->test_z >= box_z) 	p_traj->test_z -= box_z;
      while (p_traj->test_z < 0) 	p_traj->test_z += box_z;

      pthread_mutex_lock(&mutex);
      printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf", p_traj->thread_id, 
                                            p_traj->test_x, 
                                            p_traj->test_y, 
                                            p_traj->test_z, 
                                            p_traj->diameter,
                                            drift);
      if (include_center_energy) printf("\t%lf", calculateEnergy(p_traj, p_traj->diameter));
      if (show_steps) printf("\t%d", p_traj->attempts);
      printf("\n");
      fflush(stdout);
      pthread_mutex_unlock(&mutex);
    }
  }
  
  free(p_traj);
  sem_post(&semaphore);
  sem_post(&completion_semaphore);
  pthread_exit(NULL);
}

void generateTestPoint(Trajectory *p_traj)
{
  p_traj->test_x = vertex_x[p_traj->thread_id];
  p_traj->test_y = vertex_y[p_traj->thread_id];
  p_traj->test_z = vertex_z[p_traj->thread_id];

  makeVerletList(p_traj);
}

void makeVerletList(Trajectory *p_traj)
{
  int i;
  double dx, dy, dz, dd;
  double shift_x, shift_y, shift_z;

  while (p_traj->test_x > box_x) p_traj->test_x -= box_x;
  while (p_traj->test_y > box_y) p_traj->test_y -= box_y;
  while (p_traj->test_z > box_z) p_traj->test_z -= box_z;

  while (p_traj->test_x < 0) p_traj->test_x += box_x;
  while (p_traj->test_y < 0) p_traj->test_y += box_y;
  while (p_traj->test_z < 0) p_traj->test_z += box_z;

  p_traj->verlet_center_x=p_traj->test_x;
  p_traj->verlet_center_y=p_traj->test_y;
  p_traj->verlet_center_z=p_traj->test_z;

  p_traj->close_molecules=0;
  for (i=0; i<number_of_molecules; i++)
  {
    for (shift_x = -box_x; shift_x <= box_x; shift_x += box_x)
    for (shift_y = -box_y; shift_y <= box_y; shift_y += box_y)
    for (shift_z = -box_z; shift_z <= box_z; shift_z += box_z)
    {
      dx = shift_x + x[i] - p_traj->test_x;
      dy = shift_y + y[i] - p_traj->test_y;
      dz = shift_z + z[i] - p_traj->test_z;

      dd = dx*dx + dy*dy + dz*dz;

      if (dd < verlet_cutoff) 
      { 
        p_traj->close_x[p_traj->close_molecules] = shift_x + x[i];
        p_traj->close_y[p_traj->close_molecules] = shift_y + y[i];
        p_traj->close_z[p_traj->close_molecules] = shift_z + z[i];
        p_traj->close_sigma[p_traj->close_molecules] = sigma[i];
        p_traj->close_sigma6[p_traj->close_molecules] = sigma[i]*sigma[i]*sigma[i]*sigma[i]*sigma[i]*sigma[i];
        p_traj->close_sigma12[p_traj->close_molecules] = p_traj->close_sigma6[p_traj->close_molecules]*p_traj->close_sigma6[p_traj->close_molecules];
        p_traj->close_epsilon[p_traj->close_molecules] = epsilon[i];

        p_traj->close_molecules++;
        assert(p_traj->close_molecules < MAX_CLOSE);
      }
    }
  }
}

// return value is the drift from insertion point
double findEnergyMinimum(Trajectory *p_traj)
{
  double dx, dy, dz, dd, d6, d14;
  double factor;
  double old_energy;
  double new_energy;
  double grad_x, grad_y, grad_z;
  double step_x, step_y, step_z;
  int i;
  double drift_sq;
  //int attempts;

  makeVerletList(p_traj);

  double cumulative_drift_x = 0.0;
  double cumulative_drift_y = 0.0;
  double cumulative_drift_z = 0.0;
  double cumulative_drift_sq;

  double drift_x;
  double drift_y;
  double drift_z;

  // begin loop to iterate until minimum found
  for (p_traj->attempts=0; p_traj->attempts<n_steps; p_traj->attempts++)
  {
    drift_x = (p_traj->test_x-p_traj->verlet_center_x);
    drift_y = (p_traj->test_y-p_traj->verlet_center_y);
    drift_z = (p_traj->test_z-p_traj->verlet_center_z);
    drift_sq = drift_x * drift_x + drift_y * drift_y + drift_z * drift_z; 

    if (drift_sq > .01 * verlet_cutoff)
    {
      makeVerletList(p_traj);
      cumulative_drift_x += drift_x;
      cumulative_drift_y += drift_y;
      cumulative_drift_z += drift_z;
    }

    // find the gradient at test_x, test_y, test_Z using the derivative of energy
    grad_x=0; grad_y=0; grad_z=0;

    for (i=0; i<p_traj->close_molecules; i++)
    {
      dx = p_traj->test_x - p_traj->close_x[i];
      dy = p_traj->test_y - p_traj->close_y[i];
      dz = p_traj->test_z - p_traj->close_z[i];
      dd = dx*dx + dy*dy + dz*dz;
      d6 = dd*dd*dd;
      d14 = d6*d6*dd;

      // The analytical expression for the gradient contribution contains a factor of -48.0.
      // The minus is reflected in the sense of the step taken.  The factor of 48 is factored out in the normalization.
      factor = p_traj->close_epsilon[i] * p_traj->close_sigma12[i] / d14;

      grad_x += dx * factor;
      grad_y += dy * factor;
      grad_z += dz * factor;
    }

    double grad_sq = grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
    if (grad_sq < 1.0e-12) break; // If the gradient reaches zero, not going anywhere...
    // if modulus of gradient is > characteristic_energy/characteristic_length, normalize the gradient (too big means unstable)
    // gradient modulus is guaranteed <= characteristic_energy / characteristic_length
    if (grad_sq * characteristic_length * characteristic_length > characteristic_energy * characteristic_energy){
      double grad_modulus = sqrt(grad_sq);
      grad_x /= grad_modulus;
      grad_y /= grad_modulus;
      grad_z /= grad_modulus;
    }

    // should be three reasons to break loop... either the gradient is zero -or- some tolerance test -or- number of steps exceeded.

    old_energy = calculateRepulsion(p_traj);
    step_x = grad_x * characteristic_length;
    step_y = grad_y * characteristic_length;
    step_z = grad_z * characteristic_length;   
    p_traj->test_x += step_x;
    p_traj->test_y += step_y;
    p_traj->test_z += step_z;

    // if new point has higher energy, keep working back midway to a new point with lower energy
    double step_size_factor = 0.5;
    for(new_energy = calculateRepulsion(p_traj); new_energy > old_energy; step_size_factor *= 0.5) {
      p_traj->test_x -= step_size_factor * step_x;
      p_traj->test_y -= step_size_factor * step_y;
      p_traj->test_z -= step_size_factor * step_z;
      new_energy = calculateRepulsion(p_traj);
      if (step_size_factor == 0) break;
    }
  } // end for atempts

  // No more moves. Add in the final drift since last verlet update:
  cumulative_drift_x += drift_x;
  cumulative_drift_y += drift_y;
  cumulative_drift_z += drift_z;

  cumulative_drift_sq = cumulative_drift_x * cumulative_drift_x
                      + cumulative_drift_y * cumulative_drift_y
                      + cumulative_drift_z * cumulative_drift_z;

  return sqrt(cumulative_drift_sq);
}

double calculateRepulsion(Trajectory *p_traj)
{
  double repulsion=0;
  double dx, dy, dz, dd, d6, d12;
  int i;

  for (i=0; i<p_traj->close_molecules; i++)
  {
    dx = p_traj->close_x[i] - p_traj->test_x;
    dy = p_traj->close_y[i] - p_traj->test_y;
    dz = p_traj->close_z[i] - p_traj->test_z;
    dd = dx*dx + dy*dy + dz*dz;
    d6 = dd*dd*dd;
    d12 = d6*d6;

    repulsion += p_traj->close_epsilon[i] * p_traj->close_sigma12[i] / d12;
  }
 
  return 4.0 * repulsion;
}

double calculateEnergy(Trajectory *p_traj, double test_diameter)
{
  double repulsion=0;
  double attraction=0;
  double dx, dy, dz, dd, d6, d12;
  double sigma, sigma6, sigma12;
  int i;

  for (i=0; i<p_traj->close_molecules; i++)
  {
    dx = p_traj->close_x[i] - p_traj->test_x;
    dy = p_traj->close_y[i] - p_traj->test_y;
    dz = p_traj->close_z[i] - p_traj->test_z;
    dd = dx*dx + dy*dy + dz*dz;
    d6 = dd*dd*dd;
    d12 = d6*d6;

    sigma = 0.5 * (p_traj->close_sigma[i] + test_diameter);
    sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
    sigma12 = sigma6*sigma6;

    repulsion += p_traj->close_epsilon[i] * sigma12/d12;
    attraction += p_traj->close_epsilon[i] * sigma6/d6;
  }

  return 4.0 * (repulsion - attraction);
}

void expandTestParticle(Trajectory *p_traj)
{
  double slope;
  double step_size;
  double energy, old_energy;
  double e0, e1, r0, r1;

  double diameter = min_diameter; 

  // improved initial guess
  old_energy = calculateEnergy(p_traj, diameter);
  if (old_energy > 0) return;
  while (diameter += .1)
  {
    energy = calculateEnergy(p_traj, diameter);
    if (energy > old_energy) break;
    old_energy = energy;
  }

  while(1) // Newton's method
  {
    r0 = diameter - .001;
    r1 = diameter + .001;
    
    e0 = calculateEnergy(p_traj, r0);
    e1 = calculateEnergy(p_traj, r1);
    energy = calculateEnergy(p_traj, diameter);

    slope = (e1-e0)/(r1-r0);
    step_size = -energy/slope;

    diameter = diameter + step_size;

    if (step_size*step_size < .00000001) break;
  }

  // copy diameter to trajectory data
  p_traj->diameter = diameter;
}

