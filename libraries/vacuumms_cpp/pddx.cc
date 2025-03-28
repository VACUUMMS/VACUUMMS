/* vacuumms_cpp/pddx.cc */

#include <vacuumms/pddx.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>

#include <vacuumms/limits.h>

#include <math.h>
#include <pthread.h>

#include <cassert>


PDDX::PDDX(Configuration c, Parameters p) : 
    c{c}, p{p} 
{
}

PDDX::PDDX()
{
}


void PDDX::setParameters(Parameters _p)
{
    p = _p;
}

void PDDX::setConfiguration(Configuration _c)
{
    c = _c;
}

Configuration PDDX::getConfiguration()
{
    return c;
}

#ifdef BUILD_PYBIND_BINDINGS
pybind11::str PDDX::__repr__()
{
    pybind11::str retval;
    retval += c.__repr__();
    retval += p.__repr__();
    retval += result.__repr__();
    return retval;
}
#endif

CavityConfiguration PDDX::getResult()
{
    return result;   
}

void PDDX::printUsage()
{
    printf("\nPDDX usage:\t-box [ 6.0 6.0 6.0 ]\n");
    printf("\t\t-seed [ 1 ]\n");
    printf("\t\t-n_threads[ 1 ]\n");
    printf("\t\t-randomize \n");
    printf("\t\t-characteristic_length [ 1.0 ]\n");
    printf("\t\t-characteristic_energy [ 1.0 ]\n");
    printf("\t\t-precision_parameter [ 0.001 ]\n");
    printf("\t\t-n_steps [ 1000 ] (roughly reciprocal of precision parameter)\n");
    printf("\t\t-show_steps (includes steps taken as final column)\n");
    printf("\t\t-verlet_cutoff [ 100.0 ]\n");
    printf("\t\t-n [ 1 ]\n");
    printf("\t\t-volume_sampling \n");
    printf("\t\t-include_center_energy \n");
    printf("\t\t-min_diameter [ 0.0 ]");
    printf("\n");
}


/* Helper struct and function for threading of member function */

struct ThreadArgs
{
    PDDX* calling_instance;
    int thread_id;
};

void* threadEntry(void* arg)
{
    ThreadArgs* ta = static_cast<ThreadArgs*>(arg);
    ta->calling_instance->ThreadMain((void*)(long)(ta->thread_id));
//    delete ta;
    return nullptr;
}


// set up, create, run, and join threads
void PDDX::execute()
{
  double sq_distance_from_initial_pt;

  verbose = p.getFlagParam((char*)"-verbose");
  p.getIntParam((char*)"-seed", &seed);
  p.getIntParam((char*)"-n_threads", &n_threads);
    // replaced rng with prng/Mersenne
  //if (p.getFlagParam((char*)"-randomize")) randomize();
  //else initializeRandomNumberGeneratorTo(seed);

  p.getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
  p.getDoubleParam((char*)"-characteristic_length", &characteristic_length);
  p.getDoubleParam((char*)"-characteristic_energy", &characteristic_energy);
  p.getDoubleParam((char*)"-precision_parameter", &precision_parameter);
  p.getDoubleParam((char*)"-verlet_cutoff", &verlet_cutoff);
  p.getIntParam((char*)"-n", &number_of_samples);
  p.getIntParam((char*)"-n_steps", &n_steps);
  volume_sampling = p.getFlagParam((char*)"-volume_sampling");
  include_center_energy = p.getFlagParam((char*)"-include_center_energy");
  show_steps = p.getFlagParam((char*)"-show_steps");
  p.getDoubleParam((char*)"-min_diameter", &min_diameter);
  if (p.getFlagParam((char*)"-usage")) printUsage();

  // load configuration
  
  for (int i=0; i<c.getSize(); i++) 
  {
    ConfigurationRecord r = c.recordAt(i);
    x[i] = r.x;
    y[i] = r.y;
    z[i] = r.z;
    sigma[i] = r.sigma;
    epsilon[i] = r.epsilon;
  }
  
  number_of_molecules = c.getSize();
  
  // make and verify all the threads and resources
  
  threads = (pthread_t*)malloc(sizeof(pthread_t) * number_of_samples);
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
  for (thread_idx=0; thread_idx<number_of_samples; thread_idx++) {
    sem_wait(&semaphore); // thread waits to become eligible
    int rc;

    ThreadArgs passval = {this, thread_idx};
    rc = pthread_create(&threads[thread_idx], &thread_attr, threadEntry, &passval);

    assert(rc == 0);
  }

  //  spinlock to wait for completion
  while(complete < number_of_samples) sem_getvalue(&completion_semaphore, &complete);

  free(threads);

} // end execute()
 

//------
// replaces execute() of ddx serial version
//void *PDDX::ThreadMain(void *threadID)
//{
//------
//
void *PDDX::ThreadMain(void* passval) 
{
    Trajectory *p_traj = (Trajectory*)malloc(sizeof(Trajectory));
    assert(p_traj);
    /* passval is just an int wrapped as void*, need to cast down to int via long to match type size */
    p_traj->thread_id = (int)(long)passval;
    MersenneInitialize(&(p_traj->rng), seed + p_traj->thread_id);

    generateTestPoint(p_traj);

    while (calculateEnergy(p_traj, 0.0) > 0) generateTestPoint(p_traj);
    
    findEnergyMinimum(p_traj);
    
    p_traj->sq_distance_from_initial_pt = 
        (p_traj->test_x-p_traj->test_x0)*(p_traj->test_x-p_traj->test_x0) + 
        (p_traj->test_y-p_traj->test_y0)*(p_traj->test_y-p_traj->test_y0) + 
        (p_traj->test_z-p_traj->test_z0)*(p_traj->test_z-p_traj->test_z0);
    if (!volume_sampling || (p_traj->sq_distance_from_initial_pt < .25 * p_traj->diameter * p_traj->diameter))
    {
      makeVerletList(p_traj);
      expandTestParticle(p_traj);
      if (p_traj->diameter > min_diameter) 
      {

        // correct for box edges...
        while (p_traj->test_x >= box_x) p_traj->test_x -= box_x;
        while (p_traj->test_x < 0) p_traj->test_x += box_x;
        while (p_traj->test_y >= box_y) p_traj->test_y -= box_y;
        while (p_traj->test_y < 0) p_traj->test_y += box_y;
        while (p_traj->test_z >= box_z) p_traj->test_z -= box_z;
        while (p_traj->test_z < 0) p_traj->test_z += box_z;

//        printf("%lf\t%lf\t%lf\t%lf", test_x, test_y, test_z, diameter);
        result.pushBack(Cavity(p_traj->test_x, p_traj->test_y, p_traj->test_z, p_traj->diameter));
//        if (include_center_energy) printf("\t%lf", calculateEnergy(diameter));
//        if (show_steps) printf("\t%d", attempts);
//        printf("\n");
//        number_of_samples--;
      }
    }

free(p_traj);
  sem_post(&semaphore);
  sem_post(&completion_semaphore);
  pthread_exit(NULL);

    return nullptr;

} // end PDDX::ThreadMain()


void PDDX::generateTestPoint(Trajectory* p_traj)
{
  p_traj->test_x = p_traj->test_x0 = prnd(&(p_traj->rng)) * box_x;
  p_traj->test_y = p_traj->test_y0 = prnd(&(p_traj->rng)) * box_y;
  p_traj->test_z = p_traj->test_z0 = prnd(&(p_traj->rng)) * box_z;

  makeVerletList(p_traj);

} // end PDDX::generateTestPoint()

void PDDX::makeVerletList(Trajectory* p_traj)
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
} // end PDDX::makeVerletList()

void PDDX::findEnergyMinimum(Trajectory *p_traj)
{
  double dx, dy, dz, dd, d6, d14;
  double factor;
  double old_energy;
  double new_energy;
  double grad_x, grad_y, grad_z;
  double step_x, step_y, step_z;
  int i;
  double drift_sq;

  makeVerletList(p_traj);

  // begin loop to iterate until minimum found
  for (p_traj->attempts=0; p_traj->attempts<n_steps; p_traj->attempts++)
  {
    drift_sq = (p_traj->test_x-p_traj->verlet_center_x)*(p_traj->test_x-p_traj->verlet_center_x) 
             + (p_traj->test_y-p_traj->verlet_center_y)*(p_traj->test_y-p_traj->verlet_center_y) 
             + (p_traj->test_z-p_traj->verlet_center_z)*(p_traj->test_z-p_traj->verlet_center_z);

    if (drift_sq > .01 * verlet_cutoff) makeVerletList(p_traj);

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

    // normalize the gradient
    double grad_sq = grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
    double grad_modulus = sqrt(grad_sq);
    grad_x /= grad_modulus;
    grad_y /= grad_modulus;
    grad_z /= grad_modulus;

    old_energy = calculateRepulsion(p_traj);
//printf("xyz:: %12lf\t%12lf\t%12lf ::%12lf::\t  %12lf\t%12lf\t%12lf\n", test_x, test_y, test_z, old_energy, grad_x, grad_y, grad_z);

    step_x = grad_x * characteristic_energy * characteristic_length; while (step_x * step_x > characteristic_length * characteristic_length * precision_parameter * precision_parameter) {step_x *=.5;}
    step_y = grad_y * characteristic_energy * characteristic_length; while (step_y * step_y > characteristic_length * characteristic_length * precision_parameter * precision_parameter) {step_y *=.5;}
    step_z = grad_z * characteristic_energy * characteristic_length; while (step_z * step_z > characteristic_length * characteristic_length * precision_parameter * precision_parameter) {step_z *=.5;}

//    removed this criteria for assessing minima.... step size no longer shrinks because gradient is now normalized
//    step_sq = step_x * step_x + step_y * step_y + step_z * step_z;
//    if (step_sq < (characteristic_length * characteristic_length * precision_parameter * precision_parameter)) break;  // close enough, exit loop

    p_traj->test_x += step_x;
    p_traj->test_y += step_y;
    p_traj->test_z += step_z;
 
    // check repulsion at new location
    new_energy = calculateRepulsion(p_traj);
    // if the energy fluctuates up by a fraction of the characteristic energy, call it
    if (new_energy - old_energy > precision_parameter * characteristic_energy)
    {
//      printf("exiting from energy increase %lf\n", new_energy);
      break;
    }
  }
} // end PDDX::findEnergyMinimum()

double PDDX::calculateRepulsion(Trajectory* p_traj)
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
} // end PDDX::calculateRepulsion()

double PDDX::calculateEnergy(Trajectory* p_traj, double test_diameter)
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
} // end PDDX::calculateEnergy()

void PDDX::expandTestParticle(Trajectory* p_traj)
{
  double slope;
  double step_size;
  double energy, old_energy;
  double e0, e1, r0, r1;

  vacuumms_float diameter = 0.0;

  // improved initial guess
  old_energy = calculateEnergy(p_traj, diameter);
//printf("energy at sigma = 0:  %lf\n", old_energy);
  if (old_energy > 0) return;
  while (diameter += .1)
  {
    energy = calculateEnergy(p_traj, diameter);
    if (energy > old_energy) break;
    old_energy = energy;
  }

//printf("starting w/diameter=%lf\n", diameter);
    
  // Newton's method

  while(1)
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

    p_traj->diameter = diameter;
  }
} // end PDDX::expandTestParticle()

