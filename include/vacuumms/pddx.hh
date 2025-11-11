/* vacuumms/pddx.hh */

#pragma once

#include <semaphore.h>

#include <vacuumms/limits.h>
#include <vacuumms/types.h>

#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/parameters.hh>
#include <vacuumms/operations.hh>
//#include <vacuumms/prng.hh>

#include <vacuumms/exports.hh>


void* threadEntry(void* arg);

class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
PDDX : public Operation
{
    public:

        PDDX(Configuration c, Parameters p);
        PDDX();
        void execute();
//        CavityConfiguration getOutput();
/* gutting this so only constructor remains, to provide warning/redirect to DDX
        void setParameters(Parameters p);
        void setConfiguration(Configuration c);
        Configuration getConfiguration();
        static void printUsage();
        CavityConfiguration getResult();

        // This needs to be public so helper function can access it... ick.
        void *ThreadMain(void *threadID);

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

*/
    private:

        // All the data unique to a particular sample is in this struct */
/*
        typedef struct 
        {
            int                       thread_id;
            int                       close_molecules;
            int                       attempts;
            double                    test_x0, test_y0, test_z0;
            double                    test_x, test_y, test_z;
            double                    verlet_center_x, verlet_center_y, verlet_center_z;
            double                    diameter;
            double                    close_x[VACUUMMS_MAX_CLOSE], close_y[VACUUMMS_MAX_CLOSE], close_z[VACUUMMS_MAX_CLOSE];
            double                    close_sigma[VACUUMMS_MAX_CLOSE];
            double                    close_sigma6[VACUUMMS_MAX_CLOSE];
            double                    close_sigma12[VACUUMMS_MAX_CLOSE];
            double                    close_epsilon[VACUUMMS_MAX_CLOSE];
            double                    sq_distance_from_initial_pt;
            struct MersenneTwister    rng;
        } Trajectory;

        double calculateRepulsion(Trajectory*);
        double calculateEnergy(Trajectory*, double test_diameter);
        void generateTestPoint(Trajectory*);
        void findEnergyMinimum(Trajectory*);
        void makeVerletList(Trajectory*);
        void expandTestParticle(Trajectory*);

        CavityConfiguration result;
*/
        Configuration c;
        Parameters p;
/*
        // Working vars from C implementation
        double x[VACUUMMS_MAX_NUMBER_OF_MOLECULES];
        double y[VACUUMMS_MAX_NUMBER_OF_MOLECULES];
        double z[VACUUMMS_MAX_NUMBER_OF_MOLECULES];
        double sigma[VACUUMMS_MAX_NUMBER_OF_MOLECULES];
        double epsilon[VACUUMMS_MAX_NUMBER_OF_MOLECULES];

        double close_x[VACUUMMS_MAX_CLOSE], close_y[VACUUMMS_MAX_CLOSE], close_z[VACUUMMS_MAX_CLOSE];
        double close_sigma[VACUUMMS_MAX_CLOSE];
        double close_sigma6[VACUUMMS_MAX_CLOSE];
        double close_sigma12[VACUUMMS_MAX_CLOSE];
        double close_epsilon[VACUUMMS_MAX_CLOSE];

        double box_x=6, box_y=6, box_z=6;
        double verlet_cutoff=100.0;

        //double step_size_factor = 1.0;
        int n_steps = 1000;
        int n_threads = 1;

        // pthreads stuff
        pthread_t*          threads;    // the threads
        void**              passvals;   // values passed to each thread
        sem_t               semaphore;  // semaphore to restrict number of threads running
        sem_t               completion_semaphore;   // semaphore to count # of completed threads
        pthread_mutex_t     mutex;
        int                 thread_idx;

        int number_of_samples = 1;
        int volume_sampling = 0;
        int include_center_energy = 0;
        int show_steps = 0;

// These vars are all localized to Trajectory thread
//        double test_x0, test_y0, test_z0;
//        double test_x, test_y, test_z;
//        double verlet_center_x, verlet_center_y, verlet_center_z;
//        double diameter = 1.0;
//        int close_molecules;
//        int attempts;

        double min_diameter = 0.0;
        double characteristic_length = 1.0;
        double characteristic_energy = 1.0;
        double precision_parameter = 0.001; // decimal 
        int seed = 1;

        int number_of_molecules = 0;

        FILE *instream;

        int verbose;
*/

}; // end class PDDX

