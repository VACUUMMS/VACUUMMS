/* vacuumms/ddx.hh */

#pragma once

// Required for sem_init on macOS/Linux
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 
#elif _XOPEN_SOURCE < 600
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif

#include <thread>
#include <semaphore.h>
#include <mutex>

//FTW can i get rid of this yet?
#include <vacuumms/limits.h>

#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/parameters.hh>
#include <vacuumms/operations.hh>
#include <vacuumms/exports.hh>
#include <vacuumms/rng.hh>


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
DDX : public Operation
{
    public:

        DDX(Configuration c, Parameters p);
        DDX(Configuration c);
        DDX();

        void execute();
        void setParameters(Parameters p);
        void setConfiguration(Configuration c);
        Configuration getConfiguration();
        void setNumberOfSamples(int);
        void setVerletCutoff(vacuumms_float);
        void setVerletExtent(int);
        void setNumberOfSteps(int);
        void setMinDiameter(vacuumms_float);
        void setLearningRate(vacuumms_float);
        void setTolerance(vacuumms_float);
        static void printUsage();
        CavityConfiguration getResult();
        void reorderResults();
        void setNumberOfThreads(int);
        void setRNGSeed(int);

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

    private:

        // Results and outer settings
        CavityConfiguration results;
        std::vector<int> results_order;
        Configuration configuration;
        Parameters parameters;

        // settable parameters
        vacuumms_float verlet_cutoff=100.0;
        int verlet_extent = 1;
        int number_of_samples = 1;
        int number_of_steps = 100;
        int rng_seed = 1;
        vacuumms_float min_diameter = 0.0;
        vacuumms_float learning_rate = 0.01f;
        vacuumms_float tolerance = 10.0f;
        int volume_sampling = 0;

        // Concurrency stuff
        int number_of_threads = 1;
        std::vector<std::thread> threads;
        sem_t semaphore;
        std::mutex results_mutex;
        void run_sample(int);

        vacuumms_float box_x=0.0, box_y=0.0, box_z=0.0; // vals pulled from Configuration c

        class Sample
        {
            // per-sample members
            vacuumms_float test_x0, test_y0, test_z0;
            vacuumms_float test_x, test_y, test_z;
            vacuumms_float verlet_center_x, verlet_center_y, verlet_center_z;
            vacuumms_float diameter = 0.0;
            Configuration verlet_list;
            MersenneTwister rng;
            pthread_t thread;
            DDX* outer;
            int id;

        public:
            Sample(DDX* _outer, int _id);

            // Per-sample operations
            static void* start(void* arg);
            vacuumms_float calculateRepulsion();
            vacuumms_float calculateEnergy(vacuumms_float test_diameter);
            void generateTestPoint();
            void findEnergyMinimum();
            void makeVerletList();
            void expandTestParticle();
        };

}; // end class DDX

