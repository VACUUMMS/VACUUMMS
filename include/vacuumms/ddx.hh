/* vacuumms/ddx.hh */

#pragma once

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
        void setRNGSeed(int);
        static void printUsage();
        CavityConfiguration getResult();

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

    private:

        vacuumms_float calculateRepulsion();
        vacuumms_float calculateEnergy(vacuumms_float test_diameter);

        void generateTestPoint();
        void findEnergyMinimum();
        void makeVerletList();
        void expandTestParticle();

        CavityConfiguration result;
        Configuration configuration;
        Configuration verlet_list;
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

/*
        // Working vars from C implementation
        vacuumms_float x[VACUUMMS_MAX_NUMBER_OF_MOLECULES];
        vacuumms_float y[VACUUMMS_MAX_NUMBER_OF_MOLECULES];
        vacuumms_float z[VACUUMMS_MAX_NUMBER_OF_MOLECULES];
        vacuumms_float sigma[VACUUMMS_MAX_NUMBER_OF_MOLECULES];
        vacuumms_float epsilon[VACUUMMS_MAX_NUMBER_OF_MOLECULES];

        vacuumms_float close_x[VACUUMMS_MAX_CLOSE], close_y[VACUUMMS_MAX_CLOSE], close_z[VACUUMMS_MAX_CLOSE];
        vacuumms_float close_sigma[VACUUMMS_MAX_CLOSE];
        vacuumms_float close_sigma6[VACUUMMS_MAX_CLOSE];
        vacuumms_float close_sigma12[VACUUMMS_MAX_CLOSE];
        vacuumms_float close_epsilon[VACUUMMS_MAX_CLOSE];
*/
        vacuumms_float box_x=0.0, box_y=0.0, box_z=0.0; // vals pulled from Configuration c


        // Operating parameters
        vacuumms_float test_x0, test_y0, test_z0;
        vacuumms_float test_x, test_y, test_z;
        vacuumms_float verlet_center_x, verlet_center_y, verlet_center_z;
        vacuumms_float diameter = 10.0;

        MersenneTwister rng;

        int number_of_molecules = 0;
        int close_molecules;

}; // end class DDX

