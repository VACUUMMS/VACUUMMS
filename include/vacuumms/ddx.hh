/* vacuumms/ddx.hh */

#pragma once

#include <vacuumms/limits.h>

#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/parameters.hh>
#include <vacuumms/operations.hh>
#include <vacuumms/exports.hh>


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
        void setSeed(int);
        void setVerletCutoff(vacuumms_float);
        void setNumberOfSteps(int);
        void setMinDiameter(vacuumms_float);
        void Randomize();
        static void printUsage();
        CavityConfiguration getResult();

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

    private:

        double calculateRepulsion();
        double calculateEnergy(double test_diameter);

        void generateTestPoint();
        void findEnergyMinimum();
        void makeVerletList();
        void expandTestParticle();

        CavityConfiguration result;
        Configuration c;
        Parameters p;

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

        double box_x=0.0, box_y=0.0, box_z=0.0;
        double verlet_cutoff=100.0;

        //double step_size_factor = 1.0;
        int n_steps = 1000;

        int number_of_samples = 1;
        int volume_sampling = 0;
        int include_center_energy = 0;
        int show_steps = 0;

        double test_x0, test_y0, test_z0;
        double test_x, test_y, test_z;
        double verlet_center_x, verlet_center_y, verlet_center_z;
        double diameter = 1.0;
        double min_diameter = 0.0;
        double characteristic_length = 1.0;
        double characteristic_energy = 1.0;
        double precision_parameter = 0.001; // decimal 
        int seed = 1;

        int number_of_molecules = 0;
        int close_molecules;
        int attempts;

        FILE *instream;

        int verbose;

}; // end class DDX


