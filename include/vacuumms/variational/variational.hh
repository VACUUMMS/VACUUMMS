// variational.hh
#include "configuration.hh"
#include <vacuumms/types.h>

extern const vacuumms_float machine_epsilon;
extern const vacuumms_float sqrt_machine_epsilon;

class Variational2D
{
    private:

        vacuumms_float start_x; 
        vacuumms_float start_y; 
        vacuumms_float end_x; 
        vacuumms_float end_y; 
        vacuumms_float sigma; 
        vacuumms_float epsilon; 
    
        int n_var_points;
        Configuration *configuration;
        vacuumms_float(*energy_function)(vacuumms_float x, vacuumms_float y);

        bool use_configuration_energy;

        vacuumms_float* var_x;
        vacuumms_float* var_y;

    public:

        void init(vacuumms_float start_x, 
             vacuumms_float start_y, 
             vacuumms_float end_x, 
             vacuumms_float end_y, 
             vacuumms_float sigma, 
             vacuumms_float epsilon, 
             int n_var_points);

        Variational2D(vacuumms_float start_x, 
                      vacuumms_float start_y, 
                      vacuumms_float end_x, 
                      vacuumms_float end_y, 
                      vacuumms_float sigma, 
                      vacuumms_float epsilon, 
                      int n_var_points, 
                      vacuumms_float(*energy_function)(vacuumms_float x, vacuumms_float y));

        Variational2D(vacuumms_float start_x, 
                      vacuumms_float start_y, 
                      vacuumms_float end_x, 
                      vacuumms_float end_y, 
                      vacuumms_float sigma, 
                      vacuumms_float epsilon, 
                      int n_var_points, 
                      Configuration *c);

        void printValues();
        void iterate();
        void iterateWork();

        ~Variational2D();

        vacuumms_float rebalancePoints2D();

};

class Variational3D
{
    private:

        // step size regulators
        vacuumms_float alpha = 0.01;
        vacuumms_float alpha_max = 0.01;
        vacuumms_float beta = 1.1; // alpha adjustment parameter
        vacuumms_float delta_max = 0.01; // maximum step size
    
        vacuumms_float start_x; 
        vacuumms_float start_y; 
        vacuumms_float start_z; 
        vacuumms_float end_x; 
        vacuumms_float end_y; 
        vacuumms_float end_z; 
        int n_var_points;

        // parameters for inserted particle
        vacuumms_float sigma; 
        vacuumms_float epsilon;

        Configuration *configuration;
        vacuumms_float(*energy_function)(vacuumms_float x, vacuumms_float y, vacuumms_float z);

//        vacuumms_float respaceKernel(/*int _n_var_points,*/ vacuumms_float _start_x, vacuumms_float _start_y, vacuumms_float _start_z, 
        int respaceKernel(vacuumms_float _start_x, vacuumms_float _start_y, vacuumms_float _start_z, 
                          vacuumms_float _end_x, vacuumms_float _end_y, vacuumms_float _end_z, vacuumms_float _var_x[], vacuumms_float _var_y[], vacuumms_float _var_z[], 
                          vacuumms_float _new_var_x[], vacuumms_float _new_var_y[], vacuumms_float _new_var_z[]);

        vacuumms_float calculateCurveLength(vacuumms_float _start_x, vacuumms_float _start_y, vacuumms_float _start_z, vacuumms_float _end_x, vacuumms_float _end_y, vacuumms_float _end_z, vacuumms_float _var_x[], vacuumms_float _var_y[], vacuumms_float _var_z[]);

        bool use_configuration_energy;

        vacuumms_float* var_x;
        vacuumms_float* var_y;
        vacuumms_float* var_z;

    public:

        void init(vacuumms_float start_x, 
             vacuumms_float start_y, 
             vacuumms_float start_z, 
             vacuumms_float end_x, 
             vacuumms_float end_y, 
             vacuumms_float end_z, 
             vacuumms_float sigma,
             vacuumms_float epsilon,
             int n_var_points);

/*
        Variational3D(vacuumms_float start_x, 
                      vacuumms_float start_y, 
                      vacuumms_float start_z, 
                      vacuumms_float end_x, 
                      vacuumms_float end_y, 
                      vacuumms_float end_z, 
                      int n_var_points, 
                      vacuumms_float(*energy_function)(vacuumms_float x, vacuumms_float y, vacuumms_float z));
*/

        Variational3D(vacuumms_float start_x, 
                      vacuumms_float start_y, 
                      vacuumms_float start_z, 
                      vacuumms_float end_x, 
                      vacuumms_float end_y, 
                      vacuumms_float end_z, 
                      vacuumms_float sigma,
                      vacuumms_float epsilon,
                      int n_var_points, 
                      Configuration *c);

        vacuumms_float* getX();
        vacuumms_float* getY();
        vacuumms_float* getZ();

        void setAlpha(vacuumms_float _alpha);
        void setAlphaMax(vacuumms_float _alpha_max);
        void setDeltaMax(vacuumms_float _delta_max);
        void printValues();
        void iterate();
        void iterateWork();
        vacuumms_float adaptiveIterateAndUpdate();

        ~Variational3D();

        vacuumms_float rebalancePoints3D();

};
