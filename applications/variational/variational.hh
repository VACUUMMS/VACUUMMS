// variational.hh
#include "configuration.hh"

extern const float machine_epsilon;
extern const float sqrt_machine_epsilon;

class Variational2D
{
    private:

        float start_x; 
        float start_y; 
        float end_x; 
        float end_y; 
        int n_var_points;
        Configuration *configuration;
        float(*energy_function)(float x, float y);

        bool use_configuration_energy;

        float* var_x;
        float* var_y;

    public:

        void init(float start_x, 
             float start_y, 
             float end_x, 
             float end_y, 
             int n_var_points);

        Variational2D(float start_x, 
                      float start_y, 
                      float end_x, 
                      float end_y, 
                      int n_var_points, 
                      float(*energy_function)(float x, float y));

        Variational2D(float start_x, 
                      float start_y, 
                      float end_x, 
                      float end_y, 
                      int n_var_points, 
                      Configuration *c);

        void printValues();
        void iterate();
        void iterateWork();

        ~Variational2D();

        float rebalancePoints2D();

};

class Variational3D
{
    private:

        // step size regulators
        float alpha = 0.01;
        float alpha_max = 0.01;
        float beta = 1.1; // alpha adjustment parameter
        float delta_max = 0.01; // maximum step size
    
        float start_x; 
        float start_y; 
        float start_z; 
        float end_x; 
        float end_y; 
        float end_z; 
        int n_var_points;

        // parameters for inserted particle
        float sigma; 
        float epsilon;

        Configuration *configuration;
        float(*energy_function)(float x, float y, float z);

//        float respaceKernel(/*int _n_var_points,*/ float _start_x, float _start_y, float _start_z, 
        int respaceKernel(float _start_x, float _start_y, float _start_z, 
                          float _end_x, float _end_y, float _end_z, float _var_x[], float _var_y[], float _var_z[], 
                          float _new_var_x[], float _new_var_y[], float _new_var_z[]);

        float calculateCurveLength(float _start_x, float _start_y, float _start_z, float _end_x, float _end_y, float _end_z, float _var_x[], float _var_y[], float _var_z[]);

        bool use_configuration_energy;

        float* var_x;
        float* var_y;
        float* var_z;

    public:

        void init(float start_x, 
             float start_y, 
             float start_z, 
             float end_x, 
             float end_y, 
             float end_z, 
             float sigma,
             float epsilon,
             int n_var_points);

/*
        Variational3D(float start_x, 
                      float start_y, 
                      float start_z, 
                      float end_x, 
                      float end_y, 
                      float end_z, 
                      int n_var_points, 
                      float(*energy_function)(float x, float y, float z));
*/

        Variational3D(float start_x, 
                      float start_y, 
                      float start_z, 
                      float end_x, 
                      float end_y, 
                      float end_z, 
                      float sigma,
                      float epsilon,
                      int n_var_points, 
                      Configuration *c);

        float* getX();
        float* getY();
        float* getZ();

        void setAlpha(float _alpha);
        void setAlphaMax(float _alpha_max);
        void setDeltaMax(float _delta_max);
        void printValues();
        void iterate();
        void iterateWork();
        float adaptiveIterateAndUpdate();

        ~Variational3D();

        float rebalancePoints3D();

};
