// variational.hh
#include "configuration.hh"

extern const float epsilon;
extern const float sqrt_epsilon;

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

        float start_x; 
        float start_y; 
        float start_z; 
        float end_x; 
        float end_y; 
        float end_z; 
        int n_var_points;
        Configuration *configuration;
        float(*energy_function)(float x, float y, float z);

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
             int n_var_points);

        Variational3D(float start_x, 
                      float start_y, 
                      float start_z, 
                      float end_x, 
                      float end_y, 
                      float end_z, 
                      int n_var_points, 
                      float(*energy_function)(float x, float y, float z));

        Variational3D(float start_x, 
                      float start_y, 
                      float start_z, 
                      float end_x, 
                      float end_y, 
                      float end_z, 
                      int n_var_points, 
                      Configuration *c);

        void printValues();
        void iterate();
        void iterateWork();

        ~Variational3D();

        float rebalancePoints3D();

};