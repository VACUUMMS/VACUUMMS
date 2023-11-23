// variational.hh
#include "configuration.hh"


class Variational2D
{
    private:

        float start_x; 
        float start_y; 
        float end_x; 
        float end_y; 
        int n_var_points;
        Configuration configuration;
        float(*energy_function)(float x, float y);

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
                      Configuration c);

        void printValues();
        void iterate();

        ~Variational2D();
};

// 1D, 2D, and 3D cases. Each returns the shrinkage, i.e. the ratio of the length of the new curve to the original
float rebalance_points_1D(float _start_x, float _end_x, float* _var_x, int _n_var_points);
float rebalance_points_2D(float _start_x, float _start_y, float _end_x, float _end_y, float* _var_x, float* _var_y, int _n_var_points);
float rebalance_points_3D(float _start_x, float _start_y, float _start_z, float _end_x, float _end_y, float _end_z, float* _var_x, float* _var_y, float* _var_z, int _n_var_points);

