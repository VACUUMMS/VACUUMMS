// 1D, 2D, and 3D cases. Each returns the shrinkage, i.e. the ratio of the length of the new curve to the original
float rebalance_points_1D(float _start_x, float _end_x, float* _var_x, int _n_var_points);
float rebalance_points_2D(float _start_x, float _start_y, float _end_x, float _end_y, float* _var_x, float* _var_y, int _n_var_points);
float rebalance_points_3D(float _start_x, float _start_y, float _start_z, float _end_x, float _end_y, float _end_z, float* _var_x, float* _var_y, float* _var_z, int _n_var_points);

