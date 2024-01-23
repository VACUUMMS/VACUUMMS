#include <iostream>
#include <stdio.h>
#include <math.h>
#include "variational.hh"
#include "quaternion.hh"

// These are provided and initialized in constants.o
//extern float epsilon;
//extern float sqrt_epsilon;


// prototype only
//float (*energy)(float x, float y);


void Variational3D::init(float _start_x, 
                    float _start_y, 
                    float _start_z, 
                    float _end_x, 
                    float _end_y, 
                    float _end_z, 
                    float _sigma, 
                    float _epsilon, 
                    int _n_var_points)
{
    // grab the passed values
    start_x = _start_x; 
    start_y = _start_y; 
    start_z = _start_z; 
    end_x = _end_x; 
    end_y = _end_y; 
    end_z = _end_z; 
    sigma = _sigma; 
    epsilon = _epsilon; 
    n_var_points = _n_var_points;

    // start constructing 

    var_x = new float[n_var_points];
    var_y = new float[n_var_points];
    var_z = new float[n_var_points];

    // initialize set of points
    
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = start_x + (i + 1) * (end_x - start_x) / (n_var_points + 1);
        var_y[i] = start_y + (i + 1) * (end_y - start_y) / (n_var_points + 1);
        var_z[i] = start_z + (i + 1) * (end_z - start_z) / (n_var_points + 1);
    }
}

/*
Variational3D::Variational3D(float _start_x, 
                             float _start_y, 
                             float _start_z, 
                             float _end_x, 
                             float _end_y, 
                             float _end_z, 
                             int _n_var_points, 
                             float(*_energy_function)(float x, float y, float z))
{
    init(_start_x, _start_y, _start_z, _end_x, _end_y, _end_z, _n_var_points);
    use_configuration_energy = false;
    configuration = nullptr;
    energy_function = _energy_function;
}
*/

Variational3D::Variational3D(float _start_x, 
                             float _start_y, 
                             float _start_z, 
                             float _end_x, 
                             float _end_y, 
                             float _end_z, 
                             float _sigma, 
                             float _epsilon, 
                             int _n_var_points, 
                             Configuration *_configuration)
{
    init(_start_x, _start_y, _start_z, _end_x, _end_y, _end_z, _sigma, _epsilon, _n_var_points);
    use_configuration_energy = true;
    configuration = _configuration;
    energy_function = nullptr;
}


Variational3D::~Variational3D()
{
    // free any resources allocated
    // std::cout << "destructing " << this << std::endl;
    delete var_x, var_y, var_z;
}

void Variational3D::printValues()
{
    std::cout << "#####" << start_x << " " << start_y << " " << start_z << std::endl;
    for (int i=0; i<n_var_points; i++)
        std::cout << "#####" << var_x[i] << " " << var_y[i] << " " << var_z[i] << std::endl;
    std::cout << "#####" << end_x << " " << end_y << " " << end_z << std::endl;
}
    
void Variational3D::setAlpha(float _alpha)
{
    alpha = _alpha;
}

void Variational3D::setAlphaMax(float _alpha_max)
{
    alpha_max = _alpha_max;
}

void Variational3D::iterate()
{
    // FTW: need to make a full pass before updating!!!
    float new_x[n_var_points], new_y[n_var_points], new_z[n_var_points];
    // looking along direction of line, point before and point after
    float fore_x, aft_x;
    float fore_y, aft_y;
    float fore_z, aft_z;

    //float alpha = 0.01; // nudge size?
    float beta = 0.1; // when step size is too large, multiply by beta

    for (int i=0; i<n_var_points; i++)
    {
        if (i>0)
        {
            fore_x = var_x[i-1];
            fore_y = var_y[i-1];
            fore_z = var_z[i-1];
        }
        else 
        {
            fore_x = start_x;
            fore_y = start_y;
            fore_z = start_z;
        }

        if (i+1<n_var_points) 
        {
            aft_x = var_x[i+1];
            aft_y = var_y[i+1];
            aft_z = var_z[i+1];
        }
        else
        {
            aft_x = end_x;
            aft_y = end_y;
            aft_z = end_z;
        }

        // tangent_* here is the vector tangent to the curve S at the point of interest
        float tangent_x = aft_x - fore_x;
        float tangent_y = aft_y - fore_y;
        float tangent_z = aft_z - fore_z;

//printf("got tangent vector (%f, %f, %f)\n", tangent_x, tangent_y, tangent_z);

        float u_x, u_y, u_z; // unit axis of rotation, to be extracted from tangent
        float theta = extract_axis(0, 0, 1, tangent_x, tangent_y, tangent_z, &u_x, &u_y, &u_z);
//printf("got theta = %f and axis u = (%f, %f, %f)\n ", theta, u_x, u_y, u_z);

        // variables to receive the values of the rotated i and j unit vectors
        float i_x, i_y, i_z; 
        rotate_vector(1, 0, 0, theta, u_x, u_y, u_z, &i_x, &i_y, &i_z);
//printf("i rotates to (%f, %f, %f)\n", i_x, i_y, i_z);
        float j_x, j_y, j_z; 
        rotate_vector(0, 1, 0, theta, u_x, u_y, u_z, &j_x, &j_y, &j_z);
//printf("j rotates to (%f, %f, %f)\n", j_x, j_y, j_z);
        // the directional vectors are ostensibly normalized

        // sample directional derivatives (directions perpendicular to curve)

        // get direction perpendicular, negative reciprocal of slope
        // float directional_y = -(aft_x - fore_x);
        // float directional_x = (aft_y - fore_y);

//        std::cout << "direction for var point " << i << ": (" << var_x[i] << ", " << var_y[i] << "):(" << directional_x << ", " << directional_y << ")" << "\tmagnitude: " << sqrt(directional_x * directional_x + directional_y * directional_y) << std::endl;

        // resize the direction vector to machine epsilon
        // directional_x *= sqrt_epsilon;
        // directional_y *= sqrt_epsilon;

        // resize the directionals as deltas

        i_x *= sqrt_epsilon;
        i_y *= sqrt_epsilon;
        i_z *= sqrt_epsilon;
//printf("i resized to (%f, %f, %f)\n", i_x, i_y, i_z);

        j_x *= sqrt_epsilon;
        j_y *= sqrt_epsilon;
        j_z *= sqrt_epsilon;
//printf("j resized to (%f, %f, %f)\n", j_x, j_y, j_z);

//printf("using resized directional x,y: %f, %f\n", directional_x, directional_y);
//printf("sqrt_epsilon = %0.012f\n", sqrt_epsilon);

        // energy at var[i], and at directional points
        float energy_center, energy_i, energy_j;
        
/*
        // sample energy to evaluate derivative
        float sample_left_x = var_x[i] - directional_x;
        float sample_left_y = var_y[i] - directional_y;
        float sample_left_z = var_z[i] - directional_z;
        float sample_right_x = var_x[i] + directional_x;
        float sample_right_y = var_y[i] + directional_y;
        float sample_right_z = var_z[i] + directional_z;
*/
        //float energy_left;
        //float energy_right;

        if (use_configuration_energy)
        {
            // energy_left = configuration->insertionEnergy2D(sample_left_x, sample_left_y);
            // energy_right = configuration->insertionEnergy2D(sample_right_x, sample_right_y);
            energy_center = configuration->insertionEnergy(var_x[i], var_y[i], var_z[i], sigma, epsilon);
            energy_i = configuration->insertionEnergy(var_x[i] + i_x, var_y[i] + i_y, var_z[i] + i_z, sigma, epsilon);
            energy_j = configuration->insertionEnergy(var_x[i] + j_x, var_y[i] + j_y, var_z[i] + j_z, sigma, epsilon);
        }
        else 
        {
            printf("using non-configuration energy function not implemented.\n");
            exit(1);
        }

//printf("energy left/right: %.012f <--> %.012f\n", energy_left, energy_right);
        // Not using alpha step size, just nudging. may need to normalize step size somehow
        // dE = (dE/dx)dx + (dE/dy)dy
//        float dE = energy_right - energy_left;
//std::cout << "dE = " << dE << "\n";

        float dE_i = energy_i - energy_center;
        float dE_j = energy_j - energy_center;

        // update position
//        new_x[i] = var_x[i] - alpha * dE * directional_x;
//        new_y[i] = var_y[i] - alpha * dE * directional_y;

        float delta_x = - alpha * (dE_i * i_x + dE_j * j_x);
        float delta_y = - alpha * (dE_i * i_y + dE_j * j_y);
        float delta_z = - alpha * (dE_i * i_z + dE_j * j_z);
        float delta_sq = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
        float delta = sqrt(delta_sq);
printf("Got dE_i = %f, dE_j = %f, delta = (%f, %f, %f) |delta| = %f\n", dE_i, dE_j, delta_x, delta_y, delta_z, delta);

// try rescaling if delta is too big?
//        if (delta_sq > alpha) 
//        {
//            delta_x *= beta;
//            delta_y *= beta;
//            delta_z *= beta;
//        }

        new_x[i] = var_x[i] + delta_x; 
        new_y[i] = var_y[i] + delta_y;
        new_z[i] = var_z[i] + delta_z;
    }

    // FTW: After full pass, copy the updates back.
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = new_x[i];
        var_y[i] = new_y[i];
        var_z[i] = new_z[i];
    }
}

//########################################## rebalance ###########################################//
float Variational3D::rebalancePoints3D()
{
    float original_curve_length = calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, var_x, var_y, var_z);

    float new_forward_var_x[n_var_points];
    float new_forward_var_y[n_var_points];
    float new_forward_var_z[n_var_points];

    int status;
    status = respaceKernel(/*n_var_points,*/ start_x, start_y, start_z, end_x, end_y, end_z, var_x, var_y, var_z, new_forward_var_x, new_forward_var_y, new_forward_var_z);
    if (status != 0) exit(status);

    float forward_curve_length = calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, new_forward_var_x, new_forward_var_y, new_forward_var_z);

    // now run it backward

    float backward_var_x_in[n_var_points];
    float backward_var_y_in[n_var_points];
    float backward_var_z_in[n_var_points];
    
    float backward_var_x_out[n_var_points];
    float backward_var_y_out[n_var_points];
    float backward_var_z_out[n_var_points];
    
/*
    //for (int i=0; i<n_var_points; i++) printf("before: backward_var_in = (%f, %f, %f)\n", backward_var_x_in[i], backward_var_y_in[i], backward_var_z_in[i]);
printf("--------------------------\n");
    for (int i=0; i<n_var_points; i++) printf("before: var[%d] = (%f, %f, %f)\n", i, var_x[i], var_y[i], var_z[i]);
printf("--------------------------\n");
*/
    // flip the curve
    for (int i=0; i<n_var_points; i++)
    {
        backward_var_x_in[i] = var_x[n_var_points - i - 1];
        backward_var_y_in[i] = var_y[n_var_points - i - 1];
        backward_var_z_in[i] = var_z[n_var_points - i - 1];
    }

/*
printf("--------------------------\n");
    for (int i=0; i<n_var_points; i++) printf("backward_var_in = (%f, %f, %f)\n", backward_var_x_in[i], backward_var_y_in[i], backward_var_z_in[i]);
printf("--------------------------\n");
*/

    status = respaceKernel(/*n_var_points,*/ end_x, end_y, end_z, start_x, start_y, start_z, backward_var_x_in, backward_var_y_in, backward_var_z_in, backward_var_x_out, backward_var_y_out, backward_var_z_out);
    if (status != 0) exit(status);

    // and flip it back
    float new_backward_var_x[n_var_points];
    float new_backward_var_y[n_var_points];
    float new_backward_var_z[n_var_points];
    
    for (int i=0; i<n_var_points; i++)
    {
        new_backward_var_x[i] = backward_var_x_out[n_var_points - i - 1];
        new_backward_var_y[i] = backward_var_y_out[n_var_points - i - 1];
        new_backward_var_z[i] = backward_var_z_out[n_var_points - i - 1];
    }

    float backward_curve_length = calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, new_forward_var_x, new_forward_var_y, new_forward_var_z);

// should match
printf("original backward curve_length = %f\n", backward_curve_length);

    // combine the results
    float new_var_x[n_var_points];
    float new_var_y[n_var_points];
    float new_var_z[n_var_points];

    for (int i=0; i<n_var_points; i++) 
    {
/*
        new_var_x[i] = new_forward_var_x[i];
        new_var_y[i] = new_forward_var_y[i];
        new_var_z[i] = new_forward_var_z[i];

        new_var_x[i] = new_backward_var_x[i];
        new_var_y[i] = new_backward_var_y[i];
        new_var_z[i] = new_backward_var_z[i];
*/
        new_var_x[i] = 0.5 * (new_forward_var_x[i] + new_backward_var_x[i]);
        new_var_y[i] = 0.5 * (new_forward_var_y[i] + new_backward_var_y[i]);
        new_var_z[i] = 0.5 * (new_forward_var_z[i] + new_backward_var_z[i]);
    }

    // check new combined curve length, add distance to each variational point to get shrinkage
    float new_curve_length = calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, new_var_x, new_var_y, new_var_z);

/*
    float prev_x = start_x;
    float prev_y = start_y;
    float prev_z = start_z;
    
    for (int i=0; i < n_var_points; i++)
    {
        new_curve_length += sqrt(
                            (new_var_x[i] - prev_x) * (new_var_x[i] - prev_x) 
                          + (new_var_y[i] - prev_y) * (new_var_y[i] - prev_y)
                          + (new_var_z[i] - prev_z) * (new_var_z[i] - prev_z)
                            );

        // update for next iteration
        prev_x = new_var_x[i]; 
        prev_y = new_var_y[i];
        prev_z = new_var_z[i];
    }

    // add the last piece
    new_curve_length += sqrt(
                        (end_x - prev_x) * (end_x - prev_x) 
                      + (end_y - prev_y) * (end_y - prev_y)
                      + (end_z - prev_z) * (end_z - prev_z)
                        );
 */   
    float shrinkage = new_curve_length / original_curve_length;

    // copy out
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = new_var_x[i];
        var_y[i] = new_var_y[i];
        var_z[i] = new_var_z[i];
    }
   
    return shrinkage;
}

//#################################################### factoring to Kernel #################################//
// return value is 0 if successful, 1 if not successful
int Variational3D::respaceKernel(float _start_x, float _start_y, float _start_z, 
                                 float _end_x, float _end_y, float _end_z,
                                 float _var_x[], float _var_y[], float _var_z[], 
                                 float _new_var_x[], float _new_var_y[], float _new_var_z[])
{
    float curve_length = calculateCurveLength(_start_x, _start_y, _start_z, _end_x, _end_y, _end_z, _var_x, _var_y, _var_z);
    // Cut each new segment along old path to this length
    float new_segment_length = curve_length / (/*_*/n_var_points + 1.0);

    // now the respace
    float cursor_x = _start_x;
    float cursor_y = _start_y;
    float cursor_z = _start_z;
    float remainder_x = 0;         
    float remainder_y = 0;         
    float remainder_z = 0;         
    float remainder = 0;
    int old_point = 0;

    int all_segments_added = 0;
    for (int new_point=0; new_point</*_*/n_var_points; ) 
    {
int something_happened = 0;
        // grab old segments until there's enough curve to cut at least one new segment
        while ((remainder < new_segment_length) && (old_point <= /*_*/n_var_points))
        {
            float segment_x, segment_y, segment_z;
            if (old_point < /*_*/n_var_points) 
            {
                segment_x = _var_x[old_point] - cursor_x;
                segment_y = _var_y[old_point] - cursor_y;
                segment_z = _var_z[old_point] - cursor_z;
            }
            else if (old_point == /*_*/n_var_points)
            {
                segment_x = _end_x - cursor_x;
                segment_y = _end_y - cursor_y;
                segment_z = _end_z - cursor_z;
                all_segments_added = 1;
            }
            remainder_x += segment_x;
            remainder_y += segment_y;
            remainder_z += segment_z;
            remainder = sqrt(remainder_x*remainder_x +
                             remainder_y*remainder_y +
                             remainder_z*remainder_z);
            cursor_x = _var_x[old_point];
            cursor_y = _var_y[old_point];
            cursor_z = _var_z[old_point];
            old_point++;
            something_happened++;
        }
        
        // now cut new segments until there's not enough left
        //while ((new_segment_length <= remainder) && (i <= n_var_points)) // this should run stop adding points, even with remaining segment. 
        while ((new_segment_length <= remainder) && (new_point < n_var_points)) // this should run stop adding points, even with remaining segment. 
        {
            // get projections of new_segment_length onto x, y directions
            float delta_x = (remainder_x / remainder) * new_segment_length;
            float delta_y = (remainder_y / remainder) * new_segment_length;
            float delta_z = (remainder_z / remainder) * new_segment_length;
            if (new_point==0) 
            {
                _new_var_x[new_point] = _start_x + delta_x;
                _new_var_y[new_point] = _start_y + delta_y;
                _new_var_z[new_point] = _start_z + delta_z;
            }
            else if (new_point < n_var_points)
            {
                _new_var_x[new_point] = _new_var_x[new_point-1] + delta_x;
                _new_var_y[new_point] = _new_var_y[new_point-1] + delta_y;
                _new_var_z[new_point] = _new_var_z[new_point-1] + delta_z;
            }

            // subtract the last piece from remainder since it's been cut
            remainder_x -= delta_x;
            remainder_y -= delta_y;
            remainder_z -= delta_z;
            remainder = sqrt(remainder_x*remainder_x 
                        + remainder_y*remainder_y
                        + remainder_z*remainder_z);

            new_point++; // move to next new point
            something_happened++;
        }            

        if (!something_happened)
        {
            printf("Curve cannot be set. Wedged, exiting.\n");
            printf("new_point = %d, old_point = %d, all_segments_added = %d\n", new_point, old_point, all_segments_added);
            printf("remainder = %f, new_segment_length = %f\n", remainder, new_segment_length);
            printf("consider checking / reducing the value of alpha = %f\n", alpha);
            return 1;
        }
    }

    return 0;
} // respaceKernel()


float Variational3D::calculateCurveLength(float _start_x, float _start_y, float _start_z, float _end_x, float _end_y, float _end_z, float _var_x[], float _var_y[], float _var_z[])
{
    // get the total length of curve
    float prev_x = _start_x;
    float prev_y = _start_y;
    float prev_z = _start_z;
    float curve_length = 0;

    // add distance to each variational point
    for (int i=0; i < n_var_points; i++)
    {
        curve_length += sqrt(
                        (_var_x[i] - prev_x) * (_var_x[i] - prev_x) 
                      + (_var_y[i] - prev_y) * (_var_y[i] - prev_y)
                      + (_var_z[i] - prev_z) * (_var_z[i] - prev_z)
                        );

        // update for next iteration
        prev_x = _var_x[i]; 
        prev_y = _var_y[i];
        prev_z = _var_z[i];
    }

    // add the last piece
    curve_length += sqrt(
                    (_end_x - prev_x) * (_end_x - prev_x) 
                  + (_end_y - prev_y) * (_end_y - prev_y)
                  + (_end_z - prev_z) * (_end_z - prev_z)
                    );

    return curve_length;
}


float Variational3D::adaptiveIterateAndUpdate()
{
    float curve_length = Variational3D::calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, var_x, var_y, var_z);
    float new_x[n_var_points], new_y[n_var_points], new_z[n_var_points];

attempt_iteration:

fprintf(stdout, "Attempting iteration with alpha = %f\n", alpha);
    
    float fore_x, aft_x;
    float fore_y, aft_y;
    float fore_z, aft_z;

    float beta = 1.1; // adaptive parameter, increase/decrease alpha by this factor

    for (int i=0; i<n_var_points; i++)
    {
        if (i>0)
        {
            fore_x = var_x[i-1];
            fore_y = var_y[i-1];
            fore_z = var_z[i-1];
        }
        else 
        {
            fore_x = start_x;
            fore_y = start_y;
            fore_z = start_z;
        }

        if (i+1<n_var_points) 
        {
            aft_x = var_x[i+1];
            aft_y = var_y[i+1];
            aft_z = var_z[i+1];
        }
        else
        {
            aft_x = end_x;
            aft_y = end_y;
            aft_z = end_z;
        }

        // tangent_* here is the vector tangent to the curve S at the point of interest
        float tangent_x = aft_x - fore_x;
        float tangent_y = aft_y - fore_y;
        float tangent_z = aft_z - fore_z;

        float u_x, u_y, u_z; // unit axis of rotation, to be extracted from tangent
        float theta = extract_axis(0, 0, 1, tangent_x, tangent_y, tangent_z, &u_x, &u_y, &u_z);

        // variables to receive the values of the rotated i and j unit vectors
        float i_x, i_y, i_z; 
        rotate_vector(1, 0, 0, theta, u_x, u_y, u_z, &i_x, &i_y, &i_z);
        float j_x, j_y, j_z; 
        rotate_vector(0, 1, 0, theta, u_x, u_y, u_z, &j_x, &j_y, &j_z);

        // resize the directionals as deltas

        i_x *= sqrt_epsilon;
        i_y *= sqrt_epsilon;
        i_z *= sqrt_epsilon;

        j_x *= sqrt_epsilon;
        j_y *= sqrt_epsilon;
        j_z *= sqrt_epsilon;

        // energy at var[i], and at directional points
        float energy_center, energy_i, energy_j;
        
        if (use_configuration_energy)
        {
            energy_center = configuration->insertionEnergy(var_x[i], var_y[i], var_z[i], sigma, epsilon);
            energy_i = configuration->insertionEnergy(var_x[i] + i_x, var_y[i] + i_y, var_z[i] + i_z, sigma, epsilon);
            energy_j = configuration->insertionEnergy(var_x[i] + j_x, var_y[i] + j_y, var_z[i] + j_z, sigma, epsilon);
        }
        else 
        {
            printf("using non-configuration energy function not implemented.\n");
            exit(1);
        }

        float dE_i = energy_i - energy_center;
        float dE_j = energy_j - energy_center;

        float delta_x = - alpha * (dE_i * i_x + dE_j * j_x);
        float delta_y = - alpha * (dE_i * i_y + dE_j * j_y);
        float delta_z = - alpha * (dE_i * i_z + dE_j * j_z);
        float delta_sq = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
        float delta = sqrt(delta_sq);

fprintf(stdout, "Got dE_i = %f, dE_j = %f, delta = (%f, %f, %f) |delta| = %f\n", dE_i, dE_j, delta_x, delta_y, delta_z, delta); //FTW
float delta_max = 0.01;

        if (delta > delta_max)
        {
            printf("delta = %f > delta_max = %f, rescaling.\n", delta, delta_max);
            delta_x *= (delta_max/delta);
            delta_y *= (delta_max/delta);
            delta_z *= (delta_max/delta);
        }

        new_x[i] = var_x[i] + delta_x; 
        new_y[i] = var_y[i] + delta_y;
        new_z[i] = var_z[i] + delta_z;
    }

/* hold off on update until the rebalance is OK
    // FTW: After full pass, copy the updates back.
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = new_x[i];
        var_y[i] = new_y[i];
        var_z[i] = new_z[i];
    }
*/

// end of iterate part, have new_x, new_y, and new_z to carry forward

//########################################## rebalance ###########################################//

    // run it forward

    float new_forward_var_x[n_var_points];
    float new_forward_var_y[n_var_points];
    float new_forward_var_z[n_var_points];

    int forward_status = respaceKernel(start_x, start_y, start_z, end_x, end_y, end_z, new_x, new_y, new_z, new_forward_var_x, new_forward_var_y, new_forward_var_z);

    if (forward_status != 0) 
    {
//        alpha /= beta;
        fprintf(stdout, "forward respace failed, shrinking alpha = %f to %f and re-running.\n", alpha, alpha/=beta);
        goto attempt_iteration;
    }

    // now run it backward

    float backward_var_x_in[n_var_points];
    float backward_var_y_in[n_var_points];
    float backward_var_z_in[n_var_points];
    
    float backward_var_x_out[n_var_points];
    float backward_var_y_out[n_var_points];
    float backward_var_z_out[n_var_points];
    
    // flip the curve
    for (int i=0; i<n_var_points; i++)
    {
//        backward_var_x_in[i] = var_x[n_var_points - i - 1];
//        backward_var_y_in[i] = var_y[n_var_points - i - 1];
//        backward_var_z_in[i] = var_z[n_var_points - i - 1];
        backward_var_x_in[i] = new_x[n_var_points - i - 1];
        backward_var_y_in[i] = new_y[n_var_points - i - 1];
        backward_var_z_in[i] = new_z[n_var_points - i - 1];
    }

    int backward_status = respaceKernel(end_x, end_y, end_z, start_x, start_y, start_z, backward_var_x_in, backward_var_y_in, backward_var_z_in, backward_var_x_out, backward_var_y_out, backward_var_z_out);

    if (backward_status != 0) 
    {
        fprintf(stdout, "backward respace failed, shrinking alpha = %f to %f and re-running.\n", alpha, alpha/=beta);
        goto attempt_iteration;
    }

    // and flip it back
    float new_backward_var_x[n_var_points];
    float new_backward_var_y[n_var_points];
    float new_backward_var_z[n_var_points];
    
    for (int i=0; i<n_var_points; i++)
    {
        new_backward_var_x[i] = backward_var_x_out[n_var_points - i - 1];
        new_backward_var_y[i] = backward_var_y_out[n_var_points - i - 1];
        new_backward_var_z[i] = backward_var_z_out[n_var_points - i - 1];
    }

//assert(forward_curve_length ~ backward_curve_length) ?? FTW

    // combine the forward and backward results
    float respace_var_x[n_var_points];
    float respace_var_y[n_var_points];
    float respace_var_z[n_var_points];

    for (int i=0; i<n_var_points; i++) 
    {
        respace_var_x[i] = 0.5 * (new_forward_var_x[i] + new_backward_var_x[i]);
        respace_var_y[i] = 0.5 * (new_forward_var_y[i] + new_backward_var_y[i]);
        respace_var_z[i] = 0.5 * (new_forward_var_z[i] + new_backward_var_z[i]);
    }

    float new_curve_length = calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, respace_var_x, respace_var_y, respace_var_z);

    float shrinkage = new_curve_length / curve_length;

// FTW make sure this part is ready
    // copy out
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = respace_var_x[i];
        var_y[i] = respace_var_y[i];
        var_z[i] = respace_var_z[i];
    }

//float alpha_max = 0.1;
    if (alpha < alpha_max) fprintf(stdout, "successful update and rebalance. increasing alpha: %f * %f = %f\n", alpha, beta, alpha *= beta); // success, so increase alpha
   
    return shrinkage;

} // adaptiveIterateAndUpdate()

