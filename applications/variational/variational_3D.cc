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
                    int _n_var_points)
{
    // grab the passed values
    start_x = _start_x; 
    start_y = _start_y; 
    start_z = _start_z; 
    end_x = _end_x; 
    end_y = _end_y; 
    end_z = _end_z; 
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

Variational3D::Variational3D(float _start_x, 
                             float _start_y, 
                             float _start_z, 
                             float _end_x, 
                             float _end_y, 
                             float _end_z, 
                             int _n_var_points, 
                             Configuration *_configuration)
{
    init(_start_x, _start_y, _start_z, _end_x, _end_y, _end_z, _n_var_points);
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
    

void Variational3D::iterate()
{
    // FTW: need to make a full pass before updating!!!
    float new_x[n_var_points], new_y[n_var_points], new_z[n_var_points];
    // looking along direction of line, point before and point after
    float fore_x, aft_x;
    float fore_y, aft_y;
    float fore_z, aft_z;

    float alpha = 0.1; // nudge size?

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
            energy_center = configuration->insertionEnergy(var_x[i], var_y[i], var_z[i]);
            energy_i = configuration->insertionEnergy(var_x[i] + i_x, var_y[i] + i_y, var_z[i] + i_z);
            energy_j = configuration->insertionEnergy(var_x[i] + j_x, var_y[i] + j_y, var_z[i] + j_z);
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
printf("Got delta = (%f, %f, %f)\n", delta_x, delta_y, delta_z);

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

float Variational3D::rebalancePoints3D()
{

//printf("FTW: entering rebalancePoints3D()\n");

    // get the total length of curve
    float prev_x = start_x;
    float prev_y = start_y;
    float prev_z = start_z;
    float curve_length = 0;

    // add distance to each variational point
    for (int i=0; i < n_var_points; i++)
    {
        curve_length += sqrt(
                        (var_x[i] - prev_x) * (var_x[i] - prev_x) 
                      + (var_y[i] - prev_y) * (var_y[i] - prev_y)
                      + (var_z[i] - prev_z) * (var_z[i] - prev_z)
                        );

        // update for next iteration
        prev_x = var_x[i]; 
        prev_y = var_y[i];
        prev_z = var_z[i];
    }

    // add the last piece
    curve_length += sqrt(
                    (end_x - prev_x) * (end_x - prev_x) 
                  + (end_y - prev_y) * (end_y - prev_y)
                  + (end_z - prev_z) * (end_z - prev_z)
                    );
    // Cut each new segment along old path to this length
    float new_segment_length = curve_length / (n_var_points + 1.0);

//printf("FTW: starting respace in rebalancePoints3D()\n");
    // now the respace
    float cursor_x = start_x;
    float cursor_y = start_y;
    float cursor_z = start_z;
    float remainder_x = 0;         
    float remainder_y = 0;         
    float remainder_z = 0;         
    float remainder = 0;
    int old_point = 0;

    float new_var_x[n_var_points];
    float new_var_y[n_var_points];
    float new_var_z[n_var_points];

// something in this loop is corrupting the stack, probably overstepping array bounds (i or old_point?)
    for (int i=0; i<n_var_points; ) 
    {
        // grab old segments until there's enough to cut at least one new one
        // and while not EOL (end of line)
        // while (remainder_x < new_segment_length)
        while ((remainder < new_segment_length) && (old_point <= n_var_points))
        {
            float segment_x, segment_y, segment_z;
            if (old_point < n_var_points) 
            {
                segment_x = var_x[old_point] - cursor_x;
                segment_y = var_y[old_point] - cursor_y;
                segment_z = var_z[old_point] - cursor_z;
            }
            else if (old_point == n_var_points)
            {
                segment_x = end_x - cursor_x;
                segment_y = end_y - cursor_y;
                segment_z = end_z - cursor_z;
            }
            else
            {
                // This should never be reached.
                printf("rebalancePoints(): no more old segments to add, but not enough to cut new one.\n");
                exit(1);
            }
            remainder_x += segment_x;
            remainder_y += segment_y;
            remainder_z += segment_z;
            remainder = sqrt(remainder_x*remainder_x 
                        + remainder_y*remainder_y 
                        + remainder_z*remainder_z);
// this is accessing past array bounds? moving cursor to bad point, but not using it?
// but only for old_point == n_var_points, which is handled by delivering the end_x, end_y, end_z, so should be OK
//printf("adding old point %d to grow segment by (%f, %f, %f)\n", old_point, segment_x, segment_y, segment_z);
//printf("remainder is now = %f\n", remainder);
            cursor_x = var_x[old_point];
            cursor_y = var_y[old_point];
            cursor_z = var_z[old_point];
            old_point++;
        }
        
        // now cut new segments until there's not enough left
        // need to cover case of equality of segment lengths, so adding it here or above
        // while (new_segment_length <= remainder)
        while ((new_segment_length <= remainder) && (i <= n_var_points)) // this should run stop adding points, even with remaining segment. 
        {
            // get projections of new_segment_length onto x, y directions
            float delta_x = (remainder_x / remainder) * new_segment_length;
            float delta_y = (remainder_y / remainder) * new_segment_length;
            float delta_z = (remainder_z / remainder) * new_segment_length;
            if (i==0) 
            {
                new_var_x[i] = start_x + delta_x;
                new_var_y[i] = start_y + delta_y;
                new_var_z[i] = start_z + delta_z;
            }
            else if (i < n_var_points)
            {
                new_var_x[i] = new_var_x[i-1] + delta_x;
                new_var_y[i] = new_var_y[i-1] + delta_y;
                new_var_z[i] = new_var_z[i-1] + delta_z;
            }
            else if (i == n_var_points)
            {
                // last piece ends at end_x, end_y, end_z
            }
            remainder_x -= delta_x;
            remainder_y -= delta_y;
            remainder_z -= delta_z;
            remainder = sqrt(remainder_x*remainder_x 
                        + remainder_y*remainder_y
                        + remainder_z*remainder_z);
//printf("cutting new segment to generate point %d at (%f, %f, %f)\n", i, new_var_x[i], new_var_y[i], new_var_z[i]);
//printf("after cutting, remainder is now = %f\n", remainder);
            i++; // move to next new point

//if (i >= n_var_points) printf("index i=%d going out of range\n", i);
//if (old_point >= n_var_points) printf("index old_point=%d going out of range\n", old_point);

        }            
    }

//printf("after adding and cutting all points, segment remainder = %f\n", remainder);
 
//FTW pseudo
//printf("FTW checking new curve\n");

    // check new curve length, add distance to each variational point
    float new_curve_length = 0;
    prev_x = start_x;
    prev_y = start_y;
    prev_z = start_z;
    
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
//printf("FTW adding last piece \n");

    // add the last piece
    new_curve_length += sqrt(
                        (end_x - prev_x) * (end_x - prev_x) 
                      + (end_y - prev_y) * (end_y - prev_y)
                      + (end_z - prev_z) * (end_z - prev_z)
                        );
    
    float shrinkage = new_curve_length / curve_length;

//printf("FTW preparing to copy out \n");

    // copy out
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = new_var_x[i];
        var_y[i] = new_var_y[i];
        var_z[i] = new_var_z[i];
    }
   
//printf("FTW returning \n");
    return shrinkage;
}
