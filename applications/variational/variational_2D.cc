#include <iostream>
#include <stdio.h>
#include <math.h>
#include "variational.hh"

// These are provided and initialized in constants.o
//extern float epsilon;
//extern float sqrt_epsilon;


// prototype only
//float (*energy)(float x, float y);


void Variational2D::init(float _start_x, 
                    float _start_y, 
                    float _end_x, 
                    float _end_y, 
                    int _n_var_points)
{
    // grab the passed values
    start_x = _start_x; 
    start_y = _start_y; 
    end_x = _end_x; 
    end_y = _end_y; 
    n_var_points = _n_var_points;

    // std::cout << "constructing/initing Variational2D: " << this << std::endl;

    // start constructing 

    var_x = new float[n_var_points];
    var_y = new float[n_var_points];

    // initialize set of points
    
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = start_x + (i + 1) * (end_x - start_x) / (n_var_points + 1);
        var_y[i] = start_y + (i + 1) * (end_y - start_y) / (n_var_points + 1);
        // printf("var[i] = {%f,%f}\n", var_x[i], var_y[i]);
    }
}

Variational2D::Variational2D(float _start_x, 
                             float _start_y, 
                             float _end_x, 
                             float _end_y, 
                             int _n_var_points, 
                             float(*_energy_function)(float x, float y))
{
    init(_start_x, _start_y, _end_x, _end_y, _n_var_points);
    use_configuration_energy = false;
    configuration = nullptr;
    energy_function = _energy_function;
}

Variational2D::Variational2D(float _start_x, 
                             float _start_y, 
                             float _end_x, 
                             float _end_y, 
                             int _n_var_points, 
                             Configuration *_configuration)
{
    init(_start_x, _start_y, _end_x, _end_y, _n_var_points);
    use_configuration_energy = true;
    configuration = _configuration;
    energy_function = nullptr;
}


Variational2D::~Variational2D()
{
    // free any resources allocated
    // std::cout << "destructing " << this << std::endl;
    delete var_x, var_y;
}

void Variational2D::printValues()
{
    std::cout << "#####" << start_x << " " << start_y << std::endl;
    for (int i=0; i<n_var_points; i++)
        std::cout << "#####" << var_x[i] << " " << var_y[i] << std::endl;
    std::cout << "#####" << end_x << " " << end_y << std::endl;
}
    

void Variational2D::iterate()
{
    // FTW: need to make a full pass before updating!!!
    float new_x[n_var_points], new_y[n_var_points];
    // looking along direction of line, point before and point after
    float fore_x, aft_x, fore_y, aft_y;

    float alpha = 100.0; // nudge size?

    for (int i=0; i<n_var_points; i++)
    {
        if (i>0) fore_x = var_x[i-1];
        else fore_x = start_x;
        if (i+1<n_var_points) aft_x = var_x[i+1];
        else aft_x = end_x;
        if (i>0) fore_y = var_y[i-1];
        else fore_y = start_y;
        if (i+1<n_var_points) aft_y = var_y[i+1];
        else aft_y = end_y;

        // sample directional derivative (direction perpendicular to curve)

        // get direction perpendicular, negative reciprocal of slope
        float directional_y = -(aft_x - fore_x);
        float directional_x = (aft_y - fore_y);

        // normalize the direction vector
        float directional_magnitude = sqrt(directional_x*directional_x + directional_y*directional_y);
        directional_x /= directional_magnitude;
        directional_y /= directional_magnitude;

//FTW
        std::cout << "direction for var point " << i << ": (" << var_x[i] << ", " << var_y[i] << "):(" << directional_x << ", " << directional_y << ")" << "\tmagnitude: " << sqrt(directional_x * directional_x + directional_y * directional_y) << std::endl;

        // resize the direction vector to machine epsilon
        directional_x *= sqrt_epsilon;
        directional_y *= sqrt_epsilon;

printf("using resized directional x,y: %f, %f\n", directional_x, directional_y);
printf("sqrt_epsilon = %0.012f\n", sqrt_epsilon);

        // sample energy to evaluate derivative
        float sample_left_x = var_x[i] - directional_x;
        float sample_left_y = var_y[i] - directional_y;
        float sample_right_x = var_x[i] + directional_x;
        float sample_right_y = var_y[i] + directional_y;

        float energy_left;
        float energy_right;

        if (use_configuration_energy)
        {
            //energy_left = configuration->insertionEnergy2D(sample_left_x, sample_left_y);
            //energy_right = configuration->insertionEnergy2D(sample_right_x, sample_right_y);
            energy_left = configuration->insertionEnergy(sample_left_x, sample_left_y, 0);
            energy_right = configuration->insertionEnergy(sample_right_x, sample_right_y, 0);
        }
        else 
        {
            energy_left = energy_function(sample_left_x, sample_left_y);
            energy_right = energy_function(sample_right_x, sample_right_y);
        }

printf("energy left/right: %.012f <--> %.012f\n", energy_left, energy_right);
        // Not using alpha step size, just nudging. may need to normalize step size somehow
        // dE = (dE/dx)dx + (dE/dy)dy
        float dE = energy_right - energy_left;
//std::cout << "dE = " << dE << "\n";

        // update position
        new_x[i] = var_x[i] - alpha * dE * directional_x;
        new_y[i] = var_y[i] - alpha * dE * directional_y;
    }

    // FTW: After full pass, copy the updates back.
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = new_x[i];
        var_y[i] = new_y[i];
    }

    /* FTW
    printf("%f %f\n", start_x, start_y);
    for (int i=0; i<n_var_points; i++)
    {
        printf("%f %f\n", var_x[i], var_y[i]);
    }
    printf("%f %f\n", end_x, end_y);
    */
}


/* Minimize the *work* over the path 
void Variational2D::iterateWork()
{
    // FTW: need to make a full pass before updating!!!
    float new_x[n_var_points], new_y[n_var_points];
    // looking along direction of line, point before and point after
    float fore_x, aft_x, fore_y, aft_y;

    float alpha = 100.0; // nudge size?

    for (int i=0; i<n_var_points; i++)
    {
        if (i>0) fore_x = var_x[i-1];
        else fore_x = start_x;
        if (i+1<n_var_points) aft_x = var_x[i+1];
        else aft_x = end_x;
        if (i>0) fore_y = var_y[i-1];
        else fore_y = start_y;
        if (i+1<n_var_points) aft_y = var_y[i+1];
        else aft_y = end_y;

        // sample directional derivative (direction perpendicular to curve)

        // get direction perpendicular, negative reciprocal of slope
        float directional_y = -(aft_x - fore_x);
        float directional_x = (aft_y - fore_y);

        // normalize the direction vector
        float directional_magnitude = sqrt(directional_x*directional_x + directional_y*directional_y);
        directional_x /= directional_magnitude;
        directional_y /= directional_magnitude;

//FTW
        std::cout << "direction for var point " << i << ": (" << var_x[i] << ", " << var_y[i] << "):(" << directional_x << ", " << directional_y << ")" << "\tmagnitude: " << sqrt(directional_x * directional_x + directional_y * directional_y) << std::endl;

        // resize the direction vector to machine epsilon
        directional_x *= sqrt_epsilon;
        directional_y *= sqrt_epsilon;

printf("using resized directional x,y: %f, %f\n", directional_x, directional_y);
printf("sqrt_epsilon = %0.012f\n", sqrt_epsilon);

        // sample energy to evaluate derivative
        float sample_left_x = var_x[i] - directional_x;
        float sample_left_y = var_y[i] - directional_y;
        float sample_right_x = var_x[i] + directional_x;
        float sample_right_y = var_y[i] + directional_y;

        float energy_left;
        float energy_right;

//FTW this part is same as original case
        if (use_configuration_energy)
        {
            energy_left = configuration->insertionEnergy2D(sample_left_x, sample_left_y);
            energy_right = configuration->insertionEnergy2D(sample_right_x, sample_right_y);
        }
        else 
        {
            energy_left = energy_function(sample_left_x, sample_left_y);
            energy_right = energy_function(sample_right_x, sample_right_y);
        }

// Below is where we need to account for the change in work

// calculate the gradient of the directional derivative of energy, instead of the gradient of energy

printf("energy left/right: %.012f <--> %.012f\n", energy_left, energy_right);
        // Not using alpha step size, just nudging. may need to normalize step size somehow
        // dE = (dE/dx)dx + (dE/dy)dy
//        float dE = energy_right - energy_left;


//FTW       work_left = function of sample_left_xy, delta_xy, etc. for var_xy[n-1], [n+1]; 
        // get energy at fore, mid, and aft points to calculate work

        float fore_x, fore_y;
        if (i == 0) // handle special case of first point
        {
            fore_x = start_x;
            fore_y = start_y;
        }
        else 
        {
            fore_x = var_x[i-1];
            fore_y = var_y[i-1];
        }

        float aft_x, aft_y;
        if (i == (n_var_points-1)) // handle special case of last point
        {
            aft_x = end_x;
            aft_y = end_y;
        }
        else 
        {
            aft_x = var_x[i+1];
            aft_y = var_y[i+1];
        }
        
        float energy_mid_left = energy_function(sample_left_x, sample_left_y);
        float energy_mid_right = energy_function(sample_right_x, sample_right_y);
        float energy_mid_mid = energy_function(var_x[i], var_y[i]);
        float energy_aft = energy_function(aft_x, aft_y);
        float energy_fore = energy_function(fore_x, fore_y);
        float work_left = energy_mid_left - energy_fore;
        float work_right = energy_aft - energy_mid_right;

        // dW = (dW/dSx)dSx + (dW/dSy)dSy
        float dW = work_right - work_left;





        // update position
        new_x[i] = var_x[i] - alpha * dE * directional_x;
        new_y[i] = var_y[i] - alpha * dE * directional_y;
    }

    // FTW: After full pass, copy the updates back.
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = new_x[i];
        var_y[i] = new_y[i];
    }

} // iterateWork()
*/ 

float Variational2D::rebalancePoints2D()
{
    // get the total length of curve
    float prev_x = start_x, prev_y = start_y;
    float curve_length = 0;

    // add distance to each variational point
    for (int i=0; i < n_var_points; i++)
    {
        curve_length += sqrt((var_x[i] - prev_x) * (var_x[i] - prev_x) + (var_y[i] - prev_y) * (var_y[i] - prev_y));

        // update for next iteration
        prev_x = var_x[i]; 
        prev_y = var_y[i];
    }

    // add the last piece
    curve_length += sqrt((end_x - prev_x) * (end_x - prev_x) + (end_y - prev_y) * (end_y - prev_y));
    // Cut each new segment along old path to this length
    float new_segment_length = curve_length / (n_var_points + 1);


    // now the respace
    float cursor_x = start_x;
    float cursor_y = start_y;
    float remainder_x = 0;         
    float remainder_y = 0;         
    float remainder = 0;
    int old_point = 0;

    float new_var_x[n_var_points];
    float new_var_y[n_var_points];

    for (int i=0; i<n_var_points; ) 
    {
        // grab old segments until there's enough to cut at least one new one
        // while (remainder_x < new_segment_length)
        while (remainder < new_segment_length)
        {
            float segment_x, segment_y;
            if (old_point < n_var_points) 
            {
                segment_x = var_x[old_point] - cursor_x;
                segment_y = var_y[old_point] - cursor_y;
            }
            else 
            {
                segment_x = end_x - cursor_x;
                segment_y = end_y - cursor_y;
            }
            remainder_x += segment_x;
            remainder_y += segment_y;
            remainder = sqrt(remainder_x*remainder_x + remainder_y*remainder_y);
            cursor_x = var_x[old_point];
            cursor_y = var_y[old_point];
            old_point++;
        }
        
        // now cut new segments until there's not enough left
        // need to cover case equality of segment lengths, so adding it here
        while (new_segment_length <= remainder)
        {
            // get projections of new_segment_length onto x, y directions
            float delta_x = (remainder_x / remainder) * new_segment_length;
            float delta_y = (remainder_y / remainder) * new_segment_length;
            if (i==0) 
            {
                new_var_x[i] = start_x + delta_x;
                new_var_y[i] = start_y + delta_y;
            }
            else
            {
                new_var_x[i] = new_var_x[i-1] + delta_x;
                new_var_y[i] = new_var_y[i-1] + delta_y;
            }
            remainder_x -= delta_x;
            remainder_y -= delta_y;
            remainder = sqrt(remainder_x*remainder_x + remainder_y*remainder_y);
            i++; // move to next new point
        }            
    }

    // check new curve length, add distance to each variational point
    float new_curve_length = 0;
    prev_x = start_x, prev_y = start_y;
    
    for (int i=0; i < n_var_points; i++)
    {
        new_curve_length += sqrt((new_var_x[i] - prev_x) * (new_var_x[i] - prev_x) + (new_var_y[i] - prev_y) * (new_var_y[i] - prev_y));

        // update for next iteration
        prev_x = new_var_x[i]; 
        prev_y = new_var_y[i];
    }

    // add the last piece
    new_curve_length += sqrt((end_x - prev_x) * (end_x - prev_x) + (end_y - prev_y) * (end_y - prev_y));
    
    float shrinkage = new_curve_length / curve_length;

    // copy out
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = new_var_x[i];
        var_y[i] = new_var_y[i];
    }
   
    return shrinkage;
}
