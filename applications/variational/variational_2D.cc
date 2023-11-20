#include <iostream>
#include <stdio.h>
#include <math.h>
#include "variational.hh"

float epsilon;
float sqrt_epsilon;


// prototype only
float energy(float x, float y);


// constructor definition
Variational2D::Variational2D(float _start_x, float _start_y, float _end_x, float _end_y, int _n_iter, int _update, int _n_var_points, float(*_energy_function)(float x, float y))
{
    // grab the passed values
    start_x = _start_x; 
    start_y = _start_y; 
    end_x = _end_x; 
    end_y = _end_y; 
    n_iter  = _n_iter;
    update  = _update;;
    n_var_points = _n_var_points;
    energy_function = _energy_function;

    std::cout << "constructing " << this << std::endl;

    // start constructing 

    var_x = new float[n_var_points];
    var_y = new float[n_var_points];

    // initialize set of points
    
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = start_x + (i + 1) * (end_x - start_x) / (n_var_points + 1);
        var_y[i] = start_y + (i + 1) * (end_y - start_y) / (n_var_points + 1);
        printf("var[i] = {%f,%f}\n", var_x[i], var_y[i]);
    }
}

Variational2D::~Variational2D()
{
    // free any resources allocated
    std::cout << "destructing " << this << std::endl;
    delete var_x, var_y;
}

void Variational2D::printValues()
{
    std::cout << start_x << " " << start_y << std::endl;
    for (int i=0; i<n_var_points; i++)
        std::cout << var_x[i] << " " << var_y[i] << std::endl;
    std::cout << end_x << " " << end_y << std::endl;
}
    


/* everything from here to bottom is the meat

// start of old method. how do I want to return the curve? 
{

    // initialize set of points
    
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = start_x + (i + 1) * (end_x - start_x) / (n_var_points + 1);
        var_y[i] = start_y + (i + 1) * (end_y - start_y) / (n_var_points + 1);
        printf("var[i] = {%f,%f}\n", var_x[i], var_y[i]);
    }

^^^moved to constructor
*/

void Variational2D::iterate()
{    // now iterate
    for (int iter=0; iter < n_iter; iter++)
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

            std::cout << "direction for var point " << i << ": (" << directional_x << ", " << directional_y << ")" << "\tmagnitude: " << sqrt(directional_x * directional_x + directional_y * directional_y) << std::endl;

            // resize the direction vector to machine epsilon
            directional_x *= sqrt_epsilon;
            directional_y *= sqrt_epsilon;

            // sample energy to evaluate derivative
            float sample_left_x = var_x[i] - directional_x;
            float sample_left_y = var_y[i] - directional_y;
            float sample_right_x = var_x[i] + directional_x;
            float sample_right_y = var_y[i] + directional_y;
            float energy_left = energy(sample_left_x, sample_left_y);
            float energy_right = energy(sample_right_x, sample_right_y);

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

        printf("iteration: %d\n", iter);

        printf("%f %f\n", start_x, start_y);
        for (int i=0; i<n_var_points; i++)
        {
        printf("%f %f\n", var_x[i], var_y[i]);
        }
        printf("%f %f\n", end_x, end_y);

    }    
}


/* more of the original
    // now iterate
    for (int iter=0; iter < n_iter; iter++)
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

            std::cout << "direction for var point " << i << ": (" << directional_x << ", " << directional_y << ")" << "\tmagnitude: " << sqrt(directional_x * directional_x + directional_y * directional_y) << std::endl;

            // resize the direction vector to machine epsilon
            directional_x *= sqrt_epsilon;
            directional_y *= sqrt_epsilon;

            // sample energy to evaluate derivative
            float sample_left_x = var_x[i] - directional_x;
            float sample_left_y = var_y[i] - directional_y;
            float sample_right_x = var_x[i] + directional_x;
            float sample_right_y = var_y[i] + directional_y;
            float energy_left = energy(sample_left_x, sample_left_y);
            float energy_right = energy(sample_right_x, sample_right_y);

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

        printf("iteration: %d\n", iter);

        printf("%f %f\n", start_x, start_y);
        for (int i=0; i<n_var_points; i++)
        {
        printf("%f %f\n", var_x[i], var_y[i]);
        }
        printf("%f %f\n", end_x, end_y);

>>>end of iterate function
*/



/* call the rebalance function

        // call the update to re-space points:
        if (update && (iter % update == 0))
        {
            std::cout << "rebalancing at iter = " << iter << std::endl;
            float shrinkage = rebalance_points_2D(start_x, start_y, end_x, end_y, var_x, var_y, n_var_points);
printf("shrinkage: %f\n", shrinkage);
            for (int i=0; i<n_var_points; i++) printf("%f %f\n", var_x[i], var_y[i]);
        }
    } // next iter
//end of method
}

end of meat */ 
