#include <iostream>
#include <stdio.h>
#include <math.h>
#include "rebalance.hh"

extern "C"
{
#include <ftw_param.h>
}; 

float epsilon =  5.96e-08;
float sqrt_epsilon = 0.000244131112315;

float points_x[] = {-0.0149969,0.0149969};
float points_y[] = {1.0,-1.0};
int n_points = 2;
int n_iter = 1;
int n_var_points = 25;
int update = 1;
//float var_x[n_var_points], var_y[n_var_points];
float *var_x, *var_y;

float end_x = 2;
float end_y = 0;
float start_x = -2;
float start_y = 0;

// Calculate energy at this point
float energy(float x, float y)
{
    float total = 0.0;
    for (int i=0; i<n_points; i++)
    {
        float r_sq = (points_x[i] - x) * (points_x[i] - x) + (points_y[i] - y) * (points_y[i] - y);
        float r_6 = r_sq * r_sq * r_sq;
        float r_12 = r_6 * r_6;
        total += (1.0/r_12 - 1.0/r_6);
    }
    return total;    
}

int main(int argc, char** argv)
{
    setCommandLineParameters(argc, argv);
    getIntParam((char*)"-n_iter", &n_iter);
    getIntParam((char*)"-update", &update);
    getIntParam((char*)"-n_var_points", &n_var_points);

    var_x = new float[n_var_points];
    var_y = new float[n_var_points];

    // initialize set of points
    
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = start_x + (i + 1) * (end_x - start_x) / (n_var_points + 1);
        var_y[i] = start_y + (i + 1) * (end_y - start_y) / (n_var_points + 1);
        printf("var[i] = {%f,%f}\n", var_x[i], var_y[i]);
    }

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
            std::cout << "dE = " << dE << "\n";

            // update position
            new_x[i] = var_x[i] - alpha * dE * directional_x;
            new_y[i] = var_y[i] - alpha * dE * directional_y;

/*
            if (dE < 0)
            {
                std::cout << "nudging left.\n";
                new_x[i] = var_x[i] - directional_x;
                new_y[i] = var_y[i] - directional_y;
            }
            else if (dE > 0)
            {
                std::cout << "nudging right.\n";
                new_x[i] = var_x[i] + directional_x;
                new_y[i] = var_y[i] + directional_y;
            }
            else std::cout << "no change.\n";
*/

            // generate perturbation and calculate energy at perturbed points

/*
            float perturb_x = var_x[i] - delta_y;
            float perturb_y = var_y[i] + delta_x;

            float energy_0 = energy(var_x[i], var_y[i]);
            float energy_1 = energy(var_x[i] + perturb_x, var_y[i] + perturb_y);
     
            // descent?
            // or bump by epsilon in direction of normal

            float delta_energy = energy_1 - energy_0; 
            
    //        std::cout << "energy_0 = " << energy_0 << "\tenergy_1 = " << energy_1 << std::endl;
            std::cout << "delta_energy = " << delta_energy << std::endl;

            // now make the update
            var_x[i] -= delta_y * delta_energy ; 
            var_y[i] += delta_x * delta_energy ; 
*/
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

        // call the update to re-space points:
        if (update && (iter % update == 0))
        {
            std::cout << "rebalancing at iter = " << iter << std::endl;
            rebalance_points(start_x, start_y, end_x, end_y, var_x, var_y, n_var_points);
            for (int i=0; i<n_var_points; i++) printf("%f %f\n", var_x[i], var_y[i]);
        }
    }
}
