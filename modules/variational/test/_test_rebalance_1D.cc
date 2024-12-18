#include <iostream>
#include "rebalance.hh"


void initialize(float var_x[], int n_var_points)
{
    for (int i=0; i<n_var_points; i++) var_x[i] = (float)(i+1)*float((i+1));
}

void print_points(float var_x[], int n_var_points)
{
    for (int i=0; i<n_var_points; i++) std::cout << var_x[i] << std::endl;
}

/*
void respace_points(float var_x[], int n_var_points)
{
    float new_var_x[n_var_points];

    // get length of the curve by summing segments
    float curve_length=0;
    float cursor_x = start_x;
    for (int i=0; i<n_var_points; i++) 
    {
        curve_length += var_x[i] - cursor_x;
        cursor_x = var_x[i];
    }
    curve_length += end_x - cursor_x;
    float new_segment_length = curve_length / (n_var_points +1);
    
    std::cout << "got curve_length = " << curve_length << std::endl;

    // now the respace
    cursor_x = start_x;
    float remainder_x = 0;         
    int old_point = 0;

    for (int i=0; i<n_var_points; ) 
    {
        // grab old segments until there's enough to cut at least one new one
        while (remainder_x < new_segment_length)
        {
            float segment_x;
            if (old_point < n_var_points) segment_x = var_x[old_point] - cursor_x;
            else segment_x = end_x - cursor_x;
            remainder_x += segment_x;
            cursor_x = var_x[old_point++];
        }
        
        // now cut new segments until there's not enough left
        while (new_segment_length < remainder_x)
        {
            if (i==0) new_var_x[i] = start_x + new_segment_length;
            else new_var_x[i] = new_var_x[i-1] + new_segment_length;
            remainder_x -= new_segment_length;
            i++; // move to next new point
        }            
    }

    // copy out
    for (int i=0; i<n_var_points; i++) var_x[i] = new_var_x[i];

}
*/


int main()
{
    float start_x = 0.0;
    float end_x = 100.0;
    int n_var_points = 9;
    float var_x[n_var_points];

    initialize(var_x, n_var_points);
    print_points(var_x, n_var_points);
    printf("rebalancing...\n");
    float shrinkage = rebalance_points_1D(start_x, end_x, var_x, n_var_points);
    print_points(var_x, n_var_points);
    
    printf("shrinkage: %f\n", shrinkage);
}

      
