#include <math.h>
#include <iostream>
#include "rebalance.hh"

void rebalance_points(float _start_x, float _start_y, float _end_x, float _end_y, float* _var_x, float* _var_y, int _n_var_points)
{
    // get the total length of curve
    float prev_x = _start_x, prev_y = _start_y;
    float curve_length = 0;

    // add distance to each variational point
    for (int i=0; i < _n_var_points; i++)
    {
        curve_length += sqrt((_var_x[i] - prev_x) * (_var_x[i] - prev_x) + (_var_y[i] - prev_y) * (_var_y[i] - prev_y));

        // update for next iteration
        prev_x = _var_x[i]; 
        prev_y = _var_y[i];
    }

    // add the last piece
    curve_length += sqrt((_end_x - prev_x) * (_end_x - prev_x) + (_end_y - prev_y) * (_end_y - prev_y));


    // now redistribute the points 
    float new_segment_length = curve_length / (_n_var_points + 1);

    // old_segment vector is relative
    float old_segment_x = 0;
    float old_segment_y = 0;

    // new_segment vector is relative
    float new_segment_x = 0;
    float new_segment_y = 0;

    // cursor sits on old path at last cut, value is absolute
    float cursor_x = _start_x;
    float cursor_y = _start_y;

    float new_var_x[_n_var_points], new_var_y[_n_var_points]; // new points

    int new_point_index = 0;

    // outer for loop over old points/segments, inner while loop over old segment
    for (int old_point_index = 0; old_point_index < _n_var_points; old_point_index++)
    {
        // grab the segment and start working on it.
        old_segment_x += (_var_x[old_point_index] - cursor_x);
        old_segment_y += (_var_y[old_point_index] - cursor_y);
        float old_segment_length = sqrt(old_segment_x*old_segment_x + old_segment_y*old_segment_y);

        // cut new segments from old until there's not enough left; zero, once, or multiple
        while (old_segment_length > new_segment_length)
        {
            // Cut the piece from the old segment, and add to the new
            float delta_x = new_segment_length * (old_segment_x / old_segment_length);
            float delta_y = new_segment_length * (old_segment_y / old_segment_length);
            old_segment_x -= delta_x;
            old_segment_y -= delta_y;
            // update the vector length so we keep getting the right value when we normalize
            old_segment_length = sqrt(old_segment_x*old_segment_x + old_segment_y*old_segment_y);

            // add the remainder to the new segment
            new_segment_x += delta_x;
            new_segment_y += delta_y;

            // mark the new point
            new_var_x[new_point_index] = cursor_x + new_segment_x;
            new_var_y[new_point_index] = cursor_y + new_segment_y;

            printf("%d: %f, %f\n", new_point_index, new_var_x[new_point_index], new_var_y[new_point_index]);

            // increment the new point counter and move the cursor, check bounds
            if (new_point_index++ > _n_var_points)
            {
                std::cout << "too many (" << new_point_index << ") new points! exiting.\n";
                exit(1);
            }
            cursor_x += new_segment_x;
            cursor_y += new_segment_y;
            // new_segment_x = 0;
            // new_segment_y = 0;
        }

    }

    // copy out
    for (int i=0; i<_n_var_points; i++)
    {
        _var_x[i] = new_var_x[i];
        _var_y[i] = new_var_y[i];
    }
   
}

