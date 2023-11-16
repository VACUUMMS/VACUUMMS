#include <math.h>
#include "rebalance.hh"

float rebalance_points_3D(float _start_x, float _start_y, float _start_z, float _end_x, float _end_y, float _end_z, float* _var_x, float* _var_y, float* _var_z, int _n_var_points)
{
    // get the total length of curve
    float prev_x = _start_x, prev_y = _start_y, prev_z = _start_z;
    float curve_length = 0;

    // add distance to each variational point
    for (int i=0; i < _n_var_points; i++)
    {
        curve_length += sqrt((_var_x[i] - prev_x) * (_var_x[i] - prev_x) + (_var_y[i] - prev_y) * (_var_y[i] - prev_y) + (_var_z[i] - prev_z) * (_var_z[i] - prev_z));

        // update for next iteration
        prev_x = _var_x[i]; 
        prev_y = _var_y[i];
        prev_z = _var_z[i];
    }

    // add the last piece
    curve_length += sqrt((_end_x - prev_x) * (_end_x - prev_x + (_end_y - prev_y) * (_end_y - prev_y)) + (_end_z - prev_z) * (_end_z - prev_z));
    // Cut each new segment along old path to this length
    float new_segment_length = curve_length / (_n_var_points + 1);

    // now the respace
    float cursor_x = _start_x;
    float cursor_y = _start_y;
    float cursor_z = _start_z;
    float remainder_x = 0;         
    float remainder_y = 0;         
    float remainder_z = 0;         
    float remainder = 0;
    int old_point = 0;

    float new_var_x[_n_var_points];
    float new_var_y[_n_var_points];
    float new_var_z[_n_var_points];

    for (int i=0; i<_n_var_points; ) 
    {
        // grab old segments until there's enough to cut at least one new one
        while (remainder < new_segment_length)
        {
            float segment_x, segment_y, segment_z;
            if (old_point < _n_var_points) 
            {
                segment_x = _var_x[old_point] - cursor_x;
                segment_y = _var_y[old_point] - cursor_y;
                segment_z = _var_z[old_point] - cursor_z;
            }
            else 
            {
                segment_x = _end_x - cursor_x;
                segment_y = _end_y - cursor_y;
                segment_z = _end_z - cursor_z;
            }
            remainder_x += segment_x;
            remainder_y += segment_y;
            remainder_z += segment_z;
            remainder = sqrt(remainder_x*remainder_x + remainder_y*remainder_y + remainder_z*remainder_z);
            cursor_x = _var_x[old_point];
            cursor_y = _var_y[old_point];
            cursor_z = _var_z[old_point];
            old_point++;
        }
        
        // now cut new segments until there's not enough left
        // need to cover case equality of segment lengths, so adding it here
        while (new_segment_length <= remainder)
        {
            // get projections of new_segment_length onto x, y directions
            float delta_x = (remainder_x / remainder) * new_segment_length;
            float delta_y = (remainder_y / remainder) * new_segment_length;
            float delta_z = (remainder_z / remainder) * new_segment_length;
            if (i==0) 
            {
                new_var_x[i] = _start_x + delta_x;
                new_var_y[i] = _start_y + delta_y;
                new_var_z[i] = _start_z + delta_z;
            }
            else
            {
                new_var_x[i] = new_var_x[i-1] + delta_x;
                new_var_y[i] = new_var_y[i-1] + delta_y;
                new_var_z[i] = new_var_z[i-1] + delta_z;
            }
            remainder_x -= delta_x;
            remainder_y -= delta_y;
            remainder_z -= delta_z;
            remainder = sqrt(remainder_x*remainder_x + remainder_y*remainder_y + remainder_z*remainder_z);
            i++; // move to next new point
        }            
    }

    // check new curve length, add distance to each variational point
    float new_curve_length = 0;
    prev_x = _start_x, prev_y = _start_y; prev_z = _start_z;
    
    for (int i=0; i < _n_var_points; i++)
    {
        new_curve_length += sqrt((new_var_x[i] - prev_x) * (new_var_x[i] - prev_x) + (new_var_y[i] - prev_y) * (new_var_y[i] - prev_y) + (new_var_z[i] - prev_z) * (new_var_z[i] - prev_z));

        // update for next iteration
        prev_x = new_var_x[i]; 
        prev_y = new_var_y[i];
        prev_z = new_var_z[i];
    }

    // add the last piece
    new_curve_length += sqrt((_end_x - prev_x) * (_end_x - prev_x) + (_end_y - prev_y) * (_end_y - prev_y) + (_end_z - prev_z) * (_end_z - prev_z));
    
    float shrinkage = new_curve_length / curve_length;

    // copy out
    for (int i=0; i<_n_var_points; i++)
    {
        _var_x[i] = new_var_x[i];
        _var_y[i] = new_var_y[i];
        _var_z[i] = new_var_z[i];
    }
   
    return shrinkage;
}

float rebalance_points_2D(float _start_x, float _start_y, float _end_x, float _end_y, float* _var_x, float* _var_y, int _n_var_points)
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
    // Cut each new segment along old path to this length
    float new_segment_length = curve_length / (_n_var_points + 1);


    // now the respace
    float cursor_x = _start_x;
    float cursor_y = _start_y;
    float remainder_x = 0;         
    float remainder_y = 0;         
    float remainder = 0;
    int old_point = 0;

    float new_var_x[_n_var_points];
    float new_var_y[_n_var_points];

    for (int i=0; i<_n_var_points; ) 
    {
        // grab old segments until there's enough to cut at least one new one
        // while (remainder_x < new_segment_length)
        while (remainder < new_segment_length)
        {
            float segment_x, segment_y;
            if (old_point < _n_var_points) 
            {
                segment_x = _var_x[old_point] - cursor_x;
                segment_y = _var_y[old_point] - cursor_y;
            }
            else 
            {
                segment_x = _end_x - cursor_x;
                segment_y = _end_y - cursor_y;
            }
            remainder_x += segment_x;
            remainder_y += segment_y;
            remainder = sqrt(remainder_x*remainder_x + remainder_y*remainder_y);
            cursor_x = _var_x[old_point];
            cursor_y = _var_y[old_point];
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
                new_var_x[i] = _start_x + delta_x;
                new_var_y[i] = _start_y + delta_y;
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
    prev_x = _start_x, prev_y = _start_y;
    
    for (int i=0; i < _n_var_points; i++)
    {
        new_curve_length += sqrt((new_var_x[i] - prev_x) * (new_var_x[i] - prev_x) + (new_var_y[i] - prev_y) * (new_var_y[i] - prev_y));

        // update for next iteration
        prev_x = new_var_x[i]; 
        prev_y = new_var_y[i];
    }

    // add the last piece
    new_curve_length += sqrt((_end_x - prev_x) * (_end_x - prev_x) + (_end_y - prev_y) * (_end_y - prev_y));
    
    float shrinkage = new_curve_length / curve_length;

    // copy out
    for (int i=0; i<_n_var_points; i++)
    {
        _var_x[i] = new_var_x[i];
        _var_y[i] = new_var_y[i];
    }
   
    return shrinkage;
}

float rebalance_points_1D(float _start_x, float _end_x, float* _var_x, int _n_var_points)
{
    float new_var_x[_n_var_points];

    // get length of the curve by summing segments
    float curve_length=0;
    float cursor_x = _start_x;
    for (int i=0; i<_n_var_points; i++) 
    {
        curve_length += _var_x[i] - cursor_x;
        cursor_x = _var_x[i];
    }
    curve_length += _end_x - cursor_x;
    float new_segment_length = curve_length / (_n_var_points + 1);

    // now the respace
    cursor_x = _start_x;
    float remainder_x = 0;         
    int old_point = 0;

    for (int i=0; i<_n_var_points; ) 
    {
        // grab old segments until there's enough to cut at least one new one
        while (remainder_x < new_segment_length)
        {
            float segment_x;
            if (old_point < _n_var_points) segment_x = _var_x[old_point] - cursor_x;
            else segment_x = _end_x - cursor_x;
            remainder_x += segment_x;
            cursor_x = _var_x[old_point++];
        }
        
        // now cut new segments until there's not enough left
        while (new_segment_length < remainder_x)
        {
            if (i==0) new_var_x[i] = _start_x + new_segment_length;
            else new_var_x[i] = new_var_x[i-1] + new_segment_length;
            remainder_x -= new_segment_length;
            i++; // move to next new point
        }            
    }

    // get length of the new curve by summing segments
    float new_curve_length=0;
    cursor_x = _start_x;
    for (int i=0; i<_n_var_points; i++) 
    {
        new_curve_length += new_var_x[i] - cursor_x;
        cursor_x = new_var_x[i];
    }
    new_curve_length += _end_x - cursor_x;

    // copy out
    for (int i=0; i<_n_var_points; i++) _var_x[i] = new_var_x[i];

    float shrinkage = new_curve_length / curve_length;
    return shrinkage;
}


