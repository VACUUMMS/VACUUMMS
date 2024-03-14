#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vacuumms/types.h>

#include <vacuumms/variational/variational.hh>
#include <vacuumms/variational/quaternion.hh>

// These are provided and initialized in constants.o
//extern vacuumms_float machine_epsilon;
//extern vacuumms_float machine_sqrt_epsilon;


void Variational3D::init(vacuumms_float _start_x, 
                    vacuumms_float _start_y, 
                    vacuumms_float _start_z, 
                    vacuumms_float _end_x, 
                    vacuumms_float _end_y, 
                    vacuumms_float _end_z, 
                    vacuumms_float _sigma, 
                    vacuumms_float _epsilon, 
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

    var_x = new vacuumms_float[n_var_points];
    var_y = new vacuumms_float[n_var_points];
    var_z = new vacuumms_float[n_var_points];

    // initialize set of points
    
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = start_x + (i + 1) * (end_x - start_x) / (n_var_points + 1);
        var_y[i] = start_y + (i + 1) * (end_y - start_y) / (n_var_points + 1);
        var_z[i] = start_z + (i + 1) * (end_z - start_z) / (n_var_points + 1);
    }
}


Variational3D::Variational3D(vacuumms_float _start_x, 
                             vacuumms_float _start_y, 
                             vacuumms_float _start_z, 
                             vacuumms_float _end_x, 
                             vacuumms_float _end_y, 
                             vacuumms_float _end_z, 
                             vacuumms_float _sigma, 
                             vacuumms_float _epsilon, 
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
    delete var_x, var_y, var_z;
}


void Variational3D::printValues()
{
    std::cout << "#####" << start_x << " " << start_y << " " << start_z << std::endl;
    for (int i=0; i<n_var_points; i++)
        std::cout << "#####" << var_x[i] << " " << var_y[i] << " " << var_z[i] << std::endl;
    std::cout << "#####" << end_x << " " << end_y << " " << end_z << std::endl;
}

    
void Variational3D::setAlpha(vacuumms_float _alpha)
{
    alpha = _alpha;
}


void Variational3D::setAlphaMax(vacuumms_float _alpha_max)
{
    alpha_max = _alpha_max;
}


void Variational3D::setDeltaMax(vacuumms_float _delta_max)
{
    delta_max = _delta_max;
}


/* note that iterate() does not adjust alpha value automatically,
 * but adapativeIterateAndUpdate() does.
 */
void Variational3D::iterate()
{
    debug = getenv("VACUUMMS_DEBUG");

    // Space for the perturbed curve values
    vacuumms_float new_x[n_var_points], new_y[n_var_points], new_z[n_var_points];

    // looking along direction of line, point before and point after
    vacuumms_float fore_x, aft_x;
    vacuumms_float fore_y, aft_y;
    vacuumms_float fore_z, aft_z;

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
        vacuumms_float tangent_x = aft_x - fore_x;
        vacuumms_float tangent_y = aft_y - fore_y;
        vacuumms_float tangent_z = aft_z - fore_z;
        
        if (debug != NULL) printf("got tangent vector (%f, %f, %f)\n", tangent_x, tangent_y, tangent_z);

        vacuumms_float u_x, u_y, u_z; // unit axis of rotation, to be extracted from tangent
        vacuumms_float theta = extract_axis(0, 0, 1, tangent_x, tangent_y, tangent_z, &u_x, &u_y, &u_z);

        if (debug != NULL) printf("got theta = %f and axis u = (%f, %f, %f)\n ", theta, u_x, u_y, u_z);

        // variables to receive the values of the rotated i and j unit vectors
        vacuumms_float i_x, i_y, i_z; 
        rotate_vector(1, 0, 0, theta, u_x, u_y, u_z, &i_x, &i_y, &i_z);

        if (debug != NULL) printf("i rotates to (%f, %f, %f)\n", i_x, i_y, i_z);

        vacuumms_float j_x, j_y, j_z; 
        rotate_vector(0, 1, 0, theta, u_x, u_y, u_z, &j_x, &j_y, &j_z);

        if (debug != NULL) printf("j rotates to (%f, %f, %f)\n", j_x, j_y, j_z);

        // The directional vectors are ostensibly normalized.
        // Resize the directionals as deltas:

        i_x *= sqrt_machine_epsilon;
        i_y *= sqrt_machine_epsilon;
        i_z *= sqrt_machine_epsilon;

        j_x *= sqrt_machine_epsilon;
        j_y *= sqrt_machine_epsilon;
        j_z *= sqrt_machine_epsilon;

        if (debug != NULL) 
        {
            printf("i resized to (%f, %f, %f)\n", i_x, i_y, i_z);
            printf("j resized to (%f, %f, %f)\n", i_x, i_y, i_z);
        }

        // energy at var[i], and at directional points
        vacuumms_float energy_center, energy_i, energy_j;
        
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

        // dE = (dE/dx)dx + (dE/dy)dy

        vacuumms_float dE_i = energy_i - energy_center;
        vacuumms_float dE_j = energy_j - energy_center;

        // update position

        vacuumms_float delta_x = - alpha * (dE_i * i_x + dE_j * j_x);
        vacuumms_float delta_y = - alpha * (dE_i * i_y + dE_j * j_y);
        vacuumms_float delta_z = - alpha * (dE_i * i_z + dE_j * j_z);
        vacuumms_float delta_sq = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
        vacuumms_float delta = sqrt(delta_sq);

        if (debug != NULL)
        {  
            printf("Got dE_i = %f, dE_j = %f, delta = (%f, %f, %f) |delta| = %f\n", 
            dE_i, dE_j, delta_x, delta_y, delta_z, delta);
        }

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

/*########################################## rebalance ###########################################*/

vacuumms_float Variational3D::rebalancePoints3D()
{
    debug = getenv("VACUUMMS_DEBUG");

    vacuumms_float original_curve_length = calculateCurveLength(
        start_x, start_y, start_z, 
        end_x, end_y, end_z, 
        var_x, var_y, var_z
    );

    vacuumms_float new_forward_var_x[n_var_points];
    vacuumms_float new_forward_var_y[n_var_points];
    vacuumms_float new_forward_var_z[n_var_points];

    int status;
    status = respaceKernel(start_x, start_y, start_z, 
                           end_x, end_y, end_z, 
                           var_x, var_y, var_z, 
                           new_forward_var_x, new_forward_var_y, new_forward_var_z);
    if (status != 0) exit(status);

    vacuumms_float forward_curve_length = calculateCurveLength(
        start_x, start_y, start_z, 
        end_x, end_y, end_z, 
        new_forward_var_x, new_forward_var_y, new_forward_var_z
    );

    // now run it backward

    vacuumms_float backward_var_x_in[n_var_points];
    vacuumms_float backward_var_y_in[n_var_points];
    vacuumms_float backward_var_z_in[n_var_points];
    
    vacuumms_float backward_var_x_out[n_var_points];
    vacuumms_float backward_var_y_out[n_var_points];
    vacuumms_float backward_var_z_out[n_var_points];
    
    // flip the curve
    for (int i=0; i<n_var_points; i++)
    {
        backward_var_x_in[i] = var_x[n_var_points - i - 1];
        backward_var_y_in[i] = var_y[n_var_points - i - 1];
        backward_var_z_in[i] = var_z[n_var_points - i - 1];
    }

    status = respaceKernel(end_x, end_y, end_z, 
                           start_x, start_y, start_z, 
                           backward_var_x_in, backward_var_y_in, backward_var_z_in, 
                           backward_var_x_out, backward_var_y_out, backward_var_z_out);
    if (status != 0) exit(status);

    // and flip it back
    vacuumms_float new_backward_var_x[n_var_points];
    vacuumms_float new_backward_var_y[n_var_points];
    vacuumms_float new_backward_var_z[n_var_points];
    
    for (int i=0; i<n_var_points; i++)
    {
        new_backward_var_x[i] = backward_var_x_out[n_var_points - i - 1];
        new_backward_var_y[i] = backward_var_y_out[n_var_points - i - 1];
        new_backward_var_z[i] = backward_var_z_out[n_var_points - i - 1];
    }

    vacuumms_float backward_curve_length = calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, new_forward_var_x, new_forward_var_y, new_forward_var_z);

    // should match
    if (debug != NULL) printf("original backward curve_length = %f\n", backward_curve_length);

    // combine the results
    vacuumms_float new_var_x[n_var_points];
    vacuumms_float new_var_y[n_var_points];
    vacuumms_float new_var_z[n_var_points];

    for (int i=0; i<n_var_points; i++) 
    {
        new_var_x[i] = 0.5 * (new_forward_var_x[i] + new_backward_var_x[i]);
        new_var_y[i] = 0.5 * (new_forward_var_y[i] + new_backward_var_y[i]);
        new_var_z[i] = 0.5 * (new_forward_var_z[i] + new_backward_var_z[i]);
    }

    // check new combined curve length, add distance to each variational point to get shrinkage
    vacuumms_float new_curve_length = calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, new_var_x, new_var_y, new_var_z);

    vacuumms_float shrinkage = new_curve_length / original_curve_length;

    // copy out
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = new_var_x[i];
        var_y[i] = new_var_y[i];
        var_z[i] = new_var_z[i];
    }
   
    return shrinkage;
}


// return value is 0 if successful, 1 if not successful
int Variational3D::respaceKernel(vacuumms_float _start_x, vacuumms_float _start_y, vacuumms_float _start_z, 
                                 vacuumms_float _end_x, vacuumms_float _end_y, vacuumms_float _end_z,
                                 vacuumms_float _var_x[], vacuumms_float _var_y[], vacuumms_float _var_z[], 
                                 vacuumms_float _new_var_x[], vacuumms_float _new_var_y[], vacuumms_float _new_var_z[])
{
    debug = getenv("VACUUMMS_DEBUG");

    vacuumms_float curve_length = calculateCurveLength(_start_x, _start_y, _start_z, _end_x, _end_y, _end_z, _var_x, _var_y, _var_z);
    // Cut each new segment along old path to this length
    vacuumms_float new_segment_length = curve_length / (/*_*/n_var_points + 1.0);

    // now the respace
    vacuumms_float cursor_x = _start_x;
    vacuumms_float cursor_y = _start_y;
    vacuumms_float cursor_z = _start_z;
    vacuumms_float remainder_x = 0;         
    vacuumms_float remainder_y = 0;         
    vacuumms_float remainder_z = 0;         
    vacuumms_float remainder = 0;
    int old_point = 0;

    int all_segments_added = 0;
    for (int new_point=0; new_point</*_*/n_var_points; ) 
    {
        int something_happened = 0;
        // grab old segments until there's enough curve to cut at least one new segment
        while ((remainder < new_segment_length) && (old_point <= /*_*/n_var_points))
        {
            vacuumms_float segment_x, segment_y, segment_z;
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
        while ((new_segment_length <= remainder) && (new_point < n_var_points)) 
        {
            // get projections of new_segment_length onto x, y directions
            vacuumms_float delta_x = (remainder_x / remainder) * new_segment_length;
            vacuumms_float delta_y = (remainder_y / remainder) * new_segment_length;
            vacuumms_float delta_z = (remainder_z / remainder) * new_segment_length;

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
            printf("or adding more variational points (current n_var_points = %d) \n", n_var_points);
            return 1;
        }
    }

    return 0;
} // respaceKernel()


vacuumms_float Variational3D::calculateCurveLength(vacuumms_float _start_x, vacuumms_float _start_y, vacuumms_float _start_z, vacuumms_float _end_x, vacuumms_float _end_y, vacuumms_float _end_z, vacuumms_float _var_x[], vacuumms_float _var_y[], vacuumms_float _var_z[])
{
    debug = getenv("VACUUMMS_DEBUG");

    // get the total length of curve
    vacuumms_float prev_x = _start_x;
    vacuumms_float prev_y = _start_y;
    vacuumms_float prev_z = _start_z;
    vacuumms_float curve_length = 0;

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


/* note that iterate() does not adjust alpha value automatically,
 * but adapativeIterateAndUpdate() does.
 */
vacuumms_float Variational3D::adaptiveIterateAndUpdate()
{
    debug = getenv("VACUUMMS_DEBUG");
    int wedge_count = 0;

    vacuumms_float curve_length = calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, var_x, var_y, var_z);
    vacuumms_float new_x[n_var_points], new_y[n_var_points], new_z[n_var_points];

attempt_iteration:

    if (debug != NULL) printf("Attempting iteration with alpha = %f\n", alpha);
    
    vacuumms_float fore_x, aft_x;
    vacuumms_float fore_y, aft_y;
    vacuumms_float fore_z, aft_z;

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
        vacuumms_float tangent_x = aft_x - fore_x;
        vacuumms_float tangent_y = aft_y - fore_y;
        vacuumms_float tangent_z = aft_z - fore_z;

        vacuumms_float u_x, u_y, u_z; // unit axis of rotation, to be extracted from tangent
        vacuumms_float theta = extract_axis(0, 0, 1, tangent_x, tangent_y, tangent_z, &u_x, &u_y, &u_z);

        // variables to receive the values of the rotated i and j unit vectors
        vacuumms_float i_x, i_y, i_z; 
        rotate_vector(1, 0, 0, theta, u_x, u_y, u_z, &i_x, &i_y, &i_z);
        vacuumms_float j_x, j_y, j_z; 
        rotate_vector(0, 1, 0, theta, u_x, u_y, u_z, &j_x, &j_y, &j_z);

        // resize the directionals as deltas

        i_x *= sqrt_machine_epsilon;
        i_y *= sqrt_machine_epsilon;
        i_z *= sqrt_machine_epsilon;

        j_x *= sqrt_machine_epsilon;
        j_y *= sqrt_machine_epsilon;
        j_z *= sqrt_machine_epsilon;

        // energy at var[i], and at directional points
        vacuumms_float energy_center, energy_i, energy_j;
        
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

        vacuumms_float dE_i = energy_i - energy_center;
        vacuumms_float dE_j = energy_j - energy_center;

        vacuumms_float delta_x = - alpha * (dE_i * i_x + dE_j * j_x);
        vacuumms_float delta_y = - alpha * (dE_i * i_y + dE_j * j_y);
        vacuumms_float delta_z = - alpha * (dE_i * i_z + dE_j * j_z);
        vacuumms_float delta_sq = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
        vacuumms_float delta = sqrt(delta_sq);

        if (debug != NULL) 
        {
            printf("Got dE_i = %f, dE_j = %f, delta = (%f, %f, %f) |delta| = %f\n", 
                   dE_i, dE_j, delta_x, delta_y, delta_z, delta); 
        }

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

// end of iterate portion; have new_x, new_y, and new_z to carry forward

//########################################## rebalance ###########################################//

    // run it forward

    vacuumms_float new_forward_var_x[n_var_points];
    vacuumms_float new_forward_var_y[n_var_points];
    vacuumms_float new_forward_var_z[n_var_points];

    int forward_status = respaceKernel(start_x, start_y, start_z, end_x, end_y, end_z, new_x, new_y, new_z, new_forward_var_x, new_forward_var_y, new_forward_var_z);

    if (forward_status != 0) 
    {
        wedge_count++;
        if (wedge_count > 100)
        {
            fprintf(stderr, "vacuumms/variational: wedge count > 100 in respaceKernel, giving up.\n");
            exit(1);
        }
        fprintf(stderr, "forward respace failed, shrinking alpha = %f to %f and re-running.\n", alpha, alpha/=beta);
        goto attempt_iteration;
    }

    // now run it backward

    vacuumms_float backward_var_x_in[n_var_points];
    vacuumms_float backward_var_y_in[n_var_points];
    vacuumms_float backward_var_z_in[n_var_points];
    
    vacuumms_float backward_var_x_out[n_var_points];
    vacuumms_float backward_var_y_out[n_var_points];
    vacuumms_float backward_var_z_out[n_var_points];
    
    // flip the curve
    for (int i=0; i<n_var_points; i++)
    {
        backward_var_x_in[i] = new_x[n_var_points - i - 1];
        backward_var_y_in[i] = new_y[n_var_points - i - 1];
        backward_var_z_in[i] = new_z[n_var_points - i - 1];
    }

    int backward_status = respaceKernel(end_x, end_y, end_z, start_x, start_y, start_z, backward_var_x_in, backward_var_y_in, backward_var_z_in, backward_var_x_out, backward_var_y_out, backward_var_z_out);

    if (backward_status != 0) 
    {
        wedge_count++;
        if (wedge_count > 100)
        {
            fprintf(stderr, "vacuumms/variational: wedge count > 100 in respaceKernel, giving up.\n");
            exit(1);
        }
        fprintf(stderr, "backward respace failed, shrinking alpha = %f to %f and re-running.\n", alpha, alpha/=beta);
        goto attempt_iteration;
    }

    // and flip it back
    vacuumms_float new_backward_var_x[n_var_points];
    vacuumms_float new_backward_var_y[n_var_points];
    vacuumms_float new_backward_var_z[n_var_points];
    
    for (int i=0; i<n_var_points; i++)
    {
        new_backward_var_x[i] = backward_var_x_out[n_var_points - i - 1];
        new_backward_var_y[i] = backward_var_y_out[n_var_points - i - 1];
        new_backward_var_z[i] = backward_var_z_out[n_var_points - i - 1];
    }

    // TO DO: Add some sort of assert(forward_curve_length ~= backward_curve_length) ?? 

    // combine the forward and backward results
    vacuumms_float respace_var_x[n_var_points];
    vacuumms_float respace_var_y[n_var_points];
    vacuumms_float respace_var_z[n_var_points];

    for (int i=0; i<n_var_points; i++) 
    {
        respace_var_x[i] = 0.5 * (new_forward_var_x[i] + new_backward_var_x[i]);
        respace_var_y[i] = 0.5 * (new_forward_var_y[i] + new_backward_var_y[i]);
        respace_var_z[i] = 0.5 * (new_forward_var_z[i] + new_backward_var_z[i]);
    }

    vacuumms_float new_curve_length = calculateCurveLength(start_x, start_y, start_z, end_x, end_y, end_z, respace_var_x, respace_var_y, respace_var_z);

    vacuumms_float shrinkage = new_curve_length / curve_length;

    // copy out
    for (int i=0; i<n_var_points; i++)
    {
        var_x[i] = respace_var_x[i];
        var_y[i] = respace_var_y[i];
        var_z[i] = respace_var_z[i];
    }

    if (debug != NULL) printf("successful update and rebalance. ");

    // With success, increase alpha by factor of beta
    vacuumms_float new_alpha = alpha * beta;

    if (new_alpha < alpha_max) // success, so increase alpha
    {
        if (debug != NULL) printf("Increasing alpha: %f * %f = %f\n", alpha, beta, new_alpha); 
        alpha = new_alpha;
    }
    if (debug != NULL) printf("\n");
   
    return shrinkage;

} // adaptiveIterateAndUpdate()


vacuumms_float* Variational3D::getX()
{
    return var_x;
}


vacuumms_float* Variational3D::getY()
{
    return var_y;
}


vacuumms_float* Variational3D::getZ()
{
    return var_z;
}
