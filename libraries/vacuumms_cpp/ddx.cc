/* vacuumms_cpp/ddx.cc */

#include <vacuumms/ddx.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>

#include <vacuumms/limits.h>
#include <vacuumms/rng.h>

#include <math.h>


DDX::DDX(Configuration c, Parameters p) : 
    c{c}, p{p} 
{
}

DDX::DDX(Configuration c) : 
    c{c} 
{
}

DDX::DDX()
{
}


void DDX::setParameters(Parameters _p)
{
    p = _p;
    Operation::setParameters(_p);
}


void DDX::setConfiguration(Configuration _c)
{
    c = _c;
}


Configuration DDX::getConfiguration()
{
    return c;
}


void DDX::setNumberOfSamples(int _number_of_samples)
{
    number_of_samples = _number_of_samples;
}


void DDX::setPrecisionParameter(vacuumms_float _precision_parameter)
{
    precision_parameter = _precision_parameter;
}


void DDX::setSeed(int _seed)
{
    seed = _seed;
}


void DDX::setVerletCutoff(vacuumms_float _verlet_cutoff)
{
    verlet_cutoff = _verlet_cutoff;
}


void DDX::setVerletExtent(int _verlet_extent)
{
    verlet_extent = _verlet_extent;
}


void DDX::setNumberOfSteps(int _n_steps)
{
    n_steps = _n_steps;
}


void DDX::setMinDiameter(vacuumms_float _min_diameter)
{
    min_diameter = _min_diameter;
}


void DDX::Randomize()
{
    seed = 0; // Zero seed will force randomization
//    int val = randomize();
//printf("initialized RNG to %d\n", val);
}


#ifdef BUILD_PYBIND_BINDINGS
pybind11::str DDX::__repr__()
{
    pybind11::str retval;
    retval += c.__repr__();
    retval += p.__repr__();
    retval += result.__repr__();
    return retval;
}
#endif


CavityConfiguration DDX::getResult()
{
    return result;   
}


void DDX::printUsage()
{
    printf("\nDDX usage:\t-box [ 6.0 6.0 6.0 ]\n");
    printf("\t\t-seed [ 1 ]\n");
    printf("\t\t-randomize \n");
    printf("\t\t-characteristic_length [ 1.0 ]\n");
    printf("\t\t-characteristic_energy [ 1.0 ]\n");
    printf("\t\t-precision_parameter [ 0.001 ]\n");
    printf("\t\t-n_steps [ 1000 ] (roughly reciprocal of precision parameter)\n");
    printf("\t\t-show_steps (includes steps taken as final column)\n");
    printf("\t\t-verlet_cutoff [ 100.0 ]\n");
    printf("\t\t-n [ 1 ]\n");
    printf("\t\t-volume_sampling \n");
    printf("\t\t-include_center_energy \n");
    printf("\t\t-min_diameter [ 0.0 ]");
    printf("\n");
}


void DDX::execute()
{
    // Clear old results, if any
    result.reset();

    if (c.getSize() == 0) 
    {
        std::cout << "no atoms in configuration, declining to execute." << std::endl;
        return;
    }
    
    // Get box_dims from config info
    std::vector<vacuumms_float> box_dims = c.getBoxDimensions();
    box_x = box_dims[0];
    box_y = box_dims[1];
    box_z = box_dims[2];

    // override if passed as params.
    p.getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
    result.setBoxDimensions(box_dims);

    if (box_x * box_y * box_z < 0.000001) 
    {
        std::cout << "box dimensions not properly set in configuration, declining to execute." << std::endl;
        return;
    }

    vacuumms_double sq_distance_from_initial_pt;

    verbose = p.getFlagParam((char*)"-verbose");
    p.getIntParam((char*)"-seed", &seed);
    if (seed == 0) seed = randomize();
    else initializeRandomNumberGeneratorTo(seed);

    if ((box_x * box_y * box_z) == 0.0) 
    {
        fprintf(stderr, "Found simulation box volume = 0.0, gracefully exiting.\n");
        return;
    }

    p.getDoubleParam((char*)"-characteristic_length", &characteristic_length);
    p.getDoubleParam((char*)"-characteristic_energy", &characteristic_energy);
    p.getDoubleParam((char*)"-precision_parameter", &precision_parameter);
    p.getDoubleParam((char*)"-verlet_cutoff", &verlet_cutoff);
    p.getIntParam((char*)"-n", &number_of_samples);
    p.getIntParam((char*)"-n_steps", &n_steps);
    volume_sampling = p.getFlagParam((char*)"-volume_sampling");
    include_center_energy = p.getFlagParam((char*)"-include_center_energy");
    show_steps = p.getFlagParam((char*)"-show_steps");
    p.getDoubleParam((char*)"-min_diameter", &min_diameter);

    if (p.getFlagParam((char*)"-usage")) printUsage();

    // was  loadConfiguration(); 
    // now just copy over the records and run the old algorithm
    for (int i=0; i<c.getSize(); i++) 
    {
        ConfigurationRecord r = c.recordAt(i);
        x[i] = r.x;
        y[i] = r.y;
        z[i] = r.z;
        sigma[i] = r.sigma;
        epsilon[i] = r.epsilon;
    }
  
    number_of_molecules = c.getSize();
  
    int remaining_samples = number_of_samples;
    while (remaining_samples-- > 0)
    {
        generateTestPoint();
        while (calculateEnergy(0.0) > 0) generateTestPoint();
    
        findEnergyMinimum();
    
        sq_distance_from_initial_pt = (test_x-test_x0)*(test_x-test_x0) + (test_y-test_y0)*(test_y-test_y0) + (test_z-test_z0)*(test_z-test_z0);
        if (!volume_sampling || (sq_distance_from_initial_pt < .25 * diameter * diameter))
        {
            makeVerletList();
            expandTestParticle();
            if (diameter > min_diameter) 
            {

                // correct for box edges...
                while (test_x >= box_x) test_x -= box_x;
                while (test_x < 0) test_x += box_x;
                while (test_y >= box_y) test_y -= box_y;
                while (test_y < 0) test_y += box_y;
                while (test_z >= box_z) test_z -= box_z;
                while (test_z < 0) test_z += box_z;

                result.pushBack(Cavity(test_x, test_y, test_z, diameter));
            }
        }
    }
  
} // end DDX::execute()


void DDX::generateTestPoint()
{
    test_x = test_x0 = rnd() * box_x;
    test_y = test_y0 = rnd() * box_y;
    test_z = test_z0 = rnd() * box_z;

    makeVerletList();
} // end DDX::generateTestPoint()


void DDX::makeVerletList()
{
    int i;
    vacuumms_double dx, dy, dz, dd;
    vacuumms_double shift_x, shift_y, shift_z;

    while (test_x > box_x) test_x -= box_x;
    while (test_y > box_y) test_y -= box_y;
    while (test_z > box_z) test_z -= box_z;

    while (test_x < 0) test_x += box_x;
    while (test_y < 0) test_y += box_y;
    while (test_z < 0) test_z += box_z;

    verlet_center_x=test_x;
    verlet_center_y=test_y;
    verlet_center_z=test_z;

    close_molecules=0;
    for (i=0; i<number_of_molecules; i++)
    {
        for (int index_x = -verlet_extent; index_x <= verlet_extent; index_x++)
        for (int index_y = -verlet_extent; index_y <= verlet_extent; index_y++)
        for (int index_z = -verlet_extent; index_z <= verlet_extent; index_z++)
        {
            shift_x = index_x * box_x;
            shift_y = index_y * box_y;
            shift_z = index_z * box_z;

            dx = shift_x + x[i] - test_x;
            dy = shift_y + y[i] - test_y;
            dz = shift_z + z[i] - test_z;

            dd = dx*dx + dy*dy + dz*dz;

            if (dd < verlet_cutoff) 
            {  
                close_x[close_molecules] = shift_x + x[i];
                close_y[close_molecules] = shift_y + y[i];
                close_z[close_molecules] = shift_z + z[i];
                close_sigma[close_molecules] = sigma[i];
                close_sigma6[close_molecules] = sigma[i]*sigma[i]*sigma[i]*sigma[i]*sigma[i]*sigma[i];
                close_sigma12[close_molecules] = close_sigma6[close_molecules]*close_sigma6[close_molecules];
                close_epsilon[close_molecules] = epsilon[i];

                close_molecules++;
            }
/*
        for (shift_x = -box_x; shift_x <= box_x; shift_x += box_x)
        for (shift_y = -box_y; shift_y <= box_y; shift_y += box_y)
        for (shift_z = -box_z; shift_z <= box_z; shift_z += box_z)
        {
            dx = shift_x + x[i] - test_x;
            dy = shift_y + y[i] - test_y;
            dz = shift_z + z[i] - test_z;

            dd = dx*dx + dy*dy + dz*dz;

            if (dd < verlet_cutoff) 
            {  
                close_x[close_molecules] = shift_x + x[i];
                close_y[close_molecules] = shift_y + y[i];
                close_z[close_molecules] = shift_z + z[i];
                close_sigma[close_molecules] = sigma[i];
                close_sigma6[close_molecules] = sigma[i]*sigma[i]*sigma[i]*sigma[i]*sigma[i]*sigma[i];
                close_sigma12[close_molecules] = close_sigma6[close_molecules]*close_sigma6[close_molecules];
                close_epsilon[close_molecules] = epsilon[i];

                close_molecules++;
            }
*/
        }
    }
} // end DDX::makeVerletList()


void DDX::findEnergyMinimum()
{
    vacuumms_double dx, dy, dz, dd, d6, d14;
    vacuumms_double factor;
    vacuumms_double old_energy;
    vacuumms_double new_energy;
    vacuumms_double grad_x, grad_y, grad_z;
//    vacuumms_double step_x, step_y, step_z;
    int i;
    vacuumms_double drift_sq;

/*
    vacuumms_float step_x = characteristic_length;
    vacuumms_float step_y = characteristic_length;
    vacuumms_float step_z = characteristic_length;
*/

    vacuumms_float alpha = 0.95f;
    vacuumms_float multiplier = 0.001f;
    vacuumms_float learning_rate = 0.01f;

    makeVerletList();

    // begin loop to iterate until minimum found
    for (attempts=0; attempts<n_steps; attempts++)
    {
        drift_sq = (test_x-verlet_center_x)*(test_x-verlet_center_x) 
                 + (test_y-verlet_center_y)*(test_y-verlet_center_y) 
                 + (test_z-verlet_center_z)*(test_z-verlet_center_z);

        if (drift_sq > .01 * verlet_cutoff) makeVerletList();

        // find the gradient at test_x, test_y, test_Z using the derivative of energy
        grad_x=0; grad_y=0; grad_z=0;

        for (i=0; i<close_molecules; i++)
        {
            dx = test_x - close_x[i];
            dy = test_y - close_y[i];
            dz = test_z - close_z[i];
            dd = dx*dx + dy*dy + dz*dz;
            d6 = dd*dd*dd;
            d14 = d6*d6*dd;

            // The analytical expression for the gradient contribution contains a factor of -48.0.
            // The minus is reflected in the sense of the step taken.  The factor of 48 is factored out in the normalization.
            factor = close_epsilon[i] * close_sigma12[i] / d14;

            grad_x += dx * factor;
            grad_y += dy * factor;
            grad_z += dz * factor;
        }

        // normalize the gradient
        vacuumms_double grad_sq = grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
        vacuumms_double grad_modulus = sqrt(grad_sq);
//std::cout << "grad_modulus: " << grad_modulus << std::endl;
//FTW declare convergence if grad_modulus is small
if (grad_modulus < 10.0f) 
{
//    std::cout << "grad_modulus indicates convergence, stopping after " << attempts << std::endl;
    break;
}

        grad_x /= grad_modulus;
        grad_y /= grad_modulus;
        grad_z /= grad_modulus;

        old_energy = calculateRepulsion();

/*
        vacuumms_double alpha = 0.5;
        step_x = grad_x * characteristic_energy * characteristic_length; while (step_x * step_x > characteristic_length * characteristic_length * precision_parameter * precision_parameter) {step_x *= alpha;}
        step_y = grad_y * characteristic_energy * characteristic_length; while (step_y * step_y > characteristic_length * characteristic_length * precision_parameter * precision_parameter) {step_y *= alpha;}
        step_z = grad_z * characteristic_energy * characteristic_length; while (step_z * step_z > characteristic_length * characteristic_length * precision_parameter * precision_parameter) {step_z *= alpha;}
*/

        // multiplier *= alpha;
        vacuumms_float step_x = learning_rate * grad_x; // * characteristic_length; 
        vacuumms_float step_y = learning_rate * grad_y; // * characteristic_length; 
        vacuumms_float step_z = learning_rate * grad_z; // * characteristic_length; 

//std::cout << attempts << ": " << step_x << ", " << step_y << ", " << step_z << std::endl;

        // removed this criteria for assessing minima.... step size no longer shrinks because gradient is now normalized
        // step_sq = step_x * step_x + step_y * step_y + step_z * step_z;
        // if (step_sq < (characteristic_length * characteristic_length * precision_parameter * precision_parameter)) break;  // close enough, exit loop

        test_x += step_x;
        test_y += step_y;
        test_z += step_z;
 
/*
        // check repulsion at new location
        new_energy = calculateRepulsion();
        if (new_energy - old_energy > 0.0f) 
        {
std::cout << attempts << std::endl;
            break;
        }
*/

/*
        // if the energy fluctuates up by a fraction of the characteristic energy, call it
        if (new_energy - old_energy > precision_parameter * characteristic_energy)
        {
            break;
        }
*/
    } // attempts
//std::cout << "done ." << std::endl;

} // end DDX::findEnergyMinimum()


vacuumms_double DDX::calculateRepulsion()
{
    vacuumms_double repulsion=0;
    vacuumms_double dx, dy, dz, dd, d6, d12;
    int i;

    for (i=0; i<close_molecules; i++)
    {
        dx = close_x[i] - test_x;
        dy = close_y[i] - test_y;
        dz = close_z[i] - test_z;
        dd = dx*dx + dy*dy + dz*dz;
        d6 = dd*dd*dd;
        d12 = d6*d6;

        repulsion += close_epsilon[i] * close_sigma12[i] / d12;
    }
 
    return 4.0 * repulsion;
} // end DDX::calculateRepulsion()


vacuumms_double DDX::calculateEnergy(vacuumms_double test_diameter)
{
    vacuumms_double repulsion=0;
    vacuumms_double attraction=0;
    vacuumms_double dx, dy, dz, dd, d6, d12;
    vacuumms_double sigma, sigma6, sigma12;
    int i;

    for (i=0; i<close_molecules; i++)
    {
        dx = close_x[i] - test_x;
        dy = close_y[i] - test_y;
        dz = close_z[i] - test_z;
        dd = dx*dx + dy*dy + dz*dz;
        d6 = dd*dd*dd;
        d12 = d6*d6;

        sigma = 0.5 * (close_sigma[i] + test_diameter);
        sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
        sigma12 = sigma6*sigma6;

        repulsion += close_epsilon[i] * sigma12/d12;
        attraction += close_epsilon[i] * sigma6/d6;
    }

    return 4.0 * (repulsion - attraction);
} // end DDX::calculateEnergy()


void DDX::expandTestParticle()
{
    vacuumms_double slope;
    vacuumms_double step_size;
    vacuumms_double energy, old_energy;
    vacuumms_double e0, e1, r0, r1;

    diameter = 0;

    // improved initial guess
    old_energy = calculateEnergy(diameter);
    if (old_energy > 0) return;
    while (diameter += .1)
    {
        energy = calculateEnergy(diameter);
        if (energy > old_energy) break;
        old_energy = energy;
    }
    
    // Newton's method

    while(1)
    {
        r0 = diameter - .001;
        r1 = diameter + .001;

        e0 = calculateEnergy(r0);
        e1 = calculateEnergy(r1);
        energy = calculateEnergy(diameter);

        slope = (e1-e0)/(r1-r0);
        step_size = -energy/slope;
//std::cout << "step_size" << step_size << std::endl;
        diameter = diameter + step_size;

//FTW        if (step_size*step_size < .00000001) break;
        if (step_size*step_size < 1.0e-24) break;
    }
} // end DDX::expandTestParticle()
