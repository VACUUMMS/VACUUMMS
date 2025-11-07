/* vacuumms_cpp/ddx.cc */

#include <vacuumms/ddx.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>

#include <vacuumms/limits.h>
#include <vacuumms/rng.hh>

#include <math.h>


DDX::DDX(Configuration _configuration, Parameters _parameters) //: 
//    configuration{_configuration}, parameters{_parameters} 
{
    setConfiguration(_configuration);
    setParameters(_parameters);
}

DDX::DDX(Configuration _configuration) //: 
//    configuration{_configuration} 
{
    setConfiguration(_configuration);
}

DDX::DDX()
{
}


void DDX::setParameters(Parameters _parameters)
{
    parameters = _parameters;

    parameters.getFloatParam((char*)"-verlet_cutoff", &verlet_cutoff);
    parameters.getIntParam((char*)"-verlet_extent", &verlet_extent);
    parameters.getIntParam((char*)"-number_of_samples", &number_of_samples);
    parameters.getIntParam((char*)"-number_of_steps", &number_of_steps);
    parameters.getIntParam((char*)"-rng_seed", &rng_seed);
    parameters.getFloatParam((char*)"-min_diameter", &min_diameter);
    parameters.getFloatParam((char*)"-learning_rate", &learning_rate);
    parameters.getFloatParam((char*)"-tolerance", &tolerance);
    volume_sampling = parameters.getFlagParam((char*)"-volume_sampling");
}


void DDX::setConfiguration(Configuration _configuration)
{
    configuration = _configuration;
}


Configuration DDX::getConfiguration()
{
    return configuration;
}


void DDX::setNumberOfSamples(int _number_of_samples)
{
    number_of_samples = _number_of_samples;
}


void DDX::setRNGSeed(int _rng_seed)
{
    rng_seed = _rng_seed;
    rng = MersenneTwister(rng_seed);
}


void DDX::setVerletCutoff(vacuumms_float _verlet_cutoff)
{
    verlet_cutoff = _verlet_cutoff;
}


void DDX::setVerletExtent(int _verlet_extent)
{
    verlet_extent = _verlet_extent;
}


void DDX::setNumberOfSteps(int _number_of_steps)
{
    number_of_steps = _number_of_steps;
}


void DDX::setMinDiameter(vacuumms_float _min_diameter)
{
    min_diameter = _min_diameter;
}


void DDX::setLearningRate(vacuumms_float _learning_rate)
{
    learning_rate = _learning_rate;
}


void DDX::setTolerance(vacuumms_float _tolerance)
{
    tolerance = _tolerance;
}


#ifdef BUILD_PYBIND_BINDINGS
pybind11::str DDX::__repr__()
{
    pybind11::str retval;
    retval += configuration.__repr__();
    retval += parameters.__repr__();
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
    printf("\t\t-number_of_steps [ 1000 ] (roughly reciprocal of precision parameter)\n");
    printf("\t\t-show_steps (includes steps taken as final column)\n");
    printf("\t\t-verlet_cutoff [ 100.0 ]\n");
    printf("\t\t-n [ 1 ]\n");
    printf("\t\t-volume_sampling \n");
    printf("\t\t-min_diameter [ 0.0 ]");
    printf("\n");
}


void DDX::execute()
{
    // Clear old results, if any
    result.reset();
    result.setBoxDimensions(configuration.getBoxDimensions());

    if (configuration.getSize() == 0) 
    {
        std::cout << "no atoms in configuration, declining to execute." << std::endl;
        return;
    }
    
    // Get box_dims from config info
    std::vector<vacuumms_float> box_dims = configuration.getBoxDimensions();
    box_x = box_dims[0];
    box_y = box_dims[1];
    box_z = box_dims[2];

    if (box_x * box_y * box_z < 0.000001) 
    {
        std::cout << "vanishingly small box volume set in configuration, declining to execute." << std::endl;
        return;
    }

    vacuumms_float sq_distance_from_initial_pt;

    if ((box_x * box_y * box_z) == 0.0) 
    {
        fprintf(stderr, "Found simulation box volume = 0.0, gracefully exiting.\n");
        return;
    }
  
    number_of_molecules = configuration.getSize();
  
    for (int sample_number = 0; sample_number < number_of_samples; sample_number++)
    {
        generateTestPoint();
        while (calculateEnergy(0.0) > 0.0f) generateTestPoint();
    
        findEnergyMinimum();

        makeVerletList();
        expandTestParticle();
    
        sq_distance_from_initial_pt = (test_x-test_x0)*(test_x-test_x0) + (test_y-test_y0)*(test_y-test_y0) + (test_z-test_z0)*(test_z-test_z0);
        if (!volume_sampling || (sq_distance_from_initial_pt < .25 * diameter * diameter))
        {
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
    test_x = test_x0 = rng.next_float() * box_x;
    test_y = test_y0 = rng.next_float() * box_y;
    test_z = test_z0 = rng.next_float() * box_z;

    makeVerletList();
} // end DDX::generateTestPoint()


void DDX::makeVerletList()
{
    verlet_list.clear(); // clear old list

    int i;
    vacuumms_float dx, dy, dz, dd;
    vacuumms_float shift_x, shift_y, shift_z;

    while (test_x > box_x) test_x -= box_x;
    while (test_y > box_y) test_y -= box_y;
    while (test_z > box_z) test_z -= box_z;

    while (test_x < 0) test_x += box_x;
    while (test_y < 0) test_y += box_y;
    while (test_z < 0) test_z += box_z;

    verlet_center_x=test_x;
    verlet_center_y=test_y;
    verlet_center_z=test_z;

    for (i=0; i < configuration.getSize(); i++)
    {
        for (int index_x = -verlet_extent; index_x <= verlet_extent; index_x++)
        for (int index_y = -verlet_extent; index_y <= verlet_extent; index_y++)
        for (int index_z = -verlet_extent; index_z <= verlet_extent; index_z++)
        {
            shift_x = index_x * box_x;
            shift_y = index_y * box_y;
            shift_z = index_z * box_z;

            dx = shift_x + configuration.recordAt(i).x - test_x;
            dy = shift_y + configuration.recordAt(i).y - test_y;
            dz = shift_z + configuration.recordAt(i).z - test_z;

            dd = dx*dx + dy*dy + dz*dz;

            if (dd < verlet_cutoff) 
            {  
                vacuumms_float close_x = shift_x + configuration.recordAt(i).x;
                vacuumms_float close_y = shift_y + configuration.recordAt(i).y;
                vacuumms_float close_z = shift_z + configuration.recordAt(i).z;

                vacuumms_float close_sigma = configuration.recordAt(i).sigma;
                vacuumms_float close_epsilon = configuration.recordAt(i).epsilon;
                ConfigurationRecord close_atom(close_x, close_y, close_z, close_sigma, close_epsilon);
                verlet_list.pushBack(close_atom);
            }
        }
    }
} // end DDX::makeVerletList()


void DDX::findEnergyMinimum()
{
    vacuumms_float dx, dy, dz, dd, d6, d14;
    vacuumms_float factor;
    vacuumms_float old_energy;
    vacuumms_float new_energy;
    vacuumms_float grad_x, grad_y, grad_z;
    int i;
    vacuumms_float drift_sq;

    vacuumms_float alpha = 0.95f;
    vacuumms_float multiplier = 0.001f;
    vacuumms_float learning_rate = 0.01f;

    makeVerletList();

    // begin loop to iterate until minimum found
    for (int attempts=0; attempts<number_of_steps; attempts++)
    {
        drift_sq = (test_x-verlet_center_x)*(test_x-verlet_center_x) 
                 + (test_y-verlet_center_y)*(test_y-verlet_center_y) 
                 + (test_z-verlet_center_z)*(test_z-verlet_center_z);

        if (drift_sq > .01 * verlet_cutoff) makeVerletList();

        // find the gradient at test_x, test_y, test_Z using the derivative of energy
        grad_x=0; grad_y=0; grad_z=0;

        for (int i = 0; i < verlet_list.getSize(); i++)
        {
            dx = test_x - verlet_list.recordAt(i).x;
            dy = test_y - verlet_list.recordAt(i).y;
            dz = test_z - verlet_list.recordAt(i).z;
            dd = dx*dx + dy*dy + dz*dz;
            d6 = dd*dd*dd;
            d14 = d6*d6*dd;

            // The analytical expression for the gradient contribution contains a factor of -48.0.
            // The minus is reflected in the sense of the step taken.  The factor of 48 is factored out in the normalization.
            factor = verlet_list.recordAt(i).epsilon * verlet_list.recordAt(i).sigma / d14;

            grad_x += dx * factor;
            grad_y += dy * factor;
            grad_z += dz * factor;

        }

        // normalize the gradient
        vacuumms_float grad_sq = grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
        vacuumms_float grad_modulus = sqrt(grad_sq);

        // declare convergence if grad_modulus is small
        if (grad_modulus < 10.0f) 
        {
            break;
        }

        grad_x /= grad_modulus;
        grad_y /= grad_modulus;
        grad_z /= grad_modulus;

        old_energy = calculateRepulsion();

        vacuumms_float step_x = learning_rate * grad_x;
        vacuumms_float step_y = learning_rate * grad_y;
        vacuumms_float step_z = learning_rate * grad_z;

        test_x += step_x;
        test_y += step_y;
        test_z += step_z;
 
    } // attempts

} // end DDX::findEnergyMinimum()


vacuumms_float DDX::calculateRepulsion()
{
    vacuumms_float repulsion=0;
    vacuumms_float dx, dy, dz, dd, d6, d12;
    int i;

    for (i=0; i<verlet_list.getSize(); i++)
    {
        dx = verlet_list.recordAt(i).x - test_x;
        dy = verlet_list.recordAt(i).y - test_y;
        dz = verlet_list.recordAt(i).z - test_z;
        dd = dx*dx + dy*dy + dz*dz;
        d6 = dd*dd*dd;
        d12 = d6*d6;

        vacuumms_float sigma = verlet_list.recordAt(i).sigma;
        vacuumms_float sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
        vacuumms_float sigma12 = sigma6 * sigma6;
        repulsion += verlet_list.recordAt(i).epsilon * sigma / d12;
    }
 
    return 4.0 * repulsion;
} // end DDX::calculateRepulsion()


vacuumms_float DDX::calculateEnergy(vacuumms_float test_diameter)
{
    vacuumms_float repulsion=0;
    vacuumms_float attraction=0;
    vacuumms_float dx, dy, dz, dd, d6, d12;
    vacuumms_float sigma, sigma6, sigma12;
    int i;

    for (i=0; i<verlet_list.getSize(); i++)
    {
        dx = verlet_list.recordAt(i).x - test_x;
        dy = verlet_list.recordAt(i).y - test_y;
        dz = verlet_list.recordAt(i).z - test_z;
        dd = dx*dx + dy*dy + dz*dz;
        d6 = dd*dd*dd;
        d12 = d6*d6;

        sigma = 0.5 * (verlet_list.recordAt(i).sigma + test_diameter);
        sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
        sigma12 = sigma6*sigma6;

        repulsion += verlet_list.recordAt(i).epsilon * sigma12/d12;
        attraction += verlet_list.recordAt(i).epsilon * sigma6/d6;
    }

    vacuumms_float energy = 4.0 * (repulsion - attraction);
    return energy;
} // end DDX::calculateEnergy()


void DDX::expandTestParticle()
{
    vacuumms_float step_tolerance = 1.0e-6;
    vacuumms_float h = 1.0e-6; // Finite difference step size

    // Initial guess
    diameter = 0.0f;
    vacuumms_float diameter_step = 0.01;
    while(calculateEnergy(diameter += diameter_step) < 0);
    
    //while (iteration++ < number_of_steps) 
    for (int iteration = 0; iteration < number_of_steps; iteration++) 
    {
        vacuumms_float energy = calculateEnergy(diameter);
std::cout << "energy= " << energy << std::endl;
        vacuumms_float d_energy = (calculateEnergy(diameter + h) - calculateEnergy(diameter - h)) / (2.0 * h);

        if (fabs(d_energy) < 1e-10)
        {
            printf("Error: Derivative too small: %f\n", d_energy);
            fflush(stdout);
            return;
        }

        vacuumms_float step_size = - energy / d_energy;

        if ((fabs(step_size) < step_tolerance) || (fabs(energy) < step_tolerance)) 
        {
            diameter += step_size;
            return;
        }

        diameter += step_size;
    }

    // ran out of iterations, return diameter without explicit convergence
    return;
}

