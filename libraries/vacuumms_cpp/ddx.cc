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


void DDX::setNumberOfThreads(int _number_of_threads)
{
    // Set equal to number of cores when zero is specified.
    if (_number_of_threads == 0) 
    {
        number_of_threads = std::thread::hardware_concurrency();
        std::cout << "Using hardware concurrency level = " << number_of_threads << std::endl;
    }
    else number_of_threads = _number_of_threads;
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


void DDX::reorderResults()
{
    CavityConfiguration ordered_results;
    for (int i = 0; i < results_order.size(); i++)
    for (int j = 0; j < results.getSize(); j++)
    if (results_order[i] == j) 
    {
        ordered_results.pushBack(results.recordAt(j));
        continue;
    }

    results = ordered_results;
}


#ifdef BUILD_PYBIND_BINDINGS
pybind11::str DDX::__repr__()
{
    pybind11::str retval;
    retval += configuration.__repr__();
    retval += parameters.__repr__();
    retval += results.__repr__();
    return retval;
}
#endif


CavityConfiguration DDX::getResult()
{
    return results;   
}


void DDX::printUsage()
{
    std::cout << 
        "Constructors: " << std::endl << std::endl <<
        "DDX(Configuration c);        # Construct from configuration" << std::endl <<
        "DDX();                       # construct empty DDX operation" << std::endl <<
        std::endl << 
        "Member functions:" << std::endl << 
        std::endl <<
        "DDX.setNumberOfThreads(1)    # Specify number of threads or 0 to use all available cores." << std::endl <<
        "DDX.setNumberOfSteps(50)     # Maximum number of steps before giving up and accepting result without explicit convergence" << std::endl <<
        "DDX.setLearningRate(0.01)    # Set learning rate for gradient descent of location of center." << std::endl <<
        "DDX.setTolerance(10.0)       # Set comparison value to derivative to determine convergence." << std::endl <<
        "DDX.setRNGSeed(1)            # Set random number generator seed." << std::endl <<
        "DDX.setVerletExtent(1)       # Set how many levels of mirror boxes to search when building Verlet list." << std::endl <<
        "DDX.setConfiguration();      # Set atom configuration." << std::endl <<
        "DDX.getConfiguration();      # Retrieve atom configuration." << std::endl <<
        "DDX.setNumberOfSamples(1);   # Specify number of samples to generate." << std::endl <<
        "DDX.setVerletCutoff(100.0);  # Specify square of radius to use when generating Verlet list." << std::endl <<
        "DDX.setMinDiameter(0.0);     # Specify minimum diameter of cavity to include." << std::endl <<
        "DDX.getResult();             # Retrieve results of computation/operation." << std::endl <<
        "DDX.execute();               # Execute the operation." << std::endl <<
        "DDX.printUsage();            # Generate this message." << std::endl <<
        std::endl;

/*
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
*/
}


void DDX::execute()
{
    // the new wrapper
    // Clear old results, if any
    results.reset();
    results.setBoxDimensions(configuration.getBoxDimensions());

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

    // Initialize semaphore
    if (sem_init(&semaphore, 0, number_of_threads) != 0) 
    {
        std::cout << "DDX::execute(): Could not initialize semaphore. " << std::endl;
        return;
    }

    // Now launch the workers...
    for (int sample_number = 0; sample_number < number_of_samples; sample_number++)
    {
        // Launch threads
        threads.emplace_back(&DDX::run_sample, this, sample_number);
    }

    // Wait for completion
    for (auto& t : threads) {
        if(t.joinable()) t.join();
    }

    // Cleanup
    sem_destroy(&semaphore);

} // end DDX::execute()


// Thread function to run single sample for id
void DDX::run_sample(int id)
{
    sem_wait(&semaphore);
    Sample s(this, id); 
    sem_post(&semaphore);
}


// Do everything from the constructor
DDX::Sample::Sample(DDX* _outer, int _id)
{
    outer = _outer;
    id = _id;

    // Initialize RNG based on seed and thread id
    rng = MersenneTwister(outer->rng_seed + id);

    generateTestPoint();
    while (calculateEnergy(0.0) > 0.0f) generateTestPoint();
    
    findEnergyMinimum();

    makeVerletList();
    expandTestParticle();

    // Discards point if outside of cavity when using volume sampling
    vacuumms_float sq_distance_from_initial_pt = 
            (test_x-test_x0) * (test_x-test_x0) 
          + (test_y-test_y0) * (test_y-test_y0) 
          + (test_z-test_z0) * (test_z-test_z0);
    if (!outer->volume_sampling || (sq_distance_from_initial_pt < .25 * diameter * diameter))
    {
        if (diameter > outer->min_diameter) 
        {
            // correct for box edges...
            while (test_x >= outer->box_x) test_x -= outer->box_x;
            while (test_x < 0) test_x += outer->box_x;
            while (test_y >= outer->box_y) test_y -= outer->box_y;
            while (test_y < 0) test_y += outer->box_y;
            while (test_z >= outer->box_z) test_z -= outer->box_z;
            while (test_z < 0) test_z += outer->box_z;

            // Synchronized code to write result
            {
                std::lock_guard<std::mutex> lock(outer->results_mutex);
                outer->results.pushBack(Cavity(test_x, test_y, test_z, diameter));
                outer->results_order.push_back(id);
            }
            
        }
    }

} // end Sample::Sample()


void DDX::Sample::generateTestPoint()
{
    test_x = test_x0 = rng.next_float() * outer->box_x;
    test_y = test_y0 = rng.next_float() * outer->box_y;
    test_z = test_z0 = rng.next_float() * outer->box_z;

    makeVerletList();
} // end DDX::generateTestPoint()


void DDX::Sample::makeVerletList()
{
    verlet_list.clear(); // clear old list

    int i;
    vacuumms_float dx, dy, dz, dd;
    vacuumms_float shift_x, shift_y, shift_z;

    while (test_x > outer->box_x) test_x -= outer->box_x;
    while (test_y > outer->box_y) test_y -= outer->box_y;
    while (test_z > outer->box_z) test_z -= outer->box_z;

    while (test_x < 0) test_x += outer->box_x;
    while (test_y < 0) test_y += outer->box_y;
    while (test_z < 0) test_z += outer->box_z;

    verlet_center_x=test_x;
    verlet_center_y=test_y;
    verlet_center_z=test_z;

    for (i=0; i < outer->configuration.getSize(); i++)
    {
        for (int index_x = -outer->verlet_extent; index_x <= outer->verlet_extent; index_x++)
        for (int index_y = -outer->verlet_extent; index_y <= outer->verlet_extent; index_y++)
        for (int index_z = -outer->verlet_extent; index_z <= outer->verlet_extent; index_z++)
        {
            shift_x = index_x * outer->box_x;
            shift_y = index_y * outer->box_y;
            shift_z = index_z * outer->box_z;

            dx = shift_x + outer->configuration.recordAt(i).x - test_x;
            dy = shift_y + outer->configuration.recordAt(i).y - test_y;
            dz = shift_z + outer->configuration.recordAt(i).z - test_z;

            dd = dx*dx + dy*dy + dz*dz;

            if (dd < outer->verlet_cutoff) 
            {  
                vacuumms_float close_x = shift_x + outer->configuration.recordAt(i).x;
                vacuumms_float close_y = shift_y + outer->configuration.recordAt(i).y;
                vacuumms_float close_z = shift_z + outer->configuration.recordAt(i).z;

                vacuumms_float close_sigma = outer->configuration.recordAt(i).sigma;
                vacuumms_float close_epsilon = outer->configuration.recordAt(i).epsilon;
                ConfigurationRecord close_atom(close_x, close_y, close_z, close_sigma, close_epsilon);
                verlet_list.pushBack(close_atom);
            }
        }
    }
} // end DDX::makeVerletList()


void DDX::Sample::findEnergyMinimum()
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
    for (int attempts=0; attempts<outer->number_of_steps; attempts++)
    {
        drift_sq = (test_x-verlet_center_x)*(test_x-verlet_center_x) 
                 + (test_y-verlet_center_y)*(test_y-verlet_center_y) 
                 + (test_z-verlet_center_z)*(test_z-verlet_center_z);

        if (drift_sq > .01 * outer->verlet_cutoff) makeVerletList();

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


vacuumms_float DDX::Sample::calculateRepulsion()
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


vacuumms_float DDX::Sample::calculateEnergy(vacuumms_float test_diameter)
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


void DDX::Sample::expandTestParticle()
{
    vacuumms_float step_tolerance = 1.0e-6;
    vacuumms_float h = 1.0e-6; // Finite difference step size

    // Initial guess
    diameter = 0.0f;
    vacuumms_float diameter_step = 0.01;
    while(calculateEnergy(diameter += diameter_step) < 0);
    
    //while (iteration++ < number_of_steps) 
    for (int iteration = 0; iteration < outer->number_of_steps; iteration++) 
    {
        vacuumms_float energy = calculateEnergy(diameter);
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

