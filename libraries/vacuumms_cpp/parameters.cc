#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <string>
#include <iostream>
#include <vector>

#ifdef BUILD_BOOST_PYTHON_BINDINGS
#include <boost/python.hpp>
#endif

#include <vacuumms/parameters.hh>
#include <vacuumms/types.h>


Parameters::Parameters(int argc, char **argv)
{
    parameter_argc = argc;
    for (int i=0; i<argc; i++) parameter_argv.push_back(argv[i]);
}


Parameters::Parameters()
{
    parameter_argc = 0;
}

Parameters::Parameters(std::vector<std::string> p) : parameter_argv{p}
{
    parameter_argc = parameter_argc = parameter_argv.size();
}

#ifdef BUILD_BOOST_PYTHON_BINDINGS

Parameters::Parameters(const boost::python::list& _argv)
{
    parameter_argc = boost::python::len(_argv);

    for (int i=0; i<parameter_argc; i++)
    {
        char *arg = boost::python::extract<char*>(_argv[i]);
        parameter_argv.push_back(arg);
    }
}


const char* Parameters::getStringParam(char *param_name)
{
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
	if (i+1>=parameter_argc) 
	{
	    printf("reached EOL with no value specified for %s\n", param_name);
	    exit(1);
	}
	return parameter_argv[++i].c_str();
    }
    return NULL;
}


int Parameters::getIntParam(char *param_name)
{
    int retval = 0;

    int parameter, *p_parameter;
    p_parameter = &parameter;
    *p_parameter = 0;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name)
    {
	if (i+1>=parameter_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*p_parameter = (strtol(parameter_argv[++i].c_str(), NULL, 10));
	retval = 1;
    }
    return parameter;
}

long Parameters::getLongParam(char* param_name)
{
    int retval = 0;
 
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
	if (i+1>=parameter_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	retval = (strtol(parameter_argv[++i].c_str(), NULL, 10));
    }
    return retval;
}

vacuumms_float Parameters::getFloatParam(char *param_name)
{
    vacuumms_float retval = -0.0f;
 
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name)
    {
	if (i+1>=parameter_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	retval = (vacuumms_float)strtod(parameter_argv[++i].c_str(), NULL);
    }
    return retval;
}

double Parameters::getDoubleParam(char *param_name)
{
    for(int i=0; i<parameter_argc; i++)
        if (parameter_argv[i] == param_name)
        {
            if (i+1>=parameter_argc) 
            {
	        printf("no value specified for %s\n", param_name);
                exit(1);
            }
            return strtod(parameter_argv[++i].c_str(), NULL);
        }
    return -0.0f;
}

boost::python::list Parameters::getVectorParam(char *param_name)
{
    //std::vector<double> retval;
    boost::python::list retval;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
	if (i+3>=parameter_argc) 
	{
	    printf("not enough values specified for %s\n", param_name);
	    exit(1);
	}

        retval.append(strtod(parameter_argv[++i].c_str(), NULL));
        retval.append(strtod(parameter_argv[++i].c_str(), NULL));
        retval.append(strtod(parameter_argv[++i].c_str(), NULL));
    }
    return retval;
}

boost::python::list Parameters::getVectorStringParam(char *param_name)
{
    //std::vector<std::string> retval;
    boost::python::list retval;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
	if (i+3>=parameter_argc) 
	{
	    printf("not enough values specified for %s\n", param_name);
	    exit(1);
	}

        retval.append(parameter_argv[++i].c_str());
        retval.append(parameter_argv[++i].c_str());
        retval.append(parameter_argv[++i].c_str());
    }
    return retval;
}

#endif

int Parameters::getIntParam(char *param_name)
{
    int parameter, *p_parameter;
    p_parameter = &parameter;
    *p_parameter = -1;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name)
    {
        if (i+1>=parameter_argc) 
        {
            printf("no value specified for %s\n", param_name);
            break;
        }

        *p_parameter = (strtol(parameter_argv[++i].c_str(), NULL, 10));
    }
    return parameter;
}

vacuumms_float Parameters::getFloatParam(char *param_name)
{
    vacuumms_float parameter, *p_parameter;
    p_parameter = &parameter;
    *p_parameter = -0.0;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name)
    {
        if (i+1>=parameter_argc) 
        {
            printf("no value specified for %s\n", param_name);
            break;
        }

        *p_parameter = (strtod(parameter_argv[++i].c_str(), NULL));
    }
    return parameter;
}

const char* Parameters::getStringParam(char* param_name)
{
    const char* parameter;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name)
    {
        if (i+1>=parameter_argc) 
        {
            printf("no value specified for %s\n", param_name);
            break;
        }

        parameter = parameter_argv[++i].c_str();
    }
    return parameter;
}


#ifdef BUILD_PYBIND_BINDINGS

pybind11::str Parameters::__repr__()
{
    pybind11::str retval("");

    for (int i=0; i<parameter_argc; i++)
        retval = retval + pybind11::str(parameter_argv[i]) + pybind11::str("\n");
    
    return retval;
}

pybind11::str Parameters::__str__()
{
    pybind11::str retval("");

    for (int i=0; i<parameter_argc; i++)
        retval = retval + pybind11::str(parameter_argv[i]) + pybind11::str("\n");
    return retval;
}

Parameters::Parameters(const pybind11::list& _argv)
{   
    parameter_argc = pybind11::len(_argv);

    for (const auto& item : _argv) 
    {
        parameter_argv.push_back(item.cast<std::string>().c_str());
    } 
}

/* FTW can i replace with std::vector type? yes
pybind11::list Parameters::getVectorParam(char *param_name)
{
    pybind11::list retval;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
        if (i+3>=parameter_argc) 
        {
            printf("not enough values specified for %s\n", param_name);
            break;
        }

        retval.append(strtod(parameter_argv[++i].c_str(), NULL));
        retval.append(strtod(parameter_argv[++i].c_str(), NULL));
        retval.append(strtod(parameter_argv[++i].c_str(), NULL));
    }
    return retval;
}
*/

pybind11::list Parameters::getVectorStringParam(char* param_name)
{
    pybind11::list retval;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
        if (i+3>=parameter_argc) 
        {
            printf("not enough values specified for %s\n", param_name);
            break;
        }

        retval.append(parameter_argv[++i].c_str());
        retval.append(parameter_argv[++i].c_str());
        retval.append(parameter_argv[++i].c_str());
    }
    return retval;
}


#endif // BUILD_PYBIND_BINDINGS


int Parameters::addParameter(const char* parameter)
{
    parameter_argv.push_back(parameter);
    parameter_argc = parameter_argv.size();
    return parameter_argc;
}


int Parameters::getStringParam(char *param_name, const char **parameter)
{
    int retval = 0;
 
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
	if (i+1>=parameter_argc) 
	{
	    printf("reached EOL with no value specified for %s\n", param_name);
	    exit(1);
	}
	// parameter = &parameter_argv[++i];
	//FTW*parameter = parameter_argv[++i];
	*parameter = parameter_argv[++i].c_str();
	retval = 1;
    }
    return retval;
}


int Parameters::getIntParam(char *param_name, int *p_parameter)
{
    int retval = 0;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name)
    {
	if (i+1>=parameter_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*p_parameter = (strtol(parameter_argv[++i].c_str(), NULL, 10));
	retval = 1;
    }
    return retval;
}


int Parameters::getLongParam(char *param_name, long *parameter)
{
    int retval = 0;
 
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
	if (i+1>=parameter_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (strtol(parameter_argv[++i].c_str(), NULL, 10));
	retval = 1;
    }
    return retval;
}


int Parameters::getFloatParam(char *param_name, vacuumms_float *parameter)
{
    int retval = 0;
 
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name)
    {
	if (i+1>=parameter_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (vacuumms_float)strtod(parameter_argv[++i].c_str(), NULL);
	retval = 1;
    }
    return retval;
}


int Parameters::getDoubleParam(char *param_name, double *parameter)
{
    int retval = 0;
 
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
	if (i+1>=parameter_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (strtod(parameter_argv[++i].c_str(), NULL));
	retval = 1;
    }
    return retval;
}


std::vector<vacuumms_float> Parameters::getVectorParam(const char* param_name)
{
    return Parameters::getVectorParam(std::string(param_name));
}


std::vector<vacuumms_float> Parameters::getVectorParam(char* param_name)
{
    return Parameters::getVectorParam(std::string(param_name));
}


std::vector<vacuumms_float> Parameters::getVectorParam(std::string param_name)
{
    std::vector<vacuumms_float> retval(3);

    int found = 0;

    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
        found = 1;
        if (i+3>=parameter_argc) 
        {
            fprintf(stderr, "not enough values specified for %s\n", param_name);
            return retval;
        }

        retval[0] = (strtod(parameter_argv[++i].c_str(), NULL));
        retval[1] = (strtod(parameter_argv[++i].c_str(), NULL));
        retval[2] = (strtod(parameter_argv[++i].c_str(), NULL));
    }
    if (!found) fprintf(stderr, "getVectorParam <<<%s>>> not found.\n", param_name.c_str());
    return retval;
}


int Parameters::getVectorParam(char *param_name, 
                               vacuumms_float *parameter1,  
                               vacuumms_float *parameter2, 
                               vacuumms_float *parameter3)
{
    int retval = 0;
 
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
	if (i+3>=parameter_argc) 
	{
	    printf("not enough values specified for %s\n", param_name);
	    exit(1);
	}

	*parameter1 = (strtod(parameter_argv[++i].c_str(), NULL));
	*parameter2 = (strtod(parameter_argv[++i].c_str(), NULL));
	*parameter3 = (strtod(parameter_argv[++i].c_str(), NULL));
	retval = 1;
    }
    return retval;
}


int Parameters::getVectorParam(char *param_name, 
                               double *parameter1,  
                               double *parameter2, 
                               double *parameter3)
{
    int retval = 0;
 
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name) 
    {
	if (i+3>=parameter_argc) 
	{
	    printf("not enough values specified for %s\n", param_name);
	    exit(1);
	}

	*parameter1 = (strtod(parameter_argv[++i].c_str(), NULL));
	*parameter2 = (strtod(parameter_argv[++i].c_str(), NULL));
	*parameter3 = (strtod(parameter_argv[++i].c_str(), NULL));
	retval = 1;
    }
    return retval;
}


int Parameters::getVectorStringParam(char *param_name, 
                                     const char **parameter1,  
                                     const char **parameter2, 
                                     const char **parameter3)
{
    int retval = 0;
 
    for (int i=0; i<parameter_argc; i++)
    if (parameter_argv[i] == param_name)
    {
        if (i+3>=parameter_argc) 
        {
            printf("not enough values specified for %s\n", param_name);
            exit(1);
        }

        *parameter1 = parameter_argv[++i].c_str();
        *parameter2 = parameter_argv[++i].c_str();
        *parameter3 = parameter_argv[++i].c_str();
        retval = 1;
    }
    return retval;
}


int Parameters::getFlagParam(char *param_name)
{
    for (int i=0; i<parameter_argc; i++) if (parameter_argv[i] == param_name) return 1;

    return 0;
}

