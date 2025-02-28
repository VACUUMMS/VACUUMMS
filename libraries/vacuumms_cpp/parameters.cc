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
    command_line_argc = argc;
    for (int i=0; i<argc; i++) command_line_argv.push_back(argv[i]);
}


Parameters::Parameters()
{
    command_line_argc = 0;
}

#ifdef BUILD_BOOST_PYTHON_BINDINGS

Parameters::Parameters(const boost::python::list& _argv)
{
    command_line_argc = boost::python::len(_argv);

    for (int i=0; i<command_line_argc; i++)
    {
        char *arg = boost::python::extract<char*>(_argv[i]);
        command_line_argv.push_back(arg);
    }
}


const char* Parameters::getStringParam(char *param_name)
{
    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("reached EOL with no value specified for %s\n", param_name);
	    exit(1);
	}
	return command_line_argv[++i].c_str();
    }
    return NULL;
}


int Parameters::getIntParam(char *param_name)
{
    int retval = 0;

    int parameter, *p_parameter;
    p_parameter = &parameter;
    *p_parameter = 0;

    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name)
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*p_parameter = (strtol(command_line_argv[++i].c_str(), NULL, 10));
	retval = 1;
    }
    return parameter;
}

long Parameters::getLongParam(char* param_name)
{
    int retval = 0;
 
    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	retval = (strtol(command_line_argv[++i].c_str(), NULL, 10));
    }
    return retval;
}

vacuumms_float Parameters::getFloatParam(char *param_name)
{
    vacuumms_float retval = -0.0f;
 
    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name)
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	retval = (vacuumms_float)strtod(command_line_argv[++i].c_str(), NULL);
    }
    return retval;
}

double Parameters::getDoubleParam(char *param_name)
{
    for(int i=0; i<command_line_argc; i++)
        if (command_line_argv[i] == param_name)
        {
            if (i+1>=command_line_argc) 
            {
	        printf("no value specified for %s\n", param_name);
                exit(1);
            }
            return strtod(command_line_argv[++i].c_str(), NULL);
        }
    return -0.0f;
}

boost::python::list Parameters::getVectorParam(char *param_name)
{
    //std::vector<double> retval;
    boost::python::list retval;

    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
	if (i+3>=command_line_argc) 
	{
	    printf("not enough values specified for %s\n", param_name);
	    exit(1);
	}

        retval.append(strtod(command_line_argv[++i].c_str(), NULL));
        retval.append(strtod(command_line_argv[++i].c_str(), NULL));
        retval.append(strtod(command_line_argv[++i].c_str(), NULL));
    }
    return retval;
}

boost::python::list Parameters::getVectorStringParam(char *param_name)
{
    //std::vector<std::string> retval;
    boost::python::list retval;

    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
	if (i+3>=command_line_argc) 
	{
	    printf("not enough values specified for %s\n", param_name);
	    exit(1);
	}

        retval.append(command_line_argv[++i].c_str());
        retval.append(command_line_argv[++i].c_str());
        retval.append(command_line_argv[++i].c_str());
    }
    return retval;
}

#endif

int Parameters::getIntParam(char *param_name)
{
    int parameter, *p_parameter;
    p_parameter = &parameter;
    *p_parameter = -1;

    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name)
    {
        if (i+1>=command_line_argc) 
        {
            printf("no value specified for %s\n", param_name);
            break;
        }

        *p_parameter = (strtol(command_line_argv[++i].c_str(), NULL, 10));
    }
    return parameter;
}

vacuumms_float Parameters::getFloatParam(char *param_name)
{
    vacuumms_float parameter, *p_parameter;
    p_parameter = &parameter;
    *p_parameter = -0.0;

    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name)
    {
        if (i+1>=command_line_argc) 
        {
            printf("no value specified for %s\n", param_name);
            break;
        }

        *p_parameter = (strtod(command_line_argv[++i].c_str(), NULL));
    }
    return parameter;
}

const char* Parameters::getStringParam(char* param_name)
{
    const char* parameter;

    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name)
    {
        if (i+1>=command_line_argc) 
        {
            printf("no value specified for %s\n", param_name);
            break;
        }

        parameter = command_line_argv[++i].c_str();
    }
    return parameter;
}


#ifdef BUILD_PYBIND_BINDINGS

pybind11::str Parameters::__repr__()
{
    pybind11::str retval("");

    for (int i=0; i<command_line_argc; i++)
        retval = retval + pybind11::str(command_line_argv[i]) + pybind11::str("\n");
    
    return retval;
}

pybind11::str Parameters::__str__()
{
    pybind11::str retval("");

    for (int i=0; i<command_line_argc; i++)
        retval = retval + pybind11::str(command_line_argv[i]) + pybind11::str("\n");
    return retval;
}

Parameters::Parameters(const pybind11::list& _argv)
{   
    command_line_argc = pybind11::len(_argv);

    for (const auto& item : _argv) 
    {
        command_line_argv.push_back(item.cast<std::string>().c_str());
    } 
}

pybind11::list Parameters::getVectorParam(char *param_name)
{
    pybind11::list retval;

    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
        if (i+3>=command_line_argc) 
        {
            printf("not enough values specified for %s\n", param_name);
            break;
        }

        retval.append(strtod(command_line_argv[++i].c_str(), NULL));
        retval.append(strtod(command_line_argv[++i].c_str(), NULL));
        retval.append(strtod(command_line_argv[++i].c_str(), NULL));
    }
    return retval;
}

pybind11::list Parameters::getVectorStringParam(char* param_name)
{
    pybind11::list retval;

    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
        if (i+3>=command_line_argc) 
        {
            printf("not enough values specified for %s\n", param_name);
            break;
        }

        retval.append(command_line_argv[++i].c_str());
        retval.append(command_line_argv[++i].c_str());
        retval.append(command_line_argv[++i].c_str());
    }
    return retval;
}


#endif // BUILD_PYBIND_BINDINGS


int Parameters::addParameter(const char* parameter)
{
    command_line_argv.push_back(parameter);
    command_line_argc = command_line_argv.size();
    return command_line_argc;
}


int Parameters::getStringParam(char *param_name, const char **parameter)
{
    int retval = 0;
 
    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("reached EOL with no value specified for %s\n", param_name);
	    exit(1);
	}
	// parameter = &command_line_argv[++i];
	//FTW*parameter = command_line_argv[++i];
	*parameter = command_line_argv[++i].c_str();
	retval = 1;
    }
    return retval;
}


int Parameters::getIntParam(char *param_name, int *p_parameter)
{
    int retval = 0;

    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name)
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*p_parameter = (strtol(command_line_argv[++i].c_str(), NULL, 10));
	retval = 1;
    }
    return retval;
}


int Parameters::getLongParam(char *param_name, long *parameter)
{
    int retval = 0;
 
    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (strtol(command_line_argv[++i].c_str(), NULL, 10));
	retval = 1;
    }
    return retval;
}


int Parameters::getFloatParam(char *param_name, vacuumms_float *parameter)
{
    int retval = 0;
 
    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name)
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (vacuumms_float)strtod(command_line_argv[++i].c_str(), NULL);
	retval = 1;
    }
    return retval;
}


int Parameters::getDoubleParam(char *param_name, double *parameter)
{
    int retval = 0;
 
    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (strtod(command_line_argv[++i].c_str(), NULL));
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
 
    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name) 
    {
	if (i+3>=command_line_argc) 
	{
	    printf("not enough values specified for %s\n", param_name);
	    exit(1);
	}

	*parameter1 = (strtod(command_line_argv[++i].c_str(), NULL));
	*parameter2 = (strtod(command_line_argv[++i].c_str(), NULL));
	*parameter3 = (strtod(command_line_argv[++i].c_str(), NULL));
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
 
    for (int i=0; i<command_line_argc; i++)
    if (command_line_argv[i] == param_name)
    {
        if (i+3>=command_line_argc) 
        {
            printf("not enough values specified for %s\n", param_name);
            exit(1);
        }

        *parameter1 = command_line_argv[++i].c_str();
        *parameter2 = command_line_argv[++i].c_str();
        *parameter3 = command_line_argv[++i].c_str();
        retval = 1;
    }
    return retval;
}


int Parameters::getFlagParam(char *param_name)
{
    for (int i=0; i<command_line_argc; i++) if (command_line_argv[i] == param_name) return 1;

    return 0;
}

