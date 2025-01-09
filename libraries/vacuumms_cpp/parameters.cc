#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <boost/python.hpp>

#include <string>
#include <iostream>

#include <vacuumms/parameters.hh>

namespace bp = boost::python;

Parameters::Parameters(int argc, char **argv)
{
    command_line_argc = argc;
    command_line_argv = argv;
}

Parameters::Parameters()
{

    command_line_argc = 0;
    command_line_argv = NULL;

    std::cout << "Constructing default\n";
}

Parameters::Parameters(const boost::python::list& _argv)
{

//    command_line_argc = 0;
//    command_line_argv = NULL;

    int _argc = bp::len(_argv);

    char **argv = (char**)malloc(sizeof(char*) * _argc);

    for (int i=0; i<_argc; i++)
    {
        char *arg = bp::extract<char*>(_argv[i]);
        argv[i] = arg;
        std::cout << argv[i] << std::endl;
    }

    free(argv);
}

Parameters::Parameters(std::string arg_list)
{

    command_line_argc = 0;
    command_line_argv = NULL;

    std::cout << "Constructing from arglist <" << arg_list << ">\n";
}

void Parameters::setParameters(std::string arg_list)
{

    command_line_argc = 0;
    command_line_argv = NULL;

    std::cout << "setting arglist <" << arg_list << ">\n";
}

int Parameters::getStringParam(char *param_name, char **parameter)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
    if (strcmp(command_line_argv[i], param_name) == 0) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("reached EOL with no value specified for %s\n", param_name);
	    exit(1);
	}
	// parameter = &command_line_argv[++i];
	*parameter = command_line_argv[++i];
	retval = 1;
    }
    return retval;
}

int Parameters::getIntParam(char *param_name, int *parameter)
{
    int i=0;
    int retval = 0;

    while (++i<command_line_argc)
    if (strcmp(command_line_argv[i], param_name) == 0) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (strtol(command_line_argv[++i], NULL, 10));
	retval = 1;
    }
    return retval;
}

int Parameters::getLongParam(char *param_name, long *parameter)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
    if (strcmp(command_line_argv[i], param_name) == 0) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (strtol(command_line_argv[++i], NULL, 10));
	retval = 1;
    }
    return retval;
}

int Parameters::getFloatParam(char *param_name, float *parameter)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
    if (strcmp(command_line_argv[i], param_name) == 0) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (float)strtod(command_line_argv[++i], NULL);
	retval = 1;
    }
    return retval;
}

int Parameters::getDoubleParam(char *param_name, double *parameter)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
    if (strcmp(command_line_argv[i], param_name) == 0) 
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (strtod(command_line_argv[++i], NULL));
	retval = 1;
    }
    return retval;
}

int Parameters::getVectorParam(char *param_name, double *parameter1,  double *parameter2, double *parameter3)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
    if (strcmp(command_line_argv[i], param_name) == 0) 
    {
	if (i+3>=command_line_argc) 
	{
	    printf("not enough values specified for %s\n", param_name);
	    exit(1);
	}

	*parameter1 = (strtod(command_line_argv[++i], NULL));
	*parameter2 = (strtod(command_line_argv[++i], NULL));
	*parameter3 = (strtod(command_line_argv[++i], NULL));
	retval = 1;
    }
    return retval;
}

int Parameters::getVectorStringParam(char *param_name, char **parameter1,  char **parameter2, char **parameter3)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
	if (strcmp(command_line_argv[i], param_name) == 0) 
    {
	if (i+3>=command_line_argc) 
	{
	    printf("not enough values specified for %s\n", param_name);
	    exit(1);
	}

	*parameter1 = command_line_argv[++i];
	*parameter2 = command_line_argv[++i];
	*parameter3 = command_line_argv[++i];
	retval = 1;
    }
    return retval;
}

int Parameters::getFlagParam(char *param_name)
{
    int i=0;

    while (++i<command_line_argc) if (strcmp(command_line_argv[i], param_name) == 0) return 1;

    return 0;
}
