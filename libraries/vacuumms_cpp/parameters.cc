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

//remember:        std::vector<std::string> command_line_argv;


Parameters::Parameters(int argc, char **argv)
{
    command_line_argc = argc;
    for (int i=0; i<argc; i++) command_line_argv.push_back(argv[i]);
}

Parameters::Parameters()
{
    command_line_argc = 0;
//    command_line_argv = std::vector<std::string> ();
}

#ifdef BUILD_BOOST_PYTHON_BINDINGS

Parameters::Parameters(const boost::python::list& _argv)
{
    command_line_argc = boost::python::len(_argv);

//    char **argv = (char**)malloc(sizeof(char*) * _argc);

    for (int i=0; i<command_line_argc; i++)
    {
        char *arg = boost::python::extract<char*>(_argv[i]);
        command_line_argv.push_back(arg);
std::cout << arg << std::endl;

    }

    for (int i=0; i<command_line_argc; i++) printf("%03d\t%s\n", i, command_line_argv[i].c_str());
}

#endif


//int Parameters::getStringParam(char *param_name, char **parameter)
int Parameters::getStringParam(char *param_name, const char **parameter)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
//    if (strcmp(command_line_argv[i], param_name) == 0) 
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

int Parameters::getIntParam(char *param_name, int *parameter)
{
    int i=0;
    int retval = 0;

    while (++i<command_line_argc)
//    if (strcmp(command_line_argv[i], param_name) == 0) 
    if (command_line_argv[i] == param_name)
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	//FTW*parameter = (strtol(command_line_argv[++i], NULL, 10));
	*parameter = (strtol(command_line_argv[++i].c_str(), NULL, 10));
	retval = 1;
    }
    return retval;
}

int Parameters::getLongParam(char *param_name, long *parameter)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
    //if (strcmp(command_line_argv[i], param_name) == 0) 
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

int Parameters::getFloatParam(char *param_name, float *parameter)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
    //FTWif (strcmp(command_line_argv[i], param_name) == 0) 
    if (command_line_argv[i] == param_name)
    {
	if (i+1>=command_line_argc) 
	{
	    printf("no value specified for %s\n", param_name);
	    exit(1);
	}

	*parameter = (float)strtod(command_line_argv[++i].c_str(), NULL);
	retval = 1;
    }
    return retval;
}

int Parameters::getDoubleParam(char *param_name, double *parameter)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
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

int Parameters::getVectorParam(char *param_name, double *parameter1,  double *parameter2, double *parameter3)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
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

int Parameters::getVectorStringParam(char *param_name, const char **parameter1,  const char **parameter2, const char **parameter3)
{
    int i=0;
    int retval = 0;
 
    while (++i<command_line_argc)
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
    int i=0;

    while (++i<command_line_argc) if (command_line_argv[i] == param_name) return 1;

    return 0;
}

