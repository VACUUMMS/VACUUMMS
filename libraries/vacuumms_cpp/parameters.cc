#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

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


Parameters::Parameters(const char* filename)
{
    // open the file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file" << filename << "\n";
    }

    // Read the file line by line
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> values;
        std::stringstream ss(line);
        std::string value;
        if (!line.empty() && line[0] == '#') 
        {
            // std::cout << "ignoring comment line: " << line << "\n";
            continue;
        }

        // Split the line by whitespace
        while (ss >> value) 
	{
            parameter_argv.push_back(value);
        }

        // Process the values in the line
        for (const auto& val : values) {
            std::cout << val << " ";
        }
    }

    parameter_argc = parameter_argv.size();
    std::cout << "Read " << parameter_argc << " parameters." << "\n";

    // Close the file
    file.close();
}


int Parameters::toFile(const char* filename)
{
//    std::string filename(_filename);

    // Open the output file
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file" << filename << "\n";
        return 1;
    }

    // Write strings to the file
    for (size_t i = 0; i < parameter_argv.size(); ++i) {
        // Start a new line for strings starting with '-' (except for the first string)
        if (i > 0 && parameter_argv[i][0] == '-') {
            file << "\n";
        }
        // Write the string; add a space after it unless it's the last string
        file << parameter_argv[i];
        if (i < parameter_argv.size() - 1 && (parameter_argv[i + 1][0] != '-' || i == parameter_argv.size() - 1)) {
            file << " ";
        }
    }
    file << "\n";

    // Close the file
    file.close();
    return 0;
}


int Parameters::getIntParam(const char *param_name)
{
    return Parameters::getIntParam(std::string(param_name));
}


int Parameters::getIntParam(char *param_name)
{
    return Parameters::getIntParam(std::string(param_name));
}

    
int Parameters::getIntParam(std::string param_name)
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

