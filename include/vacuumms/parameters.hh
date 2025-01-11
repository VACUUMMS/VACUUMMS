/* vacuumms/parameters.hh */

#ifdef BUILD_BOOST_PYTHON_BINDINGS 
#include <boost/python.hpp>
#endif

#include <string>
#include <vector>

class Parameters
{
    private:

        int command_line_argc;
        std::vector<std::string> command_line_argv;

    public:

        Parameters();
        Parameters(int argc, char **argv);

        // Maintain a set of methods for boost::python. 
        // Note that these return the value instead of setting pointed value.

#ifdef BUILD_BOOST_PYTHON_BINDINGS 
        Parameters(const boost::python::list&);
	int getIntParam(boost::python::object obj);
	long getLongParam(boost::python::object obj);
	float getFloatParam(char *param_name);
	double getDoubleParam(char *param_name);
        const char* getStringParam(char *param_name);
        //std::vector<double> getVectorParam(char* param_name);
        boost::python::list getVectorParam(char* param_name);
        //std::vector<std::string> getVectorStringParam(char* param_name);
        boost::python::list getVectorStringParam(char* param_name);
#endif

        int addParameter(const char* parameter);

	/* if a parameter is received, return a true value, otherwise return NULL) */
	int getIntParam(char *param_name, int *parameter);
	int getLongParam(char *param_name, long *parameter);
	int getFloatParam(char *param_name, float *parameter);
	int getDoubleParam(char *param_name, double *parameter);
	int getStringParam(char *param_name, const char **parameter);
	int getVectorParam(char *param_name, 
                                 double *parameter1,  
                                 double *parameter2, 
                                 double *parameter3);
	int getVectorStringParam(char *param_name, 
                                 const char **parameter1, 
                                 const char **parameter2, 
                                 const char **parameter3);

        int getFlagParam(char *param_name);

};

