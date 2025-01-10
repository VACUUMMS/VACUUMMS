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
        //char **command_line_argv;
        std::vector<std::string> command_line_argv;
        //std::vector<char*> command_line_argv;

    public:

        Parameters();
        Parameters(int argc, char **argv);

#ifdef BUILD_BOOST_PYTHON_BINDINGS 
        Parameters(const boost::python::list&);
#endif

        int addParameter(const char* parameter);

	/* if a parameter is received, return a true value, otherwise return NULL) */
#ifdef BUILD_BOOST_PYTHON_BINDINGS 
	int getIntParam(boost::python::object obj);
#else
	int getIntParam(char *param_name, int *parameter);
#endif
	int getLongParam(char *param_name, long *parameter);
	int getFloatParam(char *param_name, float *parameter);
	int getDoubleParam(char *param_name, double *parameter);
	int getStringParam(char *param_name, const char **parameter);
	//FTWint getStringParam(char *param_name, char **parameter);
	//void getStringParam(char *param_name, char *parameter);
	int getVectorParam(char *param_name, double *parameter1,  double *parameter2, double *parameter3);
	int getVectorStringParam(char *param_name, const char **parameter1, const char **parameter2, const char **parameter3);
	int getFlagParam(char *param_name);
};

