/* vacuumms/cuda.hh */

#include <vacuumms/operations.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/parameters.hh>

#include <vacuumms/exports.hh>

class FVIX : public Operation
{
    public:
        
        FVIX();
        FVIX(Configuration, Parameters);
        void printUsage();
        void setParameters();
        void setConfiguration();
        void execute();
        // getResult();

#ifdef BUILD_PYBIND_BINDINGS
        void __repr__();
#endif

    private:
        
        Parameters p;
        Configuration c;
};
 

