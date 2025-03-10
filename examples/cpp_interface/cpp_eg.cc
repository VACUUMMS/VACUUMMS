#include <vacuumms/types.hh>
#include <cstdio>


#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/parameters.hh>
#include <vacuumms/pddx.hh>

int main()
{
    // read file and perform CESA
    Configuration c("fcc.gfg");
    std::vector<std::string> argv = {"-n_threads", "4", "-n", "100", "-box", "4.24264", "4.24264", "4.24264"};
    Parameters p(argv);
    PDDX o(c,p);
    o.execute();
    CavityConfiguration result=o.getResult();

    // Take the result, bin it, and smooth three times
    std::vector<std::string> list = {"-width", "0.1", "-n_bins", "50"};
    Parameters pp(list); 
    CavitySizeDistribution csd(result, pp);
    csd.smooth(3);

    // Display output, and also write it to file. 
    csd.print();
    csd.writeToFile((char*)"csd.hst");
    
    return 0;
}
