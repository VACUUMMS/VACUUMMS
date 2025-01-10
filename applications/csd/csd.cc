/* csd.cc */
/* input file should be of .cav format */
/* output is of .hst format */

#include <stdio.h>
#include <string.h>

#include <vacuumms/parameters.hh>
#include <vacuumms/cavity.hh>

int n_bins = 100;
double resolution = .01;
const char *input_file_name;
int histogram[1000];

int main(int argc, char *argv[])
{
    Parameters p(argc, argv);
  
    p.getStringParam((char*)"input_file_name", &input_file_name);
    p.getIntParam((char*)"n_bins", &n_bins);
    p.getDoubleParam((char*)"resolution", &resolution);

    CavityConfiguration cc(input_file_name);

    for (int i=0; i<cc.getSize(); i++)
    {
        int which_bin = (int)(cc.recordAt(i).d / resolution);
        histogram[which_bin]++;
    }

    for (int i=0; i<n_bins; i++) printf("%lf\t%d\n", i*resolution, histogram[i]);
  
    return 0;
}
