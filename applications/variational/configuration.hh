#include <vector>

class ConfigurationRecord
{
        float x;
        float y;
        float z;
        float sigma;
        float epsilon;

    public:

        ConfigurationRecord(float _x, float _y, float _z, float _sigma, float _epsilon);
};

//class Configuration: private std::Vector<ConfigurationRecord>
class Configuration: private std::vector<ConfigurationRecord>
{
        float *x;
        float *y;
        float *z;

        float *sigma;
        float *epsilon;
    
    public:

        Configuration(char *filename);
};
