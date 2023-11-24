#include <vector>

class ConfigurationRecord
{
    public:

        float x;
        float y;
        float z;
        float sigma;
        float epsilon;
        float sigma_6;
        float sigma_12;

        ConfigurationRecord(float _x, float _y, float _z, float _sigma, float _epsilon);
};

class Configuration
{
        std::vector<ConfigurationRecord> records;
        float *x;
        float *y;
        float *z;

        float *sigma;
        float *epsilon;
    
    public:

        Configuration(char *filename);
        Configuration();
        void dumpContents();
        float insertionEnergy2D(float x, float y);
        float insertionEnergy3D(float x, float y, float z);
};

