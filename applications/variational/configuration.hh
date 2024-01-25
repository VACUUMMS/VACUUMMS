#include <vector>

class ConfigurationRecord
{
    public:

        float x;
        float y;
        float z;
        float sigma;
        float epsilon;
//        float sigma_6;
//        float sigma_12;

        ConfigurationRecord(float _x, float _y, float _z, float _sigma, float _epsilon);
};

class Configuration
{
        std::vector<ConfigurationRecord> records;

        float box_x;
        float box_y;
        float box_z;

        int mirror_depth = 1;
    
    public:

        Configuration(char *filename);
        Configuration();
        void dumpContents();
        float insertionEnergy(float x, float y, float z, float sigma, float epsilon);
        void setBoxDimensions(float _box_x, float _box_y, float _box_z);
        void setMirrorDepth(int _mirror_depth);

        ConfigurationRecord recordAt(int i);
        void deleteRecordAt(int i);
        int getSize();
};


