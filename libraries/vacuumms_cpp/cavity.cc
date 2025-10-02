/* vacuumms/cavity.cc */

#include <vacuumms/cavity.hh>
#include <vacuumms/parameters.hh>

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>
#include <iostream>


Cavity::Cavity()
{
}


Cavity::Cavity(vacuumms_float _x, vacuumms_float _y, vacuumms_float _z, vacuumms_float _d)
{
    x = _x;
    y = _y;
    z = _z;
    d = _d;
}


Cavity::Cavity(int _index, 
               vacuumms_float _x, 
               vacuumms_float _y, 
               vacuumms_float _z, 
               vacuumms_float _d, 
               vacuumms_float _drift)
{
    index = _index;
    x = _x;
    y = _y;
    z = _z;
    d = _d;
    drift = _drift;
}


std::vector<vacuumms_float> Cavity::getXYZ()
{
    return std::vector<vacuumms_float>{x,y,z};
}


void Cavity::setForeignKey(int _foreign_key)
{
    foreign_key = _foreign_key;
}


int Cavity::getForeignKey()
{
    return foreign_key;
}


CavityConfiguration::CavityConfiguration()
{
    records = std::vector<Cavity>();
}

#ifdef BUILD_PYBIND_BINDINGS

pybind11::str CavityConfiguration::__repr__()
{
    pybind11::str retval("");

    for (int i=0; i<records.size(); i++)
        retval = retval +
                 pybind11::str(std::to_string(records[i].x)) +
                 pybind11::str("\t") +
                 pybind11::str(std::to_string(records[i].y)) +
                 pybind11::str("\t") +
                 pybind11::str(std::to_string(records[i].z)) +
                 pybind11::str("\t") +
                 pybind11::str(std::to_string(records[i].d)) +
                 pybind11::str("\n");
    return retval;
}

#endif


CavityConfiguration::CavityConfiguration(const char *filename)
{
    FILE* instream=fopen(filename, "r");
    records = std::vector<Cavity>();
    vacuumms_float x, y, z, d;

    while (!feof(instream))
    {
        fscanf(instream, "%f\t%f\t%f\t%f\n", &x, &y, &z, &d);
        records.push_back(Cavity(x, y, z, d));
    }
}

    
CavityConfiguration::CavityConfiguration(FILE *instream)
{
    records = std::vector<Cavity>();
    vacuumms_float x, y, z, d;

    while (!feof(instream))
    {
        fscanf(instream, "%f\t%f\t%f\t%f\n", &x, &y, &z, &d);
        records.push_back(Cavity(x, y, z, d));
    }
}

    
void CavityConfiguration::replicate(std::vector<int> depths)
{
    // Use size of original vector
    size_t size = records.size();

    for (int r = 0; r < size; r++)
    {
        for (int i=0; i<depths[0]; i++)
        for (int j=0; j<depths[1]; j++)
        for (int k=0; k<depths[2]; k++)
        {
            // skip the center box
            if (!((i == 0) && (j == 0) && (k == 0)))
                pushBack(Cavity((box_dimensions[0] * i) + records[r].x,
                         (box_dimensions[1] * j) + records[r].y,
                         (box_dimensions[2] * k) + records[r].z,
                         records[r].d));
        }
    }

    // Now adjust box dimensions
    box_dimensions[0] *= depths[0];
    box_dimensions[1] *= depths[1];
    box_dimensions[2] *= depths[2];
}


//void CavityConfiguration::setBoxDimensions(vacuumms_float _box_x, vacuumms_float _box_y, vacuumms_float _box_z)
void CavityConfiguration::setBoxDimensions(std::vector<vacuumms_float> _dims)
{
    box_dimensions = _dims;
}


std::vector<vacuumms_float> CavityConfiguration::getBoxDimensions()
{
    return box_dimensions;
}


/* what was this even for?
void CavityConfiguration::setMirrorDepth(int _mirror_depth)
{
    mirror_depth = _mirror_depth;
}
*/


Cavity CavityConfiguration::recordAt(int i)
{
    return records[i];
}


void CavityConfiguration::deleteRecordAt(int i)
{
    records.erase(records.begin() + i);
}


int CavityConfiguration::getSize()
{
    return records.size();
}


int CavityConfiguration::pushBack(Cavity _cavity)
{
    records.push_back(_cavity);
    return records.size();
}


void CavityConfiguration::setDuplicateThreshold(vacuumms_float _duplicate_threshold)
{
    duplicate_threshold = _duplicate_threshold;
}


void CavityConfiguration::reset()
{
    records.clear();
}


void CavityConfiguration::scrubDuplicates()
{
    vacuumms_float box_x = box_dimensions[0];
    vacuumms_float box_y = box_dimensions[1];
    vacuumms_float box_z = box_dimensions[2];

    for (size_t index = 0; index < records.size(); index++)
    {
        for (size_t pairing = index + 1; pairing < records.size();)
        {
            // compare center to image in all adjacent mirror boxes

            for (vacuumms_float shift_x=-box_x; shift_x<=box_x; shift_x += box_x)
            for (vacuumms_float shift_y=-box_y; shift_y<=box_y; shift_y += box_y)
            for (vacuumms_float shift_z=-box_z; shift_z<=box_z; shift_z += box_z)
            {
                vacuumms_float dsq = (shift_x + records[pairing].x - records[index].x) * (shift_x + records[pairing].x - records[index].x)
                                   + (shift_y + records[pairing].y - records[index].y) * (shift_y + records[pairing].y - records[index].y)
                                   + (shift_z + records[pairing].z - records[index].z) * (shift_z + records[pairing].z - records[index].z);

                if (dsq < duplicate_threshold)
                {
                    records.erase(records.begin() + pairing);
                    goto mirrors_done; // no need to keep looking, we know it's a duplicate
                }
            }

            // No duplicate found so move on to next pairing;
            pairing++;

            mirrors_done:
            ;
        }
    }
}


std::vector<vacuumms_float> CavityConfiguration::getDiameters()
{
    std::vector<vacuumms_float> diameters;
    for (const auto& record : records) 
    {
        diameters.push_back(record.d);
    }

    return diameters;
}


int CavityConfiguration::checkInclusion(vacuumms_float tx, vacuumms_float ty, vacuumms_float tz)
{
    int i;
    vacuumms_float dx, dy, dz, dd;

    for (i=0; i<getSize(); i++)
    {
        vacuumms_float box_x = box_dimensions[0];
        vacuumms_float box_y = box_dimensions[1];
        vacuumms_float box_z = box_dimensions[2];

        // Check inclusion in each of the eight mirror box images:

        // (0,0,0):
        dx = records[i].x - tx;
        dy = records[i].y - ty;
        dz = records[i].z - tz;
        dd = dx*dx + dy*dy + dz*dz;
        if (4*dd < (records[i].d * records[i].d)) return 1;

        // (0,0,1):
        dx = records[i].x - tx;
        dy = records[i].y - ty;
        dz = box_z + records[i].z - tz;
        dd = dx*dx + dy*dy + dz*dz;
        if (4*dd < (records[i].d * records[i].d)) return 1;

        // (0,1,0):
        dx = records[i].x - tx;
        dy = box_y + records[i].y - ty;
        dz = records[i].z - tz;
        dd = dx*dx + dy*dy + dz*dz;
        if (4*dd < (records[i].d * records[i].d)) return 1;

        // (0,1,1):
        dx = records[i].x - tx;
        dy = box_y + records[i].y - ty;
        dz = box_z + records[i].z - tz;
        dd = dx*dx + dy*dy + dz*dz;
        if (4*dd < (records[i].d * records[i].d)) return 1;

        // (1,0,0):
        dx = box_x + records[i].x - tx;
        dy = records[i].y - ty;
        dz = records[i].z - tz;
        dd = dx*dx + dy*dy + dz*dz;
        if (4*dd < (records[i].d * records[i].d)) return 1;

        // (1,0,1):
        dx = box_x + records[i].x - tx;
        dy = records[i].y - ty;
        dz = box_z + records[i].z - tz;
        dd = dx*dx + dy*dy + dz*dz;
        if (4*dd < (records[i].d * records[i].d)) return 1;

        // (1,1,0):
        dx = box_x + records[i].x - tx;
        dy = box_y + records[i].y - ty;
        dz = records[i].z - tz;
        dd = dx*dx + dy*dy + dz*dz;
        if (4*dd < (records[i].d * records[i].d)) return 1;

        // (1,1,1):
        dx = box_x + records[i].x - tx;
        dy = box_y + records[i].y - ty;
        dz = box_z + records[i].z - tz;
        dd = dx*dx + dy*dy + dz*dz;
        if (4*dd < (records[i].d * records[i].d)) return 1;
    }

    return 0;
}


CavitySizeDistribution::CavitySizeDistribution(CavityConfiguration cc, Parameters p) 
    : cc(cc), p(p)
{
    // set up bins and sizes
    if (p.getFlagParam((char*)"-n_bins"))
    {
        setNumberOfBins(p.getIntParam((char*)"-n_bins"));
    }
    if (p.getFlagParam((char*)"-width"))
    {
        setBinWidth(p.getFloatParam((char*)"-width"));
    }

    // // implement later
    // vacuumms_float start_x = p.getFloatParam((char*)"-start_x"); 

    for (int i=0; i<cc.getSize(); i++)
    {
        bin(cc.recordAt(i).d);
//        int which_bin = (int)(cc.recordAt(i).d / width_of_bins);
//        histogram[which_bin]++;
    }
}


CavitySizeDistribution::CavitySizeDistribution(CavityConfiguration cc) 
    : cc(cc)
{
    for (int i=0; i<cc.getSize(); i++) bin(cc.recordAt(i).d);
}


/* to Histogram class
void CavitySizeDistribution::setBinWidth(vacuumms_float _width)
{
    width = _width;
}


void CavitySizeDistribution::setNumberOfBins(int _n_bins)
{
    n_bins = _n_bins;
}
*/


std::vector<std::tuple<vacuumms_float, vacuumms_float>> CavitySizeDistribution::getResult()
{
    for (int i=0; i<cc.getSize(); i++) bin(cc.recordAt(i).d);

    std::vector<std::tuple<vacuumms_float, vacuumms_float>> result;

    for (int i = 1; i <= bins.size(); ++i) 
    {
//        result.emplace_back(i, i * 2); // Tuple of (i, i*2)
        result.emplace_back(i, i * 2); // Tuple of (i, i*2)
    }
    return result;
}

        
#ifdef BUILD_PYBIND_BINDINGS

pybind11::str CavitySizeDistribution::__repr__()
{
    return Histogram::__repr__();
/*
    pybind11::str retval("");

//    for (int i=0; i<records.size(); i++)
    for (int i=0; i<number_of_bins; i++) 
//        printf("%lf\t%d\n", i*resolution, histogram[i]);
        retval = retval +
                 pybind11::str(std::to_string((vacuumms_float)(i * width_of_bins))) +
                 pybind11::str("\t") +
                 pybind11::str(std::to_string(bins[i])) +
                 pybind11::str("\n");
    return retval;
*/
}

#endif
