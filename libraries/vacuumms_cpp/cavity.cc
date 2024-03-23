/* vacuumms/cavity.cc */

#include <vacuumms/cavity.hh>

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>
#include <iostream>


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

CavityConfiguration::CavityConfiguration(char *filename)
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
    
void CavityConfiguration::setBoxDimensions(vacuumms_float _box_x, vacuumms_float _box_y, vacuumms_float _box_z)
{
    box_x = _box_x;
    box_y = _box_y;
    box_z = _box_z;
}

void CavityConfiguration::setMirrorDepth(int _mirror_depth)
{
    mirror_depth = _mirror_depth;
}

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

int CavityConfiguration::checkInclusion(vacuumms_float tx, vacuumms_float ty, vacuumms_float tz)
{
    int i;
    vacuumms_float dx, dy, dz, dd;

    for (i=0; i<getSize(); i++)
    {
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

