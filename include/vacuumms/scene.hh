// configuration.hh 
#pragma once

#include <vector>
#include <vacuumms/types.h>

#include <vacuumms/exports.hh>


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
SceneComponent
{
    public:

        int type;

        SceneComponent();

        std::string generateSDL();
};

class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
Scene
{
    public:
   
        vacuumms_float box_x;
        vacuumms_float box_y;
        vacuumms_float box_z;
    
        vacuumms_float camera_x;
        vacuumms_float camera_y;
        vacuumms_float camera_z;

        std::vector<vacuumms_float> box_dimensions;

        std::vector<SceneComponent> records;
        int mirror_depth = 1;
        int replication_depth = 0;

        Scene(const char *filename);
        Scene(FILE *pipe); // allows stdin to be used to create pipeline
        Scene();
        void dumpContents();
        void setBoxDimensions(std::vector<vacuumms_float> dims);
        std::vector<vacuumms_float> getBoxDimensions();
        void setMirrorDepth(int _mirror_depth);

        // I/O
        createSceneFile(const char* filename);  // POV file
        renderScene(const char* filename);      // PNG file

        SceneComponent recordAt(int i);
        void deleteComponentAt(int i);
        int getSize();
        int pushBack(SceneComponent);

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

};

