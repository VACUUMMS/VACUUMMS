/* vacuumms/scene.hh */

#pragma once

#include <vector>

#include <vacuumms/types.h>

#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>

#ifdef BUILD_CUDA_COMPONENTS
    #include <vacuumms/fvi.hh>
#endif

#include <vacuumms/exports.hh>


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
SceneComponent
{
    public:

        std::string getComponentSDL();
        void setClipComponent(int);
        void setPhong(vacuumms_float);
        void setTransmit(vacuumms_float);
        void setColor(std::string);

    private:

        int clip = 0;
        vacuumms_float phong = 0.0;
        vacuumms_float transmit = 0.0;

    protected:

        std::string color;
};

class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
ConfigurationComponent : public SceneComponent
{
    public:

        ConfigurationComponent();
        ConfigurationComponent(Configuration);
        std::string getComponentSDL();
    
    private:
    
        vacuumms_float transmit;
        vacuumms_float phong;
        std::string color;
        std::vector<vacuumms_float> box_dims;
        int clip; // intersect with box
        Configuration configuration;
};

class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
CavityComponent : public SceneComponent
{
    public:

        CavityComponent();
        CavityComponent(CavityConfiguration);
        std::string getComponentSDL();
    
    private:
    
        vacuumms_float transmit;
        vacuumms_float phong;
        std::string color;
        std::vector<vacuumms_float> box_dims;
        int clip; // intersect with box
        CavityConfiguration configuration;
};

#ifdef BUILD_CUDA_COMPONENTS

class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
FVIComponent : public SceneComponent
{
    public:

        FVIComponent(FVIX fvix);
};

#endif


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
Scene
{
    public:

        // I/O
        int createSceneFile(const char* filename);  // POV file
        int renderScene(const char* filename);      // PNG file
        std::string generateContainerSDL();

        SceneComponent componentAt(int i);
        size_t deleteComponentAt(int i);
        size_t getNumberOfComponents();
        size_t addSceneComponent(SceneComponent);

        void setBoxDimensions(std::vector<vacuumms_float>);
        std::vector<vacuumms_float> getBoxDimensions();

        void setBackgroundColor(std::string);
        void setCameraLocation(std::vector<vacuumms_float>);
        void setCameraLookAt(std::vector<vacuumms_float>);
        size_t addLightSource(std::vector<vacuumms_float>, std::string color);
        size_t applyStandardLight();
        void applyAmbientLight();
        void setShowBox(int);
        void setBoxColor(std::string);

    private:

        int ambient_light = 0;
        int show_box = 0;

        std::vector<std::vector<vacuumms_float>> light_sources;
        std::vector<SceneComponent> components;
        
        std::vector<vacuumms_float> camera_location = {40, 40, 40};
        std::vector<vacuumms_float> camera_look_at = {0, 0, 0};
        std::vector<vacuumms_float> box_dimensions = {10, 10, 10};

        std::string light_color = "White";
        std::string box_color = "Yellow";
        std::string background_color = "Black";

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

};

