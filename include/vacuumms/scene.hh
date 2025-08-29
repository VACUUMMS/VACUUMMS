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

        virtual std::string getComponentSDL() const;
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


// Trampoline class for SceneComponent
class PySceneComponent : public SceneComponent {
public:
    using SceneComponent::SceneComponent; // Inherit constructors
    std::string getComponentSDL() const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,    // Return type
            SceneComponent, // Parent class
            getComponentSDL // Function name
        );
    }
};


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
FTWComponent: public SceneComponent
{
    public:

        FTWComponent();
        std::string getComponentSDL();
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
        std::string getComponentSDL() const;
    
    private:
    
        vacuumms_float transmit;
        vacuumms_float phong;
        std::string color = "Red";
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
        std::string getComponentSDL() const;
    
    private:
    
        vacuumms_float transmit;
        vacuumms_float phong;
        std::string color = "White";
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
        int dumpSDL();  // dump POV source to stdout
        int createSceneFile(const char* filename);  // POV file
        int renderScene(const char* filename);      // PNG file
        std::string generateContainerSDL();

//FTW        SceneComponent componentAt(int i);
        SceneComponent* componentAt(int i);
        size_t deleteComponentAt(int i);
        size_t getNumberOfComponents();
//FTW        size_t addSceneComponent(SceneComponent);
        size_t addSceneComponent(SceneComponent*);
//FTW        size_t addSceneComponent(auto*);
//        size_t addConfigurationComponent(ConfigurationComponent*);
//        size_t addCavityComponent(CavityComponent*);
//        size_t addFTWComponent(FTWComponent*);

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
//FTW        std::vector<SceneComponent> components;
        std::vector<SceneComponent*> components;
        
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

