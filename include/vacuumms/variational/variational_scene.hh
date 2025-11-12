// variational_scene.hh

#pragma once

#include <vacuumms/types.h>

#include <vacuumms/scene.hh>
#include <vacuumms/configuration.hh>

#include <vacuumms/exports.hh>


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
VariationalComponent : public SceneComponent
{
    public:

        VariationalComponent();
        VariationalComponent(Variational3D*);
        void setDiameter(vacuumms_float);
        std::string getComponentSDL() const;
    
    private:
    
        Variational3D* variational;
        vacuumms_float cylinder_diameter = 0.01f;
};

