#include <vacuumms/variational/variational.hh>
#include <vacuumms/variational/variational_scene.hh>
#include <vacuumms/scene.hh>

VariationalComponent::VariationalComponent()
{
}


VariationalComponent::VariationalComponent(Variational3D* _variational)
{
    variational = _variational;
}

std::string VariationalComponent::getComponentSDL() const
{
    if (hidden) return "// Variational Component hidden\n\n";

    std::stringstream sdl;
    sdl << "// Variational Component\n\n";

    // std::string transmit_str = " transmit " + std::to_string(transmit);
    // std::string phong_str = " finish {phong " + std::to_string(phong) + "} ";

    vacuumms_float* var_x = variational->getX();
    vacuumms_float* var_y = variational->getY();
    vacuumms_float* var_z = variational->getZ();

    // draw as cylinders for a continuous trajectory
    // format: cylinder{<1.0493,0.0344035,0.224887>,<1.23868,0.77789,-0.232467>, 0.0025 texture{ pigment {color Yellow } }}
    
    for (int point = 0; point < variational->getNVariationalPoints() - 1; point++)
    {
        sdl << "cylinder{<" << var_x[point] << ", " << var_y[point] << ", " << var_z[point] << ">, "
            << "<" << var_x[point + 1] << ", " << var_y[point + 1] << ", " << var_z[point + 1] << ">, " 
            << cylinder_diameter << " texture{ pigment {color Yellow } }}\n";
    }

    sdl << "// end of VariationalComponent SDL\n\n";
    return sdl.str();
}


