/* libraries/vacuumms_cpp/scene.cc */

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>

#include <vacuumms/scene.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>

std::string Scene::generateContainerSDL()
{
	std::stringstream out;

	// headers
	
    out << "#include \"colors.inc\"\n";

    out << "background {color " << background_color << "}\n";

    out << "camera {location <" << camera_location[0] 
        << "," << camera_location[1] 
        << "," << camera_location[2] 
        << "> look_at <" << camera_look_at[0]
        << "," << camera_look_at[1] 
        << "," << camera_look_at[2] 
        << "> right 1.0 angle 45}\n";

    // ambient light

	if (ambient_light) out << "global_settings { ambient_light rgb <" << ambient_light << "," << ambient_light << "," << ambient_light << "> }\n"; 

    // apply other light sources
    
    std::string light_color = "White";

    for (size_t i = 0; i < light_sources.size(); i++)
    {
        std::vector<vacuumms_float> source = light_sources[i];
		out << "light_source{<" << source[0] 
			<< "," << source[1] 
			<< "," << source[2] 
			<< "> color " << light_color << "}\n";
    }
    
    // box
    
    if (show_box)
    {
		out << "cylinder { <0,0,0>, <" << box_dimensions[0] 
			<< ",0,0>, 0.1 open texture { pigment { color " << box_color 
			<< " } }}\n";
		out << "cylinder { <0,0,0>, <0," << box_dimensions[1] 
			<< ",0>, 0.1 open texture { pigment { color " << box_color 
			<< " } }}\n";
		out << "cylinder { <0,0,0>, <0,0," << box_dimensions[2] 
			<< ">, 0.1 open texture { pigment { color " << box_color 
			<< " } }}\n";
		out << "cylinder { <" << box_dimensions[0] 
			<< "," << box_dimensions[1] << ",0>, <" 
			<< box_dimensions[0] << ",0,0>, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_dimensions[0] 
			<< "," << box_dimensions[1] << ",0>, <0," 
			<< box_dimensions[1] << ",0>, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_dimensions[0] << "," 
			<< box_dimensions[1] << ",0>, <" << box_dimensions[0] << "," 
			<< box_dimensions[1] << "," << box_dimensions[2] << ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <0," << box_dimensions[1] 
			<< "," << box_dimensions[2] << ">, <0," << box_dimensions[1] 
			<< ",0>, 0.1 open texture { pigment { color " << box_color 
			<< " } }}\n";
		out << "cylinder { <0," << box_dimensions[1] 
			<< "," << box_dimensions[2] << ">, <" << box_dimensions[0] << "," << box_dimensions[1] 
			<< "," << box_dimensions[2] << ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <0," << box_dimensions[1] << "," << box_dimensions[2] 
			<< ">, <0,0," << box_dimensions[2] << ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_dimensions[0] << ",0," << box_dimensions[2] 
			<< ">, <" << box_dimensions[0] << ",0,0>, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_dimensions[0] << ",0," << box_dimensions[2] 
			<< ">, <0,0," << box_dimensions[2] << ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_dimensions[0] << ",0," << box_dimensions[2] 
			<< ">, <" << box_dimensions[0] << "," << box_dimensions[1] << "," << box_dimensions[2] 
			<< ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";

		out << "sphere{<0, 0, 0>, .1 texture{ pigment {color %s}}}\n";
		out << "sphere{<0, 0, " << box_dimensions[2] 
			<< ">, .1 texture{ pigment {color %s}}}\n";
		out << "sphere{<0, " << box_dimensions[1] 
			<< " ,0>, .1 texture{ pigment {color %s}}}\n";
		out << "sphere{<0, " << box_dimensions[1] 
			<< " ," << box_dimensions[2] << ">, .1 texture{ pigment {color %s}}}\n";
		out << "sphere{<" << box_dimensions[0] 
			<< ", 0, 0>, .1 texture{ pigment {color %s}}}\n";
		out << "sphere{<" << box_dimensions[0] << ", 0, " 
			<< box_dimensions[2] << ">, .1 texture{ pigment {color %s}}}\n";
		out << "sphere{<" << box_dimensions[0] << ", " << box_dimensions[1] 
			<< ", 0>, .1 texture{ pigment {color %s}}}\n";
		out << "sphere{<" << box_dimensions[0] << ", " << box_dimensions[1] << ", " 
			<< box_dimensions[2] << ">, .1 texture{ pigment {color %s}}}\n";

        // end of SDL comment
        
        out << "# end of sceneSDL\n\n";
	}

	return out.str();
}


void SceneComponent::setClipComponent(int _clip)
{
    clip = _clip;
}


void SceneComponent::setTransmit(vacuumms_float _transmit)
{
    transmit = _transmit;
}


void SceneComponent::setPhong(vacuumms_float _phong)
{
    phong = _phong;
}


void SceneComponent::setColor(std::string _color)
{
    color = _color;
}

std::string SceneComponent::getComponentSDL() const
{
    // Return an unit orange bubble centered at origin as default. 
    // return std::string("sphere{<0.0, 0.0, 0.0>, 1.0 texture{ pigment {color Orange  transmit 0.700000  }  finish {phong 0.700000}  } }\n");
    return std::string("// SceneComponent::getComponentSDL called on base type.\n\n");
}

ConfigurationComponent::ConfigurationComponent(){}


ConfigurationComponent::ConfigurationComponent(Configuration _configuration)
{
    configuration = _configuration;
}


std::string ConfigurationComponent::getComponentSDL() const
{
    std::stringstream sdl;
    sdl << "// Configuration Component\n\n";

    std::string transmit_str = " transmit " + std::to_string(transmit);
    std::string phong_str = " finish {phong " + std::to_string(phong) + "} ";

    for (const auto record : configuration.records) 
    {

printf("dumping configuration record: %f\t%f\t%f\t%f\t%f\n", record.x, record.y, record.z, record.sigma, record.epsilon);

        if (clip) 
        {
            sdl << "intersection {sphere{<" 
                << record.x
                << ", " << record.y 
                << ", " << record.z
                << ">, " << (record.sigma * 0.5)
                << "} box {<0,0,0><" << configuration.box_dimensions[0] 
                << ", " << configuration.box_dimensions[1]
                << ", " << configuration.box_dimensions[2]
                << ">} texture { pigment { color "
                << color << " " << transmit_str 
                << " } " << phong_str << " }}\n"
                ;
        }
        else 
        {
            sdl << "sphere{<" 
                << record.x
                << ", " << record.y
                << ", " << record.z
                << ">, " << (record.sigma *0.5) 
                << "texture{ pigment {color " 
                << color << " " << transmit_str 
                << " } " << phong_str << " } }\n"
                ;
        }
    }

    // end of SDL comment
    sdl << "// end of ConfigurationComponent SDL\n\n";
    return sdl.str();
}
    

CavityComponent::CavityComponent(){}


CavityComponent::CavityComponent(CavityConfiguration _configuration)
{
    configuration = _configuration;
}


std::string CavityComponent::getComponentSDL() const
{
//    return std::string("// SceneComponent::getComponentSDL called on CavityComponent\n\n");
    std::stringstream sdl;
    sdl << "// Cavity Component\n\n";

    std::string transmit_str = " transmit " + std::to_string(transmit);
    std::string phong_str = " finish {phong " + std::to_string(phong) + "} ";

    for (const auto record : configuration.records) 
    {

printf("dumping configuration record: %f\t%f\t%f\t%f\n", record.x, record.y, record.z, record.d);

        if (clip) 
        {
            sdl << "intersection {sphere{<" 
                << record.x
                << ", " << record.y 
                << ", " << record.z
                << ">, " << (record.d * 0.5)
                << "} box {<0,0,0><" << configuration.box_dimensions[0] 
                << ", " << configuration.box_dimensions[1]
                << ", " << configuration.box_dimensions[2]
                << ">} texture { pigment { color "
                << color << " " << transmit_str 
                << " } " << phong_str << " }}\n"
                ;
        }
        else 
        {
            sdl << "sphere{<" 
                << record.x
                << ", " << record.y
                << ", " << record.z
                << ">, " << (record.d * 0.5) 
                << "texture{ pigment {color " 
                << color << " " << transmit_str 
                << " } " << phong_str << " } }\n"
                ;
        }
    }

    // end of SDL comment
    sdl << "// end of CavityComponent SDL\n\n";
    return sdl.str();
}


#ifdef BUILD_CUDA_COMPONENTS
FVIComponent::FVIComponent(FVIX fvix)
{
}
#endif

void Scene::setBoxDimensions(std::vector<vacuumms_float> dims)
{
	box_dimensions = dims;
}

std::vector<vacuumms_float> Scene::getBoxDimensions()
{
	return box_dimensions;
}

SceneComponent* Scene::componentAt(int i)
{
	return components[i];
}

size_t Scene::deleteComponentAt(int i)
{
    components.erase(components.begin() + i);
    return components.size();
}

size_t Scene::getNumberOfComponents()
{
	return components.size();
}


size_t Scene::addSceneComponent(SceneComponent* comp)
{
	components.push_back(comp);
	return components.size();
}


/*
size_t Scene::addConfigurationComponent(ConfigurationComponent* comp)
{
	components.push_back(comp);
	return components.size();
}


size_t Scene::addCavityComponent(CavityComponent* comp)
{
	components.push_back(comp);
	return components.size();
}
*/


void Scene::setBackgroundColor(std::string color)
{
	background_color = color;
}

void Scene::setCameraLocation(std::vector<vacuumms_float> location)
{
	camera_location = location;
}

size_t Scene::addLightSource(std::vector<vacuumms_float> source, std::string color)
{
	light_sources.push_back(source);
    return light_sources.size();
}

void Scene::applyAmbientLight()
{
    ambient_light = 1;
}

size_t Scene::applyStandardLight()
{
	addLightSource({0,0,100}, "White");
	addLightSource({0,100,0}, "White");
	addLightSource({100,0,0}, "White");
	addLightSource({0,0,-100}, "White");
	addLightSource({0,-100,0}, "White");
	addLightSource({-100,0,0}, "White");
    return light_sources.size();
}

void Scene::setShowBox(int yn)
{
	show_box = yn;
}

void Scene::setBoxColor(std::string color)
{
	box_color = color;
}

int Scene::dumpSDL()  // dump POV source
{
    std::stringstream scene;

    // Container

    scene << generateContainerSDL();
    
    // Components
        
    for (const auto* obj : components) {
        scene << obj->getComponentSDL(); // Calls the appropriate version
        scene << std::endl;
    }

    std::cout << scene.str();

    return 0;
}

// I/O
int Scene::createSceneFile(const char* filename)  // POV file
{
    std::ofstream scene_file(filename);
    if (scene_file.is_open()) // write it
    {
        // Container

        scene_file << generateContainerSDL();
        
        // Components
        
        for (int i=0; i < components.size(); i++)
        {
            scene_file << "// writing component " << i << std::endl;
            scene_file << components[i]->getComponentSDL();
            scene_file << std::endl;
        }

        scene_file.close();

        return 0;
    }
    else
    {
        std::cout << "Could not write file " << filename << std::endl;
        return 1;
    }
}

int Scene::renderScene(const char* filename)      // PNG file
{
    std::string basename = std::filesystem::path(filename).stem().string();
    std::string pov_filename = basename + ".pov";

    createSceneFile(pov_filename.c_str());

    std::string command = "povray -W1920 -H1080 " + pov_filename;
    std::string rm_command = "rm -f " + pov_filename;

    // Render
    int result = std::system(command.c_str());
    if (result == 0) 
	{
        std::cout << "POV-Ray render completed successfully.\n";
    } 
	else 
	{
        std::cerr << "POV-Ray render failed with exit code: " << result << "\n";
    }

    // Clean up
    int rm_result = std::system(rm_command.c_str());
    if (rm_result == 0) 
	{
        std::cout << "POV-Ray temporary file deleted successfully.\n";
    } 
	else 
	{
        std::cerr << "POV-Ray temporary file could not be deleted, failed with exit code: " << rm_result << "\n";
    }

    return result;
}

/*
#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif
*/


