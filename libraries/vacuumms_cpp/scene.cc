/* libraries/vacuumms_cpp/scene.cc */

#include <vector>
#include <iostream>

#include <vacuumms/scene.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>

std::string Scene::generateContainerSDL()
{
	std::stringstream out;

	// headers
    out << "#include \"colors.inc\"\n";
    out << "{color " << background_color << "}\n";
    out << "{location<" << camera_location[0] << "," << camera_location[1] << "," << camera_location[2] << "> look_at <0,0,0>}\n";

	if (ambient_light) out << "global_settings { ambient_light rgb <" << ambient_light << "," << ambient_light << "," << ambient_light << "> }\n"; 

    // apply light sources
    std::string light_color = "White";
//    for (auto source = light_sources.begin(); source != light_sources.end(); ++source)
    for (size_t i = 0; i < light_sources.size(); i++)
    {
        std::vector<vacuumms_float> source = light_sources[i];
		out << "light_source{<" << source[0] 
			<< "," << source[1] 
			<< "," << source[2] 
			<< "> color " << light_color << "}\n";
    }
    
/*
    if (standard_light)
    {
		out << "light_source{<100,0,0> color " << light_color << "}\n";
		out << "light_source{<0,100,0> color " << light_color << "}\n";
		out << "light_source{<0,0,100> color " << light_color << "}\n";
		out << "light_source{<-100,0,0> color " << light_color << "}\n";
		out << "light_source{<0,-100,0> color " << light_color << "}\n";
		out << "light_source{<0,0,-100> color " << light_color << "}\n";
    }

    if (light_source) 
		out << "light_source{<" << light_source_x 
			<< "," << light_source_y 
			<< "," << light_source_z 
			<< "> color " << light_color << "}\n";
*/

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
	}

	return out.str();
}

ConfigurationComponent::ConfigurationComponent(Configuration configuration)
{
}

std::string ConfigurationComponent::getComponentSDL()
{
    return "foo";
}
    

CavityComponent::CavityComponent(CavityConfiguration configuration)
{
}

FVIComponent::FVIComponent(FVIX fvix)
{
}



void Scene::setBoxDimensions(std::vector<vacuumms_float> dims)
{
	box_dimensions = dims;
}

std::vector<vacuumms_float> Scene::getBoxDimensions()
{
	return box_dimensions;
}

SceneComponent Scene::componentAt(int i)
{
	return components[i];
}

size_t Scene::deleteComponentAt(int i)
{
    components.erase(components.begin() + i);
    return components.size();
}

size_t Scene::getSize()
{
	return components.size();
}

size_t Scene::pushBack(SceneComponent comp)
{
	components.push_back(comp);
	return components.size();
}

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

// I/O
int Scene::createSceneFile(const char* filename)  // POV file
{
    return 0;
}

int Scene::renderScene(const char* filename)      // PNG file
{
    return 0;
}

/*
#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif
*/


