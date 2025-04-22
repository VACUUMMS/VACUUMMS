/* libraries/vacuumms_cpp/scene.cc */

setBackgroundColor(std::string);
setCamera(std::vector<vacuumms_float>);
setBoxDimensions(std::vector<vacuumms_float>);
addLightSource(std::vector<vacuumms_float>, std::string color);
useAmbientLight();
useStandardLight();
setBackground(std::string);
setShowBox(int);
setBoxColor(std::string);

double camera_x=40, camera_y=40, camera_z=40;
double light_source_x=25, light_source_y=25, light_source_z=25;
double box_x=10, box_y=10, box_z=10;
double ambient_light = 3.0;

char *color = "Red";
char *light_color = "White";
char *box_color = "Yellow";
char *background = "Black";

std::string Scene::generateSDL()
{
	std::stringstream out;

	// headers
    out << "#include \"colors.inc\"\n";
    out << "{color " << background << "}\n";
    out << "{location<" << camera_x << "," << camera_y << "," << camera_z << "> look_at <0,0,0>}\n";

	if (ambient_light) out << "global_settings { ambient_light rgb <" << ambient_light << "," << ambient_light << "," << ambient_light << "> }\n"; 

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

    if (show_box)
    {
		out << "cylinder { <0,0,0>, <" << box_x 
			<< ",0,0>, 0.1 open texture { pigment { color " << box_color 
			<< " } }}\n", box_x, box_color;
		out << "cylinder { <0,0,0>, <0," << box_y 
			<< ",0>, 0.1 open texture { pigment { color " << box_color 
			<< " } }}\n";
		out << "cylinder { <0,0,0>, <0,0," << box_z 
			<< ">, 0.1 open texture { pigment { color " << box_color 
			<< " } }}\n";
		out << "cylinder { <" << box_x 
			<< "," << box_y << ",0>, <" 
			<< box_x << ",0,0>, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_x 
			<< "," << box_y << ",0>, <0," 
			<< box_y << ",0>, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_x << "," 
			<< box_y << ",0>, <" << box_x << "," 
			<< box_y << "," << box_z << ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <0," << box_y 
			<< "," << box_z << ">, <0," << box_y 
			<< ",0>, 0.1 open texture { pigment { color " << box_color 
			<< " } }}\n";
		out << "cylinder { <0," << box_y 
			<< "," << box_z << ">, <" << box_x << "," << box_y 
			<< "," << box_z << ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <0," << box_y << "," << box_z 
			<< ">, <0,0," << box_z << ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_x << ",0," << box_z 
			<< ">, <" << box_x << ",0,0>, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_x << ",0," << box_z 
			<< ">, <0,0," << box_z << ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";
		out << "cylinder { <" << box_x << ",0," << box_z 
			<< ">, <" << box_x << "," << box_y << "," << box_z 
			<< ">, 0.1 open texture { pigment { color " 
			<< box_color << " } }}\n";

		out << "sphere{<0, 0, 0>, .1 texture{ pigment {color %s}}}\n":
		out << "sphere{<0, 0, " << box_z 
			<< ">, .1 texture{ pigment {color %s}}}\n":
		out << "sphere{<0, " << box_y 
			<< " ,0>, .1 texture{ pigment {color %s}}}\n":
		out << "sphere{<0, " << box_y 
			<< " ," << box_z << ">, .1 texture{ pigment {color %s}}}\n":
		out << "sphere{<" << box_x 
			<< ", 0, 0>, .1 texture{ pigment {color %s}}}\n":
		out << "sphere{<" << box_x << ", 0, " 
			<< box_z << ">, .1 texture{ pigment {color %s}}}\n":
		out << "sphere{<" << box_x << ", " << box_y 
			<< ", 0>, .1 texture{ pigment {color %s}}}\n":
		out << "sphere{<" << box_x << ", " << box_y << ", " 
			<< box_z << ">, .1 texture{ pigment {color %s}}}\n":
	}

	return out.str();
}

