#include <stdio.h>
#include <math.h>

// Generate the camera_light-%03d.pov files which will be included in the templates
int main()
{
    char filename[256], camera[256], light[256];
    float r_camera = 0.1;
    float r_light = 176.88;

    for (int i=0; i<360; i++)
    {
        sprintf(filename, "pov/camera_light-%03d.pov", i);
        printf("%s\n", filename);
        float theta = (float)i * M_PI/180.0;

//        float x_look_at = 1.4;
//        float y_look_at = 1.4;
//        float z_look_at = 1.4;

        float x_look_at = 11.0;
        float y_look_at = 11.0;
        float z_look_at = 11.0;

        // camera
        sprintf(camera, "camera{location<%f,%f,%f> look_at <%f,%f,%f> right 1.0}\n", 
                         x_look_at + r_camera * cos(theta), 
                         y_look_at, 
                         z_look_at + r_camera * sin(theta),
                         x_look_at,
                         y_look_at,
                         z_look_at);
/*
        sprintf(camera, "camera{location<%f,1.4,%f> look_at <1.4,1.4,1.4> right 1.0}\n", 
                         1.4 + r_camera * cos(theta), 
                         1.4 + r_camera * sin(theta));
*/
        
        // light
        sprintf(light, "light_source{<%f,120.000000,%f> color White }\n", 
                        r_light * cos(theta), 
                        r_light * sin(theta));

        FILE *file = fopen(filename, "w");
        fprintf(file, "%s\n", camera);
        fprintf(file, "%s\n", light);
        fclose(file);
    }

    return 0;
}
