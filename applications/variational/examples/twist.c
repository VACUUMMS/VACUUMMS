#include <stdio.h>
#include <math.h>

// Generate the camera_light-%03d.pov files which will be included in the templates
void main()
{
    char filename[256], camera[256], light[256];
//    float r_camera = 4.47;
    float r_camera = 0.1;
    float r_light = 176.88;

    for (int i=0; i<360; i++)
    {
        sprintf(filename, "pov/camera_light-%03d.pov", i);
        printf("%s\n", filename);
        FILE *file = fopen(filename, "w");
        float theta = (float)i * M_PI/180.0;

        // camera
//        sprintf(camera, "camera{location<%f,4.0,%f> look_at <1,1,1> right 1.0}\n", 1.0 + r_camera * cos(theta), 1.0 + r_camera * sin(theta));
        sprintf(camera, "camera{location<%f,1.4,%f> look_at <1.4,1.4,1.4> right 1.0}\n", 1.4 + r_camera * cos(theta), 1.4 + r_camera * sin(theta));
        
        // light
        sprintf(light, "light_source{<%f,120.000000,%f> color White }\n", r_light * cos(theta), r_light * sin(theta));

        fprintf(file, "%s\n", camera);
        fprintf(file, "%s\n", light);
        fclose(file);
    }
}
