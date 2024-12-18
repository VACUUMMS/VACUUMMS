#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <vacuumms/version.h>

/*
#define VACUUMMS_MAJOR_VERSION 1
#define VACUUMMS_MINOR_VERSION 2
#define VACUUMMS_PATCH_VERSION 0
*/

int main(int argc, char **argv)
{
    // put these in a version file to be read
    int major = VACUUMMS_MAJOR_VERSION;
    int minor = VACUUMMS_MINOR_VERSION;
    int patch = VACUUMMS_PATCH_VERSION;

    if (argc == 1)
    {
        printf("VACUUMMS version %d.%d.%d\n\n", major, minor, patch);
        exit(0);
    }
    else // exec the next command
    {
        char which_command[256];
        sprintf(which_command, "which %s", argv[1]);
        int status = system(which_command);
        printf("status = %d\n", status);       
        if (status == 0) // exists, so run the command
        {
            printf("Not implemented.\n");
            exit(127);

/*
            FILE *fp = popen(which_command, "r");
            char path[256];
            fgets(path, 256, fp);
            pclose(fp);
            
            printf("execv-ing %s:\n", path);
            execv(path, &argv[2]);

            char command[256];
            command[0] = 0;
            sprintf(command, "%s", argv[1]);
            for (int i = 2; i < argc; i++) sprintf(command, "%s %s", command, argv[i]);
            printf("command = >>>%s<<<\n", command);

*/
        }
    }
}
