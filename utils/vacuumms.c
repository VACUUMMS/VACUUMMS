#include <stdio.h>
#include <vacuumms/version.h>

int main()
{
    // put these in a version file to be read
    int major = VACUUMMS_MAJOR_VERSION;
    int minor = VACUUMMS_MINOR_VERSION;
    int patch = VACUUMMS_PATCH_VERSION;
    printf("VACUUMMS version %d.%d.%d\n\n", major, minor, patch);
}
