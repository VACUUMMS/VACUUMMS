#include <assert.h>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv)
{
    int val = (int)strtol(argv[1], NULL, 10);
    printf("got val = %d\n", val);
    assert (val==0);
    return 0;
}
