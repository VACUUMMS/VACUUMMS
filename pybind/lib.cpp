#include "lib.h"

namespace lib
{
int f()
{
    std::cout << "f()\n";
    return 0;
}

int g(int a)
{
    std::cout << "g(int)\n";
    std::cout << a << "\n";
    return 1;
}
}
