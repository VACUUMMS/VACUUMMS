/* include/vacuumms/rng.hh */

#pragma once
#include <vacuumms/types.h>


#include <cstdint>

class MersenneTwister 
{
    private:
        static constexpr int N = 624, M = 397;
        static constexpr uint32_t UPPER_MASK = 0x80000000U, LOWER_MASK = 0x7fffffffU;
        uint32_t mt[N];
        int index;

        void twist();

    public:

        MersenneTwister(uint32_t seed);
        MersenneTwister(): MersenneTwister(5489) {};

        uint32_t next();
        vacuumms_float next_float();
};


