/* include/vacuumms/rng.hh */

#pragma once
#include <vacuumms/types.h>

/*
void initializeRandomNumberGeneratorTo(int seed);
void initializeRandomNumberGenerator();
void randinit(int seed);
vacuumms_float RND();
vacuumms_float rnd();
int randomize();
*/

#include <cstdint>

class MersenneTwister 
{
    private:
        static constexpr int N = 624, M = 397;
        static constexpr uint32_t UPPER_MASK = 0x80000000U, LOWER_MASK = 0x7fffffffU;
        uint32_t mt[N];
        int index;

        void twist();
/*
        {
            for (int i = 0; i < N; i++) {
                uint32_t y = (mt[i] & UPPER_MASK) | (mt[(i + 1) % N] & LOWER_MASK);
                mt[i] = mt[(i + M) % N] ^ (y >> 1) ^ ((y & 1) ? 0x9908b0dfU : 0);
            }
            index = 0;
        }
*/

    public:

        MersenneTwister(uint32_t seed) ;
        MersenneTwister(): MersenneTwister(5489) {};
/*
        {
            mt[0] = seed;
            for (int i = 1; i < N; i++)
                mt[i] = 0x6c078965U * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i;
        }
*/
        uint32_t next();
/*
        {
            if (index >= N) twist();
            uint32_t y = mt[index++];
            y ^= (y >> 11);
            y ^= (y << 7) & 0x9d2c5680U;
            y ^= (y << 15) & 0xefc60000U;
            y ^= (y >> 18);
            return y;
        }
*/
        vacuumms_float next_float(); // { return next() / 4294967296.0; } // [0,1)
};


