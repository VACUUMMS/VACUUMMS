/*********************** rng.cc *********************************/
#include <math.h>
#include <time.h>
#include <stdio.h>

#include <vacuumms/rng.hh>

/*
#define thirtyone 2147483648.0
#ifndef DEFAULT_SEED
#define DEFAULT_SEED 12345
#endif

int ra[256], nd;

void initializeRandomNumberGeneratorTo(int seed)
{
  randinit(seed);
}

void initializeRandomNumberGenerator()
{
  randinit(DEFAULT_SEED);
}

void randinit(int seed)
{
  vacuumms_float dubee, e;
  int trune;
  long int i;
  dubee=-1.0+1.0/thirtyone;
  e = seed/thirtyone;
  for (i=1; i<10000; i++)
  {
    e*=16807;
    trune=e;
    e+=dubee*trune;
    if (e>=1) e+= dubee;
  }
  for (nd=0;nd<=255; nd++)
  {
    e=16807*e;
    trune=e;
    e+=dubee*trune;
    if (e>=1) e+=dubee;
    ra[nd]=(thirtyone)*e;
  }
  for (i=0; i<=10000; i++)
  {
    nd=(nd+1)&255;
    ra[nd]=(ra[(nd-103)&255])^(ra[(nd-250)&255]);
  }
}

vacuumms_float RND()
{
  return rnd();
}

vacuumms_float rnd()
{
  nd=(nd+1)&255;
  ra[nd]=(ra[(nd-103)&255]^ra[(nd-250)&255]);
  return ra[nd]/thirtyone;
} 

int random_int(int max)
{
  return floor(rnd() * max);
}

vacuumms_float random_vacuumms_float(vacuumms_float max)
{
  return (rnd() * max);
}

// randomly initializes rng 
int randomize()
{
  time_t now;
  struct tm *p_tyme;
  struct tm tyme;
  int rng_seed;

  now = time(NULL);
  p_tyme = gmtime(&now);
  tyme = *p_tyme;
  rng_seed = tyme.tm_sec * tyme.tm_min * tyme.tm_hour;
  initializeRandomNumberGeneratorTo(rng_seed);
  return rng_seed;
}

// returns random seed value 
int getRandomSeed()
{
  time_t now;
  struct tm *p_tyme;
  struct tm tyme;
  int rng_seed;

  now = time(NULL);
  p_tyme = gmtime(&now);
  tyme = *p_tyme;
  rng_seed = tyme.tm_sec * tyme.tm_min * tyme.tm_hour;
//printf("rng_seed = %d", rng_seed);
  return rng_seed;
}
*/ 

// New one starts here

#include <cstdint>

void MersenneTwister::twist() 
{
    for (int i = 0; i < N; i++) 
    {
        uint32_t y = (mt[i] & UPPER_MASK) | (mt[(i + 1) % N] & LOWER_MASK);
        mt[i] = mt[(i + M) % N] ^ (y >> 1) ^ ((y & 1) ? 0x9908b0dfU : 0);
    }
    index = 0;
}

MersenneTwister::MersenneTwister(uint32_t seed) : index(N) 
{
    mt[0] = seed;
    for (int i = 1; i < N; i++)
        mt[i] = 0x6c078965U * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i;
}

uint32_t MersenneTwister::next() 
{
    if (index >= N) twist();
    uint32_t y = mt[index++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680U;
    y ^= (y << 15) & 0xefc60000U;
    y ^= (y >> 18);
    return y;
}

vacuumms_float MersenneTwister::next_float() 
{
    return next() / 4294967296.0; 
} // [0,1)


