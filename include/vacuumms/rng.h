/* rng.h */

#ifndef FTW_RNG_INCLUDE
#define FTW_RNG_INCLUDE

#ifdef __cplusplus
extern "C" {
#endif

void initializeRandomNumberGeneratorTo(int seed);
void initializeRandomNumberGenerator();
void randinit(int seed);
double RND();
double rnd();
int randomize();
double rnd_double();
int rnd_int();
int getRandomSeed();

#ifdef __cplusplus
}
#endif

#endif

