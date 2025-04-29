/* prng.hh */

/****************************************************************************/
/*                                                                          */
/*                          Parallel RNG		                            */
/*                                                                          */
/*                        (C) 2012 Frank T Willmore                         */
/*                                                                          */
/*    Contains parallel implementation of Mersenne Twister                  */
/*                                                                          */
/*                                                                          */
/*    The Texas Advanced Computing Center                                   */
/*    The National Science Foundation/NSF-Teragrid                          */
/*                                                                          */
/*    correspondence to:  frankwillmore@gmail.com                           */
/*                                                                          */
/****************************************************************************/

#define VACUUMMS_MERSENNE_NN 312
#define VACUUMMS_MERSENNE_MM 156
#define VACUUMMS_MERSENNE_MATRIX_A 0xB5026F5AA96619E9ULL
#define VACUUMMS_MERSENNE_UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define VACUUMMS_MERSENNE_LM 0x7FFFFFFFULL /* Least significant 31 bits */

struct MersenneTwister
{
  unsigned long long mt[VACUUMMS_MERSENNE_NN];
  int mti;
};

/* prototypes */

void MersenneInitialize(struct MersenneTwister* MT, int seed);
double prnd(struct MersenneTwister* MT);
