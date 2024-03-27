/* gfg2fvi.h */

#include <vacuumms/types.h>

//  These are the routines to call from outside the library
//  They are conditionally declared extern because this header is included by the .cu file 
//  which is compiled by nvcc and which mangles linkage by default.  
//  The libraries thus generated will be linkable from C as well as C++

#ifdef __cplusplus
extern "C" 
#endif
vacuumms_FVI256 *gfgToFVI256(vacuumms_GFG65536 *gfg, float sigma, float epsilon);

#ifdef __cplusplus
extern "C" 
#endif
vacuumms_EnergyArray256 *GFGToEnergyArray256(vacuumms_GFG65536 *gfg, float sigma, float epsilon);




#ifdef __cplusplus
extern "C" 
#endif
vacuumms_EnergyArray256 *GFGToRepulsion256(vacuumms_GFG65536 *gfg, float sigma, float epsilon);

#ifdef __cplusplus
extern "C" 
#endif
vacuumms_EnergyArray256 *GFGToRepulsion256_612(vacuumms_GFG65536 *gfg, float sigma, float epsilon);



#ifdef __cplusplus
extern "C" 
#endif
vacuumms_EnergyArray512 *GFGToRepulsion512(vacuumms_GFG65536 *gfg, float sigma, float epsilon);

#ifdef __cplusplus
extern "C" 
#endif
vacuumms_EnergyArray512 *GFGToRepulsion512_612(vacuumms_GFG65536 *gfg, float sigma, float epsilon);



/*
#ifdef __cplusplus
extern "C" 
#endif
vacuumms_EnergyArray1024 *GFGToEnergyArray1024_69(vacuumms_GFG65536 *gfg, float sigma, float epsilon);

#ifdef __cplusplus
extern "C" 
#endif
vacuumms_EnergyArray1024 *GFGToRepulsion1024_69(vacuumms_GFG65536 *gfg, float sigma, float epsilon);
*/

#ifdef __cplusplus
extern "C" 
#endif
vacuumms_EnergyArray1024 *GFGToEnergyArray1024_612(vacuumms_GFG65536 *gfg, float sigma, float epsilon);

#ifdef __cplusplus
extern "C" 
#endif
vacuumms_EnergyArray1024 *GFGToRepulsion1024_612(vacuumms_GFG65536 *gfg, float sigma, float epsilon);

struct Chunk
{
  float energy[256][1024][1024];
};

typedef struct Chunk vacuumms_Chunk;

#endif 

