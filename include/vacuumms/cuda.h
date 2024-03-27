/* cuda.h */

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>

#include <vacuumms/types.h>

#ifndef FTW_CUDA_H
#define FTW_CUDA_H

//  These are the routines to call from outside the library
//  They are conditionally declared extern because this header is included by the .cu file 
//  which is compiled by nvcc and which mangles linkage by default.  
//  The libraries thus generated will be linkable from C as well as C++

#ifdef __cplusplus
extern "C" 
#endif
void vacuumms_printDeviceInfo(int);

#ifdef __cplusplus
extern "C" 
#endif
int vacuumms_getNumberOfDevices();

#endif

