/* config_parser.h */

#ifndef FTW_CONFIG_PARSER_H
#define FTW_CONFIG_PARSER_H

//  A set of library routines to parse input file

//  IN:  FILE* pointer to an input stream
//  OUT:  configuration, and metadata:  number of lines

#include <stdlib.h>
#include <stdio.h>
#include <vacuumms/param.h>
#include <string.h>

#include <vacuumms/limits.h>

#include <vacuumms/types.h>

// return value is pointer to configuration

#ifdef __cplusplus
extern "C"
#endif
vacuumms_GFG65536 *readGFG65536(FILE *instream);

#ifdef __cplusplus
extern "C"
#endif
vacuumms_GFG65536 *replicateGFG65536(vacuumms_GFG65536 *in);

#endif

