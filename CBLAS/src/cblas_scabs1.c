/*
 * cblas_scabs1.c
 *
 * The program is a C interface to scabs1.
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
float API_SUFFIX(cblas_scabs1)(const void *c)
{
   float cabs1 = 0.0;
   F77_scabs1_sub(c, &cabs1);
   return cabs1; 
}

