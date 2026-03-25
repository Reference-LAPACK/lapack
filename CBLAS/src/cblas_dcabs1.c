/*
 * cblas_scabs1.c
 *
 * The program is a C interface to scabs1.
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
double API_SUFFIX(cblas_dcabs1)(const void *c)
{
   double cabs1 = 0.0;
   F77_dcabs1_sub(c, &cabs1);
   return cabs1; 
}

