/*
 * cblas_crotg.c
 *
 * The program is a C interface to crotg.
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
void API_SUFFIX(cblas_crotg)(void *a, void *b, float *c, void *s)
{
   F77_crotg(a,b,c,s);
}

