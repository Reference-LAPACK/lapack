/*
 * cblas_zrotg.c
 *
 * The program is a C interface to zrotg.
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
void cblas_zrotg(void *a, void *b, double *c, void *s)
{
   F77_zrotg(a,b,c,s);
}

