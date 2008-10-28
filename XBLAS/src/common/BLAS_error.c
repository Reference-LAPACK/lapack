#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "blas_extended.h"
#include "f2c-bridge.h"

#if defined(CONFIG_USE_XERBLA)
#define xerbla_array FC_FUNC_(xerbla_array, XERBLA_ARRAY)
extern void xerbla_array(const char *, const int *, const int *);
#endif

void BLAS_error(const char *rname, int iflag, int ival, char *form, ...)
/*
 * Argument
 * ========
 * rname     (input) routine name
 *
 * iflag     (input) a negative value indicates that parameter number -IFLAG
 *                   caused the error; a nonnegative value is an
 *                   implementation-specific error code.
 *
 * ival      (input) the value of parameter number -IFLAG.
 */
{
#if !defined(CONFIG_USE_XERBLA)
  {
    va_list argptr;
    va_start(argptr, form);
    fprintf(stderr, "Error #%d from routine %s:\n", iflag, rname);
    if (form)
      vfprintf(stderr, form, argptr);
    else if (iflag < 0)
      fprintf(stderr,
	      "  Parameter number %d to routine %s had the illegal value %d\n",
	      -iflag, rname, ival);
    else
      fprintf(stderr, "  Unknown error code %d from routine %s\n",
	      iflag, rname);
    exit(iflag);
  }
#else
  {
    int ln, argno;
    ln = strlen(rname);
    argno = -iflag;
    xerbla_array(rname, &ln, &argno);
  }
#endif
}
