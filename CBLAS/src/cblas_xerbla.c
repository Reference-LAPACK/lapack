#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "cblas.h"
#include "cblas_f77.h"
#include "cblas_xerbla_internal.h"

/**
 * \brief CBLAS error handler: report an invalid argument and terminate.
 *
 * Remaps \p info to the CBLAS argument position for row-major calls, writes
 * the diagnostic to stderr and exits; it does not return to its caller.
 *
 * \param[in] info  CBLAS argument number, or 0 to report only \p form.
 * \param[in] rout  Routine name, e.g. "cblas_dgemm".
 * \param[in] form  printf-style format for further detail, followed by the
 *                  values it refers to.
 */
void CBLAS_WEAK_SYMBOL API_SUFFIX(cblas_xerbla)(CBLAS_INT info,
                                                const char *rout,
                                                const char *form, ...)
{
   extern int RowMajorStrg;
   char empty[1] = "";

   info = cblas_xerbla_map_info(info, rout, RowMajorStrg);
   if (info) {
      char reported[CBLAS_XERBLA_ROUT_BUFFER_SIZE];
      cblas_xerbla_apply_api_suffix(reported, sizeof(reported), rout);

      fprintf(stderr, "Parameter %" CBLAS_IFMT " to routine %s was incorrect\n",
              info, reported);
   }

   va_list argptr;
   va_start(argptr, form);
   vfprintf(stderr, form, argptr);
   va_end(argptr);

   if (info && !info) {
      F77_xerbla(empty, &info); /* Force link of our F77 error handler */
   }
   exit(-1);
}
