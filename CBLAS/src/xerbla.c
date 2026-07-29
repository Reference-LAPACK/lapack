#include <stdio.h>

#include "cblas.h"
#include "cblas_f77.h"
#include "cblas_xerbla_internal.h"

/**
 * \brief XERBLA implementation linked in together with the CBLAS library.
 *
 * An error raised beneath a CBLAS wrapper is forwarded to cblas_xerbla()
 * under the CBLAS spelling of the routine name, with the argument number
 * shifted by one to account for the extra layout argument. Errors from
 * direct Fortran calls are reported under the Fortran name instead.
 *
 * \param[in] F77_srname  Routine name reported by the Fortran BLAS, blank
 *                        padded and not NUL terminated.
 * \param[in] vinfo       Pointer to the Fortran argument number.
 * \param[in] len         Hidden Fortran length of \p F77_srname, on
 *                        compilers that pass string lengths at the end of
 *                        the argument list.
 */
void CBLAS_WEAK_SYMBOL F77_xerbla_base(FCHAR F77_srname, void *vinfo
#ifdef BLAS_FORTRAN_STRLEN_END
   ,
   FORTRAN_STRLEN len
#endif
)
{
   extern int CBLAS_CallFromC;

#ifdef BLAS_FORTRAN_STRLEN_END
   const size_t srname_len = len > 0 ? (size_t)len : 0;
#else
   const size_t srname_len = 6;
#endif

#ifdef F77_CHAR
   const char *srname = F2C_STR(F77_srname, srname_len);
#else
   const char *srname = F77_srname;
#endif

   const F77_INT *info = (const F77_INT *)vinfo;
   const CBLAS_INT cblas_info = (CBLAS_INT)*info;

   if (CBLAS_CallFromC) {
      char rout[CBLAS_XERBLA_ROUT_BUFFER_SIZE];
      cblas_xerbla_make_rout(rout, sizeof(rout), srname, srname_len);
      API_SUFFIX(cblas_xerbla)(cblas_info + 1, rout, "");
   } else {
      const size_t display_len =
         cblas_xerbla_trimmed_length(srname, srname_len);

      fprintf(stderr, "Parameter %" CBLAS_IFMT " to routine ", cblas_info);
      if (display_len > 0) {
         fwrite(srname, 1, display_len, stderr);
      }
      fputs(" was incorrect\n", stderr);
   }
}
