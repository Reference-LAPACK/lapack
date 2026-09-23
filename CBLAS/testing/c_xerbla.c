#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "cblas.h"
#include "cblas_test.h"
#include "cblas_xerbla_internal.h"

/**
 * \brief Test-harness stand-in for the CBLAS error handler.
 *
 * Rather than terminating, checks the reported routine name and argument
 * number against the values the calling tester put in cblas_rout and
 * cblas_info, and records the outcome in cblas_ok and cblas_lerr.
 *
 * \param[in] info  CBLAS argument number reported by the caller.
 * \param[in] rout  Routine name, e.g. "cblas_dgemm".
 * \param[in] form  Unused; kept so the signature matches the library's.
 */
void API_SUFFIX(cblas_xerbla)(CBLAS_INT info, const char *rout,
                              const char *form, ...)
{
   extern CBLAS_INT cblas_lerr, cblas_info, cblas_ok;
   extern CBLAS_INT link_xerbla;
   extern CBLAS_INT cblas_xbad;
   extern char *cblas_rout;

   (void)form;

   /* Initially, c__3chke may call this routine with
    * global variable link_xerbla=1, and F77_xerbla will set link_xerbla=0.
    * This is done to fool the linker into loading these subroutines first
    * instead of ones in the CBLAS or the legacy BLAS library.
    */
   if (link_xerbla) return;

   /* The name is checked and reported as the user would see it, suffix and
    * all, while the remapping below keys off the unsuffixed \p rout. The
    * expected name goes through the same normalisation, since the testers
    * spell it without the extended API suffix.
    */
   char reported[CBLAS_XERBLA_ROUT_BUFFER_SIZE];
   cblas_xerbla_apply_api_suffix(reported, sizeof(reported), rout);

   char expected[CBLAS_XERBLA_ROUT_BUFFER_SIZE];
   cblas_xerbla_apply_api_suffix(expected, sizeof(expected), cblas_rout);

   if (cblas_rout != NULL && strcmp(expected, reported) != 0) {
      printf(
         "***** XERBLA WAS CALLED WITH SRNAME = <%s> INSTEAD OF <%s> *******\n",
         reported, expected);
      cblas_ok = FALSE;
      cblas_xbad = 1;
   }

   info = cblas_xerbla_map_info(info, rout, RowMajorStrg);

   if (info != cblas_info) {
      printf("***** XERBLA WAS CALLED WITH INFO = %" CBLAS_IFMT
             " INSTEAD OF %" CBLAS_IFMT " in %s *******\n",
             info, cblas_info, reported);
      cblas_lerr = PASSED;
      cblas_ok = FALSE;
   } else {
      cblas_lerr = FAILED;
   }
}

/**
 * \brief Test-harness XERBLA that redirects Fortran errors to cblas_xerbla().
 *
 * \param[in] F77_srname  Routine name reported by the Fortran BLAS, blank
 *                        padded and not NUL terminated.
 * \param[in] vinfo       Pointer to the Fortran argument number.
 * \param[in] len         Hidden Fortran length of \p F77_srname, on
 *                        compilers that pass string lengths at the end.
 */
void F77_xerbla(FCHAR F77_srname, void *vinfo
#ifdef BLAS_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN len
#endif
)
{
   /* See the comment in API_SUFFIX(cblas_xerbla)() above */
   extern CBLAS_INT link_xerbla;
   if (link_xerbla) {
      link_xerbla = 0;
      return;
   }

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

   char rout[CBLAS_XERBLA_ROUT_BUFFER_SIZE];
   cblas_xerbla_make_rout(rout, sizeof(rout), srname, srname_len);

   /* We increment *info by 1 since the CBLAS interface adds one more
    * argument to all level 2 and 3 routines.
    */
   API_SUFFIX(cblas_xerbla)(cblas_info + 1, rout, "");
}
