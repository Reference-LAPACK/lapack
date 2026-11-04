#include <stdio.h>

#include "cblas.h"
#include "cblas_f77.h"
#include "cblas_xerbla_internal.h"

static blas_fortran_xerbla_callback_t active_callback = NULL;

/* This is a separate function so it can be injected when user wants
   to override cblas_xerbla, but XERBLA can be overridden with
   linker stuff. */
void API_SUFFIX(cblas_F77_xerbla) BLAS_FORTRAN_XERBLA_ARGS
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

/* This means CBLAS calls to these on macOS/Windows then go straight to
the BLAS one, thereby influencing its operation, and defeating the
2-level namespace problem */
#if !defined(_WIN32) && !defined(__APPLE__)
/**
 * \brief XERBLA implementation linked in together with the CBLAS library.
 *
 * An error raised beneath a CBLAS wrapper is forwarded to cblas_xerbla()
 * under the CBLAS spelling of the routine name, with the argument number
 * shifted by one to account for the extra layout argument. Errors from
 * direct Fortran calls are reported under the Fortran name instead.
 * As with the Fortran XERBLA, users can replace the Fortran behaviour by
 * calling SET_BLAS_XERBLA (replicated within CBLAS) with a replacement
 * handler. See separate documentation.
 *
 * \param[in] F77_srname  Routine name reported by the Fortran BLAS, blank
 *                        padded and not NUL terminated.
 * \param[in] info        Pointer to the Fortran argument number.
 * \param[in] len         Hidden Fortran length of \p F77_srname, on
 *                        compilers that pass string lengths at the end of
 *                        the argument list.
 */
void CBLAS_WEAK_SYMBOL F77_xerbla_base BLAS_FORTRAN_XERBLA_ARGS
{
   blas_fortran_xerbla_callback_t handler = active_callback;
   if (!handler) handler = API_SUFFIX(cblas_F77_xerbla);
   handler(F77_srname, info
#ifdef BLAS_FORTRAN_STRLEN_END
      , len
#endif
   );
}

/**
 * \brief Override XERBLA handler within BLAS library.
 *
 * \param[in] fnptr       Pointer to a Fortran-compatible XERBLA replacement.
 *                        Call with NULL to restore default behaviour.
 */
void F77_set_blas_xerbla_base(blas_fortran_xerbla_callback_t cb)
{
    active_callback = cb;
}

/**
 * \brief Return the current XERBLA handler within BLAS library.
 *
 * \param[in] fnptr_ptr   Pointer to a pointer to a Fortran-compatible
 *                        XERBLA replacement, mimicking Fortran behaviour.
 */
void F77_get_blas_xerbla_base(blas_fortran_xerbla_callback_t *cb_ret)
{
    *cb_ret = active_callback;
}
#endif
