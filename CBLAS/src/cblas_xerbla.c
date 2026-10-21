#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "cblas.h"
#include "cblas_f77.h"
#include "cblas_xerbla_internal.h"

static cblas_xerbla_callback_t active_callback = NULL;

/**
 * \brief CBLAS error handler: report an invalid argument and terminate.
 *
 * Remaps \p info to the CBLAS argument position for row-major calls, writes
 * the diagnostic to stderr and exits; it does not return to its caller.
 * Users can replace this behaviour by calling cblas_set_xerbla with a
 * replacement handler. See separate documentation.
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
   info = cblas_xerbla_map_info(info, rout, RowMajorStrg);
   va_list argptr;
   va_start(argptr, form);
   if (active_callback) {
      active_callback(info, rout, form, argptr);
   } else {
      if (info) {
         char reported[CBLAS_XERBLA_ROUT_BUFFER_SIZE];
         cblas_xerbla_apply_api_suffix(reported, sizeof(reported), rout);
         fprintf(stderr, "Parameter %" CBLAS_IFMT " to routine %s was incorrect\n",
                 info, reported);
      }
      vfprintf(stderr, form, argptr);
   }
   va_end(argptr);
   if (!active_callback)
      exit(-1);
}

/**
 * \brief Override cblas_xerbla handler.
 *
 * Returns the old CBLAS handler.
 *
 * The new handler takes a va_list argument instead of the "...". The
 * info parameter will already have been mapped from its Fortran value
 * so you don't need to do that in your handler.
 * Users should usually also use LAPACKE_set_xerbla since the error
 * handlers may be separate due to the vagaries of linkage.
 *
 * If users wish to chain their handlers, calling other
 * previously-registered ones after their own, they must be careful
 * not to create loops.
 *
 * \param[in] cblas_xerbla_ptr Pointer to a cblas_xerbla replacement.
 *                             Call with NULL to restore default behaviour.
 * \param[in] F77_xerbla_ptr   Pointer to a Fortran-compatible XERBLA
 *                             replacement. Call with NULL if no control
 *                             over Fortran-originated errors is desired;
 *                             the override mechanism will do the right thing.
 */
cblas_xerbla_callback_t API_SUFFIX(cblas_set_xerbla)(cblas_xerbla_callback_t cb, blas_fortran_xerbla_callback_t F77_cb)
{
   cblas_xerbla_callback_t old_cb = active_callback;
   active_callback = cb;
   F77_set_blas_xerbla_base(
      cb ? (F77_cb ? F77_cb : API_SUFFIX(cblas_F77_xerbla)) /* set up new handler */
         : NULL /* restore default */
   );
   return old_cb;
}

/**
 * \brief Returns the current cblas_xerbla handler.
 */
cblas_xerbla_callback_t API_SUFFIX(cblas_get_xerbla)(void)
{
   return active_callback;
}
