#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "cblas.h"
#include "cblas_f77.h"

/**
 * \brief Return the API version of BLAS being wrapped.
 *
 * \param[out] major BLAS API major version.
 * \param[out] minor BLAS API minor version.
 * \param[out] patch BLAS API patch number.
 */
void API_SUFFIX(cblas_iblaver)(CBLAS_INT *major,
                               CBLAS_INT *minor,
                               CBLAS_INT *patch)
{
  F77_INT F77_major = 0, F77_minor = 0, F77_patch = 0;
  if (major == NULL) {
    API_SUFFIX(cblas_xerbla)(1, "cblas_iblaver",
       "Illegal major setting, NULL\n");
    return;
  }
  if (minor == NULL) {
    API_SUFFIX(cblas_xerbla)(2, "cblas_iblaver",
       "Illegal minor setting, NULL\n");
    return;
  }
  if (patch == NULL) {
    API_SUFFIX(cblas_xerbla)(3, "cblas_iblaver",
       "Illegal patch setting, NULL\n");
    return;
  }
  F77_iblaver(&F77_major, &F77_minor, &F77_patch);
  *major = (CBLAS_INT)F77_major;
  *minor = (CBLAS_INT)F77_minor;
  *patch = (CBLAS_INT)F77_patch;
}
