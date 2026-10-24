#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "cblas.h"
#include "cblas_f77.h"

/**
 * \brief Return CBLAS version.
 *
 * \param[out] major BLAS major version.
 * \param[out] minor BLAS minor version.
 * \param[out] patch BLAS patch number.
 */
void API_SUFFIX(cblas_iblaver)(CBLAS_INT *major,
                               CBLAS_INT *minor,
                               CBLAS_INT *patch)
{
   F77_iblaver(major, minor, patch);
}
