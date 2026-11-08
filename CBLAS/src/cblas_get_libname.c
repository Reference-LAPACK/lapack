#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "cblas.h"
#include "cblas_f77.h"

/**
 * \brief Get the BLAS library self-identifying string.
 *
 * \param[out] buffer  Buffer to write the string to. Will always
                       be NUL-terminated if buf_len > 0.
 * \param[out] out_len Number of bytes written to buffer, including any NUL.
 * \param[in]  buf_len Length of the output buffer.
 */
void API_SUFFIX(cblas_get_libname)(char *buffer, CBLAS_INT *out_len, CBLAS_INT buf_len)
{
  F77_INT F77_out_len = 0;
  if (buffer == NULL) {
      API_SUFFIX(cblas_xerbla)(1, "cblas_get_libname",
        "Illegal buffer setting, NULL\n");
      return;
  }
  if (out_len == NULL) {
      API_SUFFIX(cblas_xerbla)(2, "cblas_get_libname",
        "Illegal out_len setting, NULL\n");
      return;
  }
  if (buf_len < 0) {
    API_SUFFIX(cblas_xerbla)(3, "cblas_get_libname",
      "Illegal buf_len setting, %" CBLAS_IFMT "\n", buf_len);
    return;
  }
  F77_get_blas_libname(buffer, &F77_out_len, (FORTRAN_STRLEN)buf_len);
  *out_len = (CBLAS_INT)F77_out_len;
}
