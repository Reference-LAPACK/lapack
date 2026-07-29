#ifndef CBLAS_XERBLA_INTERNAL_H
#define CBLAS_XERBLA_INTERNAL_H

#include <ctype.h>
#include <stddef.h>
#include <string.h>

#include "cblas.h"

// Keep a bit of headroom for possible future CBLAS routines
#define CBLAS_XERBLA_MAX_ROUTINE_NAME 32u
#define CBLAS_XERBLA_API_PREFIX "cblas_"
#define CBLAS_XERBLA_API64_SUFFIX "_64"
#ifdef CBLAS_API64
#define CBLAS_XERBLA_API_SUFFIX CBLAS_XERBLA_API64_SUFFIX
#else
#define CBLAS_XERBLA_API_SUFFIX ""
#endif

#define CBLAS_XERBLA_ROUT_BUFFER_SIZE                                          \
   ((sizeof(CBLAS_XERBLA_API_PREFIX) - 1) + CBLAS_XERBLA_MAX_ROUTINE_NAME +    \
    (sizeof(CBLAS_XERBLA_API_SUFFIX) - 1) + 1)

/**
 * \brief Length of a Fortran-style name with trailing blanks removed.
 *
 * \param[in] name      Character data, not necessarily NUL terminated.
 * \param[in] name_len  Number of characters available in \p name.
 *
 * \return Length up to the first NUL, excluding trailing blanks, or 0 if
 *         \p name is NULL.
 */
static inline size_t cblas_xerbla_trimmed_length(const char *name,
                                                 const size_t name_len)
{
   if (name == NULL) return 0;

   size_t actual_len = 0;
   while (actual_len < name_len && name[actual_len] != '\0') {
      actual_len++;
   }
   while (actual_len > 0 && name[actual_len - 1] == ' ') {
      actual_len--;
   }

   return actual_len;
}

/**
 * \brief Build the "cblas_"-prefixed routine name used in error messages.
 *
 * Lowercases \p name and drops trailing blanks and, in 64-bit API builds, a
 * trailing "_64". The result is always NUL terminated, and is truncated
 * rather than allowed to overflow \p rout.
 *
 * \param[out] rout       Destination buffer, normally of size
 *                        CBLAS_XERBLA_ROUT_BUFFER_SIZE.
 * \param[in]  rout_size  Size of \p rout in bytes.
 * \param[in]  name       Fortran routine name, not necessarily NUL
 *                        terminated.
 * \param[in]  name_len   Number of characters available in \p name.
 */
static inline void cblas_xerbla_make_rout(char *rout, const size_t rout_size,
                                          const char *name, size_t name_len)
{
   if (rout_size == 0) return;

   name_len = cblas_xerbla_trimmed_length(name, name_len);

   static const char suffix[] = CBLAS_XERBLA_API_SUFFIX;
   const size_t suffix_len = sizeof(suffix) - 1;
   if (suffix_len > 0 && name_len >= suffix_len &&
       strncmp(name + name_len - suffix_len, suffix, suffix_len) == 0) {
      name_len -= suffix_len;
   }

   if (name_len > CBLAS_XERBLA_MAX_ROUTINE_NAME) {
      name_len = CBLAS_XERBLA_MAX_ROUTINE_NAME;
   }

   static const char prefix[] = CBLAS_XERBLA_API_PREFIX;
   const size_t prefix_len = sizeof(prefix) - 1;
   size_t rout_len = 0;
   for (size_t i = 0; i < prefix_len && rout_len + 1 < rout_size; i++) {
      rout[rout_len++] = prefix[i];
   }
   for (size_t i = 0; i < name_len && rout_len + 1 < rout_size; i++) {
      rout[rout_len++] = (char)tolower((unsigned char)name[i]);
   }

   rout[rout_len] = '\0';
}

/**
 * \brief Copy a routine name as it should be reported to the user.
 *
 * Appends the extended API suffix, so that a 64-bit build names
 * cblas_dgemm_64() rather than cblas_dgemm() in its diagnostics. A name that
 * already carries the suffix is copied unchanged, which keeps the call
 * idempotent whatever the caller passes. The result is always NUL terminated
 * and is truncated rather than allowed to overflow \p rout.
 *
 * \param[out] rout       Destination buffer, normally of size
 *                        CBLAS_XERBLA_ROUT_BUFFER_SIZE.
 * \param[in]  rout_size  Size of \p rout in bytes.
 * \param[in]  name       Routine name, e.g. "cblas_dgemm", or NULL.
 */
static inline void cblas_xerbla_apply_api_suffix(char *rout,
                                                 const size_t rout_size,
                                                 const char *name)
{
   if (rout_size == 0) return;

   if (name == NULL) {
      rout[0] = '\0';
      return;
   }

   static const char suffix[] = CBLAS_XERBLA_API_SUFFIX;
   size_t suffix_len = sizeof(suffix) - 1;
   const size_t name_len = strlen(name);
   if (suffix_len > 0 && name_len >= suffix_len &&
       strcmp(name + name_len - suffix_len, suffix) == 0) {
      suffix_len = 0;
   }

   size_t rout_len = 0;
   for (size_t i = 0; i < name_len && rout_len + 1 < rout_size; i++) {
      rout[rout_len++] = name[i];
   }
   for (size_t i = 0; i < suffix_len && rout_len + 1 < rout_size; i++) {
      rout[rout_len++] = suffix[i];
   }

   rout[rout_len] = '\0';
}

/**
 * \brief Reduce a routine name to its bare operation.
 *
 * Skips the "cblas_" prefix and the precision character, so that
 * "cblas_dgemm" yields "gemm".
 *
 * \param[in] rout  CBLAS routine name, or NULL.
 *
 * \return Pointer into \p rout past the prefix and precision character, or
 *         NULL if \p rout is NULL.
 */
static inline const char *cblas_xerbla_operation(const char *rout)
{
   if (rout == NULL) return NULL;

   static const char prefix[] = CBLAS_XERBLA_API_PREFIX;
   const size_t prefix_len = sizeof(prefix) - 1;
   if (strncmp(rout, prefix, prefix_len) == 0) {
      rout += prefix_len;
   }
   if ((rout[0] == 's' || rout[0] == 'd' || rout[0] == 'c' || rout[0] == 'z') &&
       rout[1] != '\0') {
      rout++;
   }

   return rout;
}

/**
 * \brief Test an operation name for equality.
 *
 * The comparison is exact, so "gemm" does not also match "gemmtr". A trailing
 * extended API suffix is tolerated, so that a name arriving already suffixed
 * still selects the right remapping rather than silently selecting none.
 *
 * \param[in] operation  Result of cblas_xerbla_operation(), or NULL.
 * \param[in] expected   Operation name to match.
 *
 * \return Nonzero when \p operation equals \p expected, ignoring any trailing
 *         extended API suffix.
 */
static inline int cblas_xerbla_operation_is(const char *operation,
                                            const char *expected)
{
   if (operation == NULL) return 0;

   const size_t expected_len = strlen(expected);
   if (strncmp(operation, expected, expected_len) != 0) return 0;

   return operation[expected_len] == '\0' ||
          strcmp(operation + expected_len, CBLAS_XERBLA_API64_SUFFIX) == 0;
}

/**
 * \brief Map a Fortran argument number onto its CBLAS position.
 *
 * Row-major calls reach the Fortran BLAS with arguments swapped or
 * transposed, so the number XERBLA reports is not that of the CBLAS
 * argument actually at fault. Column-major calls, and operations needing no
 * adjustment, return \p info unchanged.
 *
 * \param[in] info       Argument number reported by the Fortran BLAS.
 * \param[in] rout       CBLAS routine name, e.g. "cblas_dgemm".
 * \param[in] row_major  Nonzero if the call used CblasRowMajor.
 *
 * \return The corresponding CBLAS argument number.
 */
static inline CBLAS_INT cblas_xerbla_map_info(CBLAS_INT info, const char *rout,
                                              const int row_major)
{
   if (!row_major) return info;

   const char *const operation = cblas_xerbla_operation(rout);
   if (cblas_xerbla_operation_is(operation, "gemmtr")) {

      if (info == 11) info = 9;
      else if (info == 9) info = 11;

   } else if (cblas_xerbla_operation_is(operation, "gemm")) {

      if (info == 5) info = 4;
      else if (info == 4) info = 5;
      else if (info == 11) info = 9;
      else if (info == 9) info = 11;

   } else if (cblas_xerbla_operation_is(operation, "symm") ||
              cblas_xerbla_operation_is(operation, "hemm") ||
              cblas_xerbla_operation_is(operation, "skewsymm")) {

      if (info == 5) info = 4;
      else if (info == 4) info = 5;

   } else if (cblas_xerbla_operation_is(operation, "trmm") ||
              cblas_xerbla_operation_is(operation, "trsm")) {

      if (info == 7) info = 6;
      else if (info == 6) info = 7;

   } else if (cblas_xerbla_operation_is(operation, "gemv")) {

      if (info == 4) info = 3;
      else if (info == 3) info = 4;

   } else if (cblas_xerbla_operation_is(operation, "gbmv")) {

      if (info == 4) info = 3;
      else if (info == 3) info = 4;
      else if (info == 6) info = 5;
      else if (info == 5) info = 6;

   } else if (cblas_xerbla_operation_is(operation, "ger") ||
              cblas_xerbla_operation_is(operation, "geru") ||
              cblas_xerbla_operation_is(operation, "gerc")) {

      if (info == 3) info = 2;
      else if (info == 2) info = 3;
      else if (info == 8) info = 6;
      else if (info == 6) info = 8;

   } else if (cblas_xerbla_operation_is(operation, "her2") ||
              cblas_xerbla_operation_is(operation, "hpr2")) {

      if (info == 8) info = 6;
      else if (info == 6) info = 8;
   }

   return info;
}

#endif // CBLAS_XERBLA_INTERNAL_H
