/******************************************************************************
 * Shared helpers for the LAPACKE test wrapper libraries.
 * See lapacke_test_wrappers.h.
 ******************************************************************************/

#include <stdio.h>
#include <string.h>

#include "lapacke_test_wrappers.h"

/* The testing XERBLA provided by the test programs (non-halting). */
#define fortran_xerbla LAPACK_GLOBAL_SUFFIX(xerbla, XERBLA)
void fortran_xerbla(const char *srname, const lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN srname_len
#endif
);

/**
 * \brief Report an argument error through the (interposed, non-halting)
 * testing XERBLA, so that errors LAPACKE detects itself rather than the
 * underlying Fortran routine also reach the CHKXER machinery.
 *
 * \param[in] srname Fortran routine name reported to XERBLA.
 * \param[in] info   The info value, in Fortran numbering.
 */
static void lapacke_test_xerbla(const char *srname, lapack_int info)
{
    lapack_int pos = -info;
    fortran_xerbla(srname, &pos
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   (FORTRAN_STRLEN)strlen(srname)
#endif
    );
}

/**
 * \brief Map a LAPACKE info return value back to Fortran numbering and report
 * argument errors through the testing XERBLA.
 *
 * LAPACKE returns negative info shifted by one because of the leading
 * matrix_layout argument. This undoes that shift. Memory errors and an
 * unexpected matrix_layout error are passed through with a diagnostic on
 * stderr.
 *
 * \param[in] srname Fortran routine name reported to XERBLA.
 * \param[in] ret    Return value of the LAPACKE call.
 * \return The info value in Fortran numbering.
 */
lapack_int lapacke_test_info(const char *srname, lapack_int ret)
{
    if (ret == LAPACK_WORK_MEMORY_ERROR ||
        ret == LAPACK_TRANSPOSE_MEMORY_ERROR) {
        fprintf(stderr,
                "LAPACKE test wrapper %s: memory allocation error"
                " (info = %d)\n",
                srname, (int)ret);
        return ret;
    }
    if (ret < 0) {
        if (ret == -1) {
            /* Invalid matrix_layout cannot originate from these wrappers */
            fprintf(stderr,
                    "LAPACKE test wrapper %s: unexpected"
                    " matrix_layout error\n",
                    srname);
            return ret;
        }
        ret = ret + 1; /* undo the matrix_layout shift */
        lapacke_test_xerbla(srname, ret);
    }
    return ret;
}

/**
 * \brief lapacke_test_info for the routines LAPACKE declares
 * layout-independent.
 *
 * xLAMCH, xLARNV, xLACGV, xROT and the tridiagonal routines take no
 * matrix_layout argument, so their info needs no shift and argument one is
 * reportable like any other.
 *
 * \param[in] srname Fortran routine name reported to XERBLA.
 * \param[in] ret    Return value of the LAPACKE call.
 * \return The info value in Fortran numbering.
 */
lapack_int lapacke_test_info_unshifted(const char *srname, lapack_int ret)
{
    if (ret == LAPACK_WORK_MEMORY_ERROR ||
        ret == LAPACK_TRANSPOSE_MEMORY_ERROR) {
        fprintf(stderr,
                "LAPACKE test wrapper %s: memory allocation error"
                " (info = %d)\n",
                srname, (int)ret);
        return ret;
    }
    if (ret < 0) {
        lapacke_test_xerbla(srname, ret);
    }
    return ret;
}

#if LAPACKE_TEST_LAYOUT == LAPACK_ROW_MAJOR

/* Row-major shadow buffers.
 *
 * One cm_to_rm/rm_to_cm pair per conventional-storage matrix type, mirroring
 * LAPACKE's own <x><type>_trans usage: cm_to_rm allocates a row-major copy of
 * the caller's column-major matrix (NULL when the allocation fails) and
 * rm_to_cm copies the result back. Only the elements the storage scheme
 * references are copied, so the untouched parts of a triangle, a band or a
 * packed array keep whatever the caller left there.
 *
 * Wrappers seed every shadow, including the matrices LAPACKE only writes: a
 * routine that writes just one triangle would otherwise have the rest of the
 * rectangle copied back as whatever the allocation happened to hold.
 *
 * The row-major leading dimension is always MAX(1, n): that is the smallest
 * value LAPACKE accepts, so the wrappers exercise the tightest buffer the
 * interface allows.
 */

/******************************************************************************/
/*                        ge: an m-by-n general matrix.                       */
/******************************************************************************/

/** Allocate a row-major shadow copy of an m-by-n general matrix. */
float *lapacke_test_sge_cm_to_rm(lapack_int m, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr)
{
    float *r;
    *ldr = MAX(1, n);
    r = (float *)LAPACKE_malloc(sizeof(float) * MAX(1, m) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_sge_trans)(LAPACK_COL_MAJOR, m, n, a, lda, r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into an m-by-n general matrix. */
void lapacke_test_sge_rm_to_cm(lapack_int m, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_sge_trans)(LAPACK_ROW_MAJOR, m, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of an m-by-n general matrix. */
double *lapacke_test_dge_cm_to_rm(lapack_int m, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, m) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dge_trans)(LAPACK_COL_MAJOR, m, n, a, lda, r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into an m-by-n general matrix. */
void lapacke_test_dge_rm_to_cm(lapack_int m, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_dge_trans)(LAPACK_ROW_MAJOR, m, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of an m-by-n general matrix. */
lapack_complex_float *lapacke_test_cge_cm_to_rm(lapack_int m, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               MAX(1, m) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_cge_trans)(LAPACK_COL_MAJOR, m, n, a, lda, r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into an m-by-n general matrix. */
void lapacke_test_cge_rm_to_cm(lapack_int m, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_cge_trans)(LAPACK_ROW_MAJOR, m, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of an m-by-n general matrix. */
lapack_complex_double *lapacke_test_zge_cm_to_rm(lapack_int m, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                MAX(1, m) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zge_trans)(LAPACK_COL_MAJOR, m, n, a, lda, r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into an m-by-n general matrix. */
void lapacke_test_zge_rm_to_cm(lapack_int m, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_zge_trans)(LAPACK_ROW_MAJOR, m, n, r, ldr, a, lda);
}

/******************************************************************************/
/*                     gb: an m-by-n general band matrix.                     */
/******************************************************************************/

/** Allocate a row-major shadow copy of an m-by-n general band matrix. */
float *lapacke_test_sgb_cm_to_rm(lapack_int m, lapack_int n, lapack_int kl,
                                 lapack_int ku, const float *a, lapack_int lda,
                                 lapack_int *ldr)
{
    float *r;
    *ldr = MAX(1, n);
    r = (float *)LAPACKE_malloc(sizeof(float) * (MAX(1, kl + ku + 1)) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_sgb_trans)(LAPACK_COL_MAJOR, m, n, kl, ku, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into an m-by-n general band matrix. */
void lapacke_test_sgb_rm_to_cm(lapack_int m, lapack_int n, lapack_int kl,
                               lapack_int ku, const float *r, lapack_int ldr,
                               float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_sgb_trans)(LAPACK_ROW_MAJOR, m, n, kl, ku, r, ldr, a,
                                  lda);
}

/** Allocate a row-major shadow copy of an m-by-n general band matrix. */
double *lapacke_test_dgb_cm_to_rm(lapack_int m, lapack_int n, lapack_int kl,
                                  lapack_int ku, const double *a,
                                  lapack_int lda, lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * (MAX(1, kl + ku + 1)) *
                                 (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dgb_trans)(LAPACK_COL_MAJOR, m, n, kl, ku, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into an m-by-n general band matrix. */
void lapacke_test_dgb_rm_to_cm(lapack_int m, lapack_int n, lapack_int kl,
                               lapack_int ku, const double *r, lapack_int ldr,
                               double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_dgb_trans)(LAPACK_ROW_MAJOR, m, n, kl, ku, r, ldr, a,
                                  lda);
}

/** Allocate a row-major shadow copy of an m-by-n general band matrix. */
lapack_complex_float *lapacke_test_cgb_cm_to_rm(lapack_int m, lapack_int n,
                                                lapack_int kl, lapack_int ku,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               (MAX(1, kl + ku + 1)) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_cgb_trans)(LAPACK_COL_MAJOR, m, n, kl, ku, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into an m-by-n general band matrix. */
void lapacke_test_cgb_rm_to_cm(lapack_int m, lapack_int n, lapack_int kl,
                               lapack_int ku, const lapack_complex_float *r,
                               lapack_int ldr, lapack_complex_float *a,
                               lapack_int lda)
{
    API_SUFFIX(LAPACKE_cgb_trans)(LAPACK_ROW_MAJOR, m, n, kl, ku, r, ldr, a,
                                  lda);
}

/** Allocate a row-major shadow copy of an m-by-n general band matrix. */
lapack_complex_double *lapacke_test_zgb_cm_to_rm(lapack_int m, lapack_int n,
                                                 lapack_int kl, lapack_int ku,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                (MAX(1, kl + ku + 1)) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zgb_trans)(LAPACK_COL_MAJOR, m, n, kl, ku, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into an m-by-n general band matrix. */
void lapacke_test_zgb_rm_to_cm(lapack_int m, lapack_int n, lapack_int kl,
                               lapack_int ku, const lapack_complex_double *r,
                               lapack_int ldr, lapack_complex_double *a,
                               lapack_int lda)
{
    API_SUFFIX(LAPACKE_zgb_trans)(LAPACK_ROW_MAJOR, m, n, kl, ku, r, ldr, a,
                                  lda);
}

/******************************************************************************/
/*       po: the uplo triangle of a symmetric positive definite matrix.       */
/******************************************************************************/

/** Allocate a row-major shadow copy of the uplo triangle of a symmetric
 * positive definite matrix. */
float *lapacke_test_spo_cm_to_rm(char uplo, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr)
{
    float *r;
    *ldr = MAX(1, n);
    r = (float *)LAPACKE_malloc(sizeof(float) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_spo_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a symmetric
 * positive definite matrix. */
void lapacke_test_spo_rm_to_cm(char uplo, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_spo_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of the uplo triangle of a symmetric
 * positive definite matrix. */
double *lapacke_test_dpo_cm_to_rm(char uplo, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dpo_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a symmetric
 * positive definite matrix. */
void lapacke_test_dpo_rm_to_cm(char uplo, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_dpo_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of the uplo triangle of a symmetric
 * positive definite matrix. */
lapack_complex_float *lapacke_test_cpo_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_cpo_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a symmetric
 * positive definite matrix. */
void lapacke_test_cpo_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_cpo_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of the uplo triangle of a symmetric
 * positive definite matrix. */
lapack_complex_double *lapacke_test_zpo_cm_to_rm(char uplo, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zpo_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a symmetric
 * positive definite matrix. */
void lapacke_test_zpo_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_zpo_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/******************************************************************************/
/*                sy: the uplo triangle of a symmetric matrix.                */
/******************************************************************************/

/** Allocate a row-major shadow copy of the uplo triangle of a symmetric matrix.
 */
float *lapacke_test_ssy_cm_to_rm(char uplo, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr)
{
    float *r;
    *ldr = MAX(1, n);
    r = (float *)LAPACKE_malloc(sizeof(float) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_ssy_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a symmetric
 * matrix. */
void lapacke_test_ssy_rm_to_cm(char uplo, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_ssy_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of the uplo triangle of a symmetric matrix.
 */
double *lapacke_test_dsy_cm_to_rm(char uplo, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dsy_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a symmetric
 * matrix. */
void lapacke_test_dsy_rm_to_cm(char uplo, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_dsy_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of the uplo triangle of a symmetric matrix.
 */
lapack_complex_float *lapacke_test_csy_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_csy_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a symmetric
 * matrix. */
void lapacke_test_csy_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_csy_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of the uplo triangle of a symmetric matrix.
 */
lapack_complex_double *lapacke_test_zsy_cm_to_rm(char uplo, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zsy_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a symmetric
 * matrix. */
void lapacke_test_zsy_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_zsy_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/******************************************************************************/
/*                he: the uplo triangle of a Hermitian matrix.                */
/******************************************************************************/

/** Allocate a row-major shadow copy of the uplo triangle of a Hermitian matrix.
 */
lapack_complex_float *lapacke_test_che_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_che_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a Hermitian
 * matrix. */
void lapacke_test_che_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_che_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/** Allocate a row-major shadow copy of the uplo triangle of a Hermitian matrix.
 */
lapack_complex_double *lapacke_test_zhe_cm_to_rm(char uplo, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zhe_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the uplo triangle of a Hermitian
 * matrix. */
void lapacke_test_zhe_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_zhe_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/******************************************************************************/
/*             tr: the referenced triangle of a triangular matrix.            */
/******************************************************************************/

/** Allocate a row-major shadow copy of the referenced triangle of a triangular
 * matrix. */
float *lapacke_test_str_cm_to_rm(char uplo, char diag, lapack_int n,
                                 const float *a, lapack_int lda,
                                 lapack_int *ldr)
{
    float *r;
    *ldr = MAX(1, n);
    r = (float *)LAPACKE_malloc(sizeof(float) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_str_trans)(LAPACK_COL_MAJOR, uplo, diag, n, a, lda,
                                      r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the referenced triangle of a
 * triangular matrix. */
void lapacke_test_str_rm_to_cm(char uplo, char diag, lapack_int n,
                               const float *r, lapack_int ldr, float *a,
                               lapack_int lda)
{
    API_SUFFIX(LAPACKE_str_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, r, ldr, a,
                                  lda);
}

/** Allocate a row-major shadow copy of the referenced triangle of a triangular
 * matrix. */
double *lapacke_test_dtr_cm_to_rm(char uplo, char diag, lapack_int n,
                                  const double *a, lapack_int lda,
                                  lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dtr_trans)(LAPACK_COL_MAJOR, uplo, diag, n, a, lda,
                                      r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the referenced triangle of a
 * triangular matrix. */
void lapacke_test_dtr_rm_to_cm(char uplo, char diag, lapack_int n,
                               const double *r, lapack_int ldr, double *a,
                               lapack_int lda)
{
    API_SUFFIX(LAPACKE_dtr_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, r, ldr, a,
                                  lda);
}

/** Allocate a row-major shadow copy of the referenced triangle of a triangular
 * matrix. */
lapack_complex_float *lapacke_test_ctr_cm_to_rm(char uplo, char diag,
                                                lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_ctr_trans)(LAPACK_COL_MAJOR, uplo, diag, n, a, lda,
                                      r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the referenced triangle of a
 * triangular matrix. */
void lapacke_test_ctr_rm_to_cm(char uplo, char diag, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_ctr_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, r, ldr, a,
                                  lda);
}

/** Allocate a row-major shadow copy of the referenced triangle of a triangular
 * matrix. */
lapack_complex_double *lapacke_test_ztr_cm_to_rm(char uplo, char diag,
                                                 lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_ztr_trans)(LAPACK_COL_MAJOR, uplo, diag, n, a, lda,
                                      r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into the referenced triangle of a
 * triangular matrix. */
void lapacke_test_ztr_rm_to_cm(char uplo, char diag, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_ztr_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, r, ldr, a,
                                  lda);
}

/******************************************************************************/
/*               pb: a symmetric positive definite band matrix.               */
/******************************************************************************/

/** Allocate a row-major shadow copy of a symmetric positive definite band
 * matrix. */
float *lapacke_test_spb_cm_to_rm(char uplo, lapack_int n, lapack_int kd,
                                 const float *a, lapack_int lda,
                                 lapack_int *ldr)
{
    float *r;
    *ldr = MAX(1, n);
    r = (float *)LAPACKE_malloc(sizeof(float) * MAX(1, kd + 1) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_spb_trans)(LAPACK_COL_MAJOR, uplo, n, kd, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a symmetric positive definite band
 * matrix. */
void lapacke_test_spb_rm_to_cm(char uplo, lapack_int n, lapack_int kd,
                               const float *r, lapack_int ldr, float *a,
                               lapack_int lda)
{
    API_SUFFIX(LAPACKE_spb_trans)(LAPACK_ROW_MAJOR, uplo, n, kd, r, ldr, a,
                                  lda);
}

/** Allocate a row-major shadow copy of a symmetric positive definite band
 * matrix. */
double *lapacke_test_dpb_cm_to_rm(char uplo, lapack_int n, lapack_int kd,
                                  const double *a, lapack_int lda,
                                  lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, kd + 1) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dpb_trans)(LAPACK_COL_MAJOR, uplo, n, kd, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a symmetric positive definite band
 * matrix. */
void lapacke_test_dpb_rm_to_cm(char uplo, lapack_int n, lapack_int kd,
                               const double *r, lapack_int ldr, double *a,
                               lapack_int lda)
{
    API_SUFFIX(LAPACKE_dpb_trans)(LAPACK_ROW_MAJOR, uplo, n, kd, r, ldr, a,
                                  lda);
}

/** Allocate a row-major shadow copy of a symmetric positive definite band
 * matrix. */
lapack_complex_float *lapacke_test_cpb_cm_to_rm(char uplo, lapack_int n,
                                                lapack_int kd,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               MAX(1, kd + 1) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_cpb_trans)(LAPACK_COL_MAJOR, uplo, n, kd, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a symmetric positive definite band
 * matrix. */
void lapacke_test_cpb_rm_to_cm(char uplo, lapack_int n, lapack_int kd,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_cpb_trans)(LAPACK_ROW_MAJOR, uplo, n, kd, r, ldr, a,
                                  lda);
}

/** Allocate a row-major shadow copy of a symmetric positive definite band
 * matrix. */
lapack_complex_double *lapacke_test_zpb_cm_to_rm(char uplo, lapack_int n,
                                                 lapack_int kd,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                MAX(1, kd + 1) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zpb_trans)(LAPACK_COL_MAJOR, uplo, n, kd, a, lda, r,
                                      *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a symmetric positive definite band
 * matrix. */
void lapacke_test_zpb_rm_to_cm(char uplo, lapack_int n, lapack_int kd,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_zpb_trans)(LAPACK_ROW_MAJOR, uplo, n, kd, r, ldr, a,
                                  lda);
}

/******************************************************************************/
/*                        tb: a triangular band matrix.                       */
/******************************************************************************/

/** Allocate a row-major shadow copy of a triangular band matrix. */
float *lapacke_test_stb_cm_to_rm(char uplo, char diag, lapack_int n,
                                 lapack_int kd, const float *a, lapack_int lda,
                                 lapack_int *ldr)
{
    float *r;
    *ldr = MAX(1, n);
    r = (float *)LAPACKE_malloc(sizeof(float) * MAX(1, kd + 1) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_stb_trans)(LAPACK_COL_MAJOR, uplo, diag, n, kd, a,
                                      lda, r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a triangular band matrix. */
void lapacke_test_stb_rm_to_cm(char uplo, char diag, lapack_int n,
                               lapack_int kd, const float *r, lapack_int ldr,
                               float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_stb_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, kd, r, ldr,
                                  a, lda);
}

/** Allocate a row-major shadow copy of a triangular band matrix. */
double *lapacke_test_dtb_cm_to_rm(char uplo, char diag, lapack_int n,
                                  lapack_int kd, const double *a,
                                  lapack_int lda, lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, kd + 1) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dtb_trans)(LAPACK_COL_MAJOR, uplo, diag, n, kd, a,
                                      lda, r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a triangular band matrix. */
void lapacke_test_dtb_rm_to_cm(char uplo, char diag, lapack_int n,
                               lapack_int kd, const double *r, lapack_int ldr,
                               double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_dtb_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, kd, r, ldr,
                                  a, lda);
}

/** Allocate a row-major shadow copy of a triangular band matrix. */
lapack_complex_float *lapacke_test_ctb_cm_to_rm(char uplo, char diag,
                                                lapack_int n, lapack_int kd,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               MAX(1, kd + 1) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_ctb_trans)(LAPACK_COL_MAJOR, uplo, diag, n, kd, a,
                                      lda, r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a triangular band matrix. */
void lapacke_test_ctb_rm_to_cm(char uplo, char diag, lapack_int n,
                               lapack_int kd, const lapack_complex_float *r,
                               lapack_int ldr, lapack_complex_float *a,
                               lapack_int lda)
{
    API_SUFFIX(LAPACKE_ctb_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, kd, r, ldr,
                                  a, lda);
}

/** Allocate a row-major shadow copy of a triangular band matrix. */
lapack_complex_double *lapacke_test_ztb_cm_to_rm(char uplo, char diag,
                                                 lapack_int n, lapack_int kd,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                MAX(1, kd + 1) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_ztb_trans)(LAPACK_COL_MAJOR, uplo, diag, n, kd, a,
                                      lda, r, *ldr);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a triangular band matrix. */
void lapacke_test_ztb_rm_to_cm(char uplo, char diag, lapack_int n,
                               lapack_int kd, const lapack_complex_double *r,
                               lapack_int ldr, lapack_complex_double *a,
                               lapack_int lda)
{
    API_SUFFIX(LAPACKE_ztb_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, kd, r, ldr,
                                  a, lda);
}

/******************************************************************************/
/*              pp: a packed symmetric positive definite matrix.              */
/******************************************************************************/

/** Allocate a row-major shadow copy of a packed symmetric positive definite
 * matrix. */
float *lapacke_test_spp_cm_to_rm(char uplo, lapack_int n, const float *ap)
{
    float *r;
    r = (float *)LAPACKE_malloc(sizeof(float) * (MAX(1, n) * MAX(2, n + 1)) /
                                2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_spp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed symmetric positive
 * definite matrix. */
void lapacke_test_spp_rm_to_cm(char uplo, lapack_int n, const float *r,
                               float *ap)
{
    API_SUFFIX(LAPACKE_spp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed symmetric positive definite
 * matrix. */
double *lapacke_test_dpp_cm_to_rm(char uplo, lapack_int n, const double *ap)
{
    double *r;
    r = (double *)LAPACKE_malloc(sizeof(double) * (MAX(1, n) * MAX(2, n + 1)) /
                                 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dpp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed symmetric positive
 * definite matrix. */
void lapacke_test_dpp_rm_to_cm(char uplo, lapack_int n, const double *r,
                               double *ap)
{
    API_SUFFIX(LAPACKE_dpp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed symmetric positive definite
 * matrix. */
lapack_complex_float *lapacke_test_cpp_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *ap)
{
    lapack_complex_float *r;
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               (MAX(1, n) * MAX(2, n + 1)) / 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_cpp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed symmetric positive
 * definite matrix. */
void lapacke_test_cpp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r,
                               lapack_complex_float *ap)
{
    API_SUFFIX(LAPACKE_cpp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed symmetric positive definite
 * matrix. */
lapack_complex_double *
lapacke_test_zpp_cm_to_rm(char uplo, lapack_int n,
                          const lapack_complex_double *ap)
{
    lapack_complex_double *r;
    r = (lapack_complex_double *)LAPACKE_malloc(
        sizeof(lapack_complex_double) * (MAX(1, n) * MAX(2, n + 1)) / 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zpp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed symmetric positive
 * definite matrix. */
void lapacke_test_zpp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r,
                               lapack_complex_double *ap)
{
    API_SUFFIX(LAPACKE_zpp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/******************************************************************************/
/*                       sp: a packed symmetric matrix.                       */
/******************************************************************************/

/** Allocate a row-major shadow copy of a packed symmetric matrix. */
float *lapacke_test_ssp_cm_to_rm(char uplo, lapack_int n, const float *ap)
{
    float *r;
    r = (float *)LAPACKE_malloc(sizeof(float) * (MAX(1, n) * MAX(2, n + 1)) /
                                2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_ssp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed symmetric matrix. */
void lapacke_test_ssp_rm_to_cm(char uplo, lapack_int n, const float *r,
                               float *ap)
{
    API_SUFFIX(LAPACKE_ssp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed symmetric matrix. */
double *lapacke_test_dsp_cm_to_rm(char uplo, lapack_int n, const double *ap)
{
    double *r;
    r = (double *)LAPACKE_malloc(sizeof(double) * (MAX(1, n) * MAX(2, n + 1)) /
                                 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dsp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed symmetric matrix. */
void lapacke_test_dsp_rm_to_cm(char uplo, lapack_int n, const double *r,
                               double *ap)
{
    API_SUFFIX(LAPACKE_dsp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed symmetric matrix. */
lapack_complex_float *lapacke_test_csp_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *ap)
{
    lapack_complex_float *r;
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               (MAX(1, n) * MAX(2, n + 1)) / 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_csp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed symmetric matrix. */
void lapacke_test_csp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r,
                               lapack_complex_float *ap)
{
    API_SUFFIX(LAPACKE_csp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed symmetric matrix. */
lapack_complex_double *
lapacke_test_zsp_cm_to_rm(char uplo, lapack_int n,
                          const lapack_complex_double *ap)
{
    lapack_complex_double *r;
    r = (lapack_complex_double *)LAPACKE_malloc(
        sizeof(lapack_complex_double) * (MAX(1, n) * MAX(2, n + 1)) / 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zsp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed symmetric matrix. */
void lapacke_test_zsp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r,
                               lapack_complex_double *ap)
{
    API_SUFFIX(LAPACKE_zsp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/******************************************************************************/
/*                       hp: a packed Hermitian matrix.                       */
/******************************************************************************/

/** Allocate a row-major shadow copy of a packed Hermitian matrix. */
lapack_complex_float *lapacke_test_chp_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *ap)
{
    lapack_complex_float *r;
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               (MAX(1, n) * MAX(2, n + 1)) / 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_chp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed Hermitian matrix. */
void lapacke_test_chp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r,
                               lapack_complex_float *ap)
{
    API_SUFFIX(LAPACKE_chp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed Hermitian matrix. */
lapack_complex_double *
lapacke_test_zhp_cm_to_rm(char uplo, lapack_int n,
                          const lapack_complex_double *ap)
{
    lapack_complex_double *r;
    r = (lapack_complex_double *)LAPACKE_malloc(
        sizeof(lapack_complex_double) * (MAX(1, n) * MAX(2, n + 1)) / 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zhp_trans)(LAPACK_COL_MAJOR, uplo, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed Hermitian matrix. */
void lapacke_test_zhp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r,
                               lapack_complex_double *ap)
{
    API_SUFFIX(LAPACKE_zhp_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ap);
}

/******************************************************************************/
/*                       tp: a packed triangular matrix.                      */
/******************************************************************************/

/** Allocate a row-major shadow copy of a packed triangular matrix. */
float *lapacke_test_stp_cm_to_rm(char uplo, char diag, lapack_int n,
                                 const float *ap)
{
    float *r;
    r = (float *)LAPACKE_malloc(sizeof(float) * (MAX(1, n) * MAX(2, n + 1)) /
                                2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_stp_trans)(LAPACK_COL_MAJOR, uplo, diag, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed triangular matrix. */
void lapacke_test_stp_rm_to_cm(char uplo, char diag, lapack_int n,
                               const float *r, float *ap)
{
    API_SUFFIX(LAPACKE_stp_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed triangular matrix. */
double *lapacke_test_dtp_cm_to_rm(char uplo, char diag, lapack_int n,
                                  const double *ap)
{
    double *r;
    r = (double *)LAPACKE_malloc(sizeof(double) * (MAX(1, n) * MAX(2, n + 1)) /
                                 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dtp_trans)(LAPACK_COL_MAJOR, uplo, diag, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed triangular matrix. */
void lapacke_test_dtp_rm_to_cm(char uplo, char diag, lapack_int n,
                               const double *r, double *ap)
{
    API_SUFFIX(LAPACKE_dtp_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed triangular matrix. */
lapack_complex_float *lapacke_test_ctp_cm_to_rm(char uplo, char diag,
                                                lapack_int n,
                                                const lapack_complex_float *ap)
{
    lapack_complex_float *r;
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               (MAX(1, n) * MAX(2, n + 1)) / 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_ctp_trans)(LAPACK_COL_MAJOR, uplo, diag, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed triangular matrix. */
void lapacke_test_ctp_rm_to_cm(char uplo, char diag, lapack_int n,
                               const lapack_complex_float *r,
                               lapack_complex_float *ap)
{
    API_SUFFIX(LAPACKE_ctp_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, r, ap);
}

/** Allocate a row-major shadow copy of a packed triangular matrix. */
lapack_complex_double *
lapacke_test_ztp_cm_to_rm(char uplo, char diag, lapack_int n,
                          const lapack_complex_double *ap)
{
    lapack_complex_double *r;
    r = (lapack_complex_double *)LAPACKE_malloc(
        sizeof(lapack_complex_double) * (MAX(1, n) * MAX(2, n + 1)) / 2);
    if (r != NULL) {
        API_SUFFIX(LAPACKE_ztp_trans)(LAPACK_COL_MAJOR, uplo, diag, n, ap, r);
    }
    return r;
}

/** Copy a row-major shadow buffer back into a packed triangular matrix. */
void lapacke_test_ztp_rm_to_cm(char uplo, char diag, lapack_int n,
                               const lapack_complex_double *r,
                               lapack_complex_double *ap)
{
    API_SUFFIX(LAPACKE_ztp_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, r, ap);
}

/**
 * \brief Number of rows xLASWP reaches, the way LAPACKE counts them.
 *
 * xLASWP takes no row count: the rows it touches follow from K1, K2 and the
 * pivot indices. LAPACKE sizes its row-major buffer from exactly this
 * expression, so the shadow copy has to use it too.
 *
 * \param[in] k1   First element of ipiv to apply.
 * \param[in] k2   Last element of ipiv to apply.
 * \param[in] ipiv Pivot indices.
 * \param[in] incx Stride (and direction) through ipiv.
 * \return The number of rows of the matrix.
 */
lapack_int lapacke_test_laswp_rows(lapack_int k1, lapack_int k2,
                                   const lapack_int *ipiv, lapack_int incx)
{
    lapack_int rows = MAX(1, k2);
    lapack_int i;
    for (i = k1; i <= k2; i++) {
        rows = MAX(rows, ipiv[k1 + (i - k1) * ABS(incx) - 1]);
    }
    return rows;
}

/**
 * \brief Report a failed shadow buffer allocation and set info to
 * LAPACK_TRANSPOSE_MEMORY_ERROR.
 *
 * Prints a diagnostic on stderr and sets info to
 * LAPACK_TRANSPOSE_MEMORY_ERROR, matching what LAPACKE itself returns when
 * its transposition buffers cannot be allocated.
 *
 * \param[in]  srname Routine name for the diagnostic.
 * \param[out] info   The wrapper's info result.
 */
void lapacke_test_report_alloc_failure(const char *srname, lapack_int *info)
{
    *info = LAPACK_TRANSPOSE_MEMORY_ERROR;
    fprintf(stderr,
            "LAPACKE test wrapper %s: shadow buffer allocation"
            " failed\n",
            srname);
}

#endif /* LAPACK_ROW_MAJOR */
