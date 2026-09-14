/******************************************************************************
 * Shared helpers for the LAPACKE test wrapper libraries.
 *
 * The wrappers (one source file per precision, lapacke_test_wrappers_<x>.c)
 * export Fortran-callable symbols <NAME>_TEST with the exact argument list
 * of the corresponding LAPACK routine. The LAPACK testing programs are
 * rewritten by the extended-API source generator so that calls to
 * allowlisted routines resolve to these wrappers, which forward to LAPACKE.
 * This way the unmodified numerical and error-exit checks of the Fortran
 * test suite exercise LAPACKE.
 *
 * The wrapper libraries are compiled once per (layout, layer) combination:
 *
 *   LAPACKE_TEST_LAYOUT=LAPACK_COL_MAJOR -- faithful passthrough
 *   LAPACKE_TEST_LAYOUT=LAPACK_ROW_MAJOR -- caller data is transposed into
 *     row-major shadow buffers so that LAPACKE's row-major conversion
 *     machinery is exercised; results are transposed back afterwards.
 *
 *   LAPACKE_TEST_LAYER=LAPACKE_TEST_LAYER_WORK -- call the middle-level
 *     (_work) interface; the caller's WORK/LWORK arguments pass through,
 *     which keeps minimal-workspace tests and LWORK error exits faithful to
 *     the Fortran routine.
 *   LAPACKE_TEST_LAYER=LAPACKE_TEST_LAYER_HIGH -- call the high-level
 *     interface, which performs NaN checks and allocates its own workspace.
 *     The caller's WORK/LWORK are ignored, so LWORK error exits cannot fire;
 *     suite runs for this layer must disable the error-exit tests.
 *
 * Both layers serve a workspace query (LWORK or TSIZE of -1) through the
 * middle-level interface in column-major, and take no row-major copies for
 * it: a query reads no matrix, only the dimensions determine the result, and
 * the caller's arrays are only guaranteed to hold the shapes the real call
 * will use; a query is routinely made with dimensions larger than that.
 *
 * The helpers declared here are implemented in lapacke_test_wrappers.c,
 * which is compiled into each wrapper library.
 ******************************************************************************/

#ifndef LAPACKE_TEST_WRAPPERS_H
#define LAPACKE_TEST_WRAPPERS_H

#include "lapacke.h"
#include "lapacke_utils.h"

#ifndef LAPACKE_TEST_LAYOUT
#error "LAPACKE_TEST_LAYOUT must be defined as LAPACK_{COL,ROW}_MAJOR"
#endif

/** LAPACKE_TEST_LAYER value selecting the middle-level (_work) API. */
#define LAPACKE_TEST_LAYER_WORK 1
/** LAPACKE_TEST_LAYER value selecting the high-level API. */
#define LAPACKE_TEST_LAYER_HIGH 2

#ifndef LAPACKE_TEST_LAYER
#error "LAPACKE_TEST_LAYER must be defined as LAPACKE_TEST_LAYER_{WORK,HIGH}"
#endif

/** Nonzero when the wrappers should call the high-level interface. */
#define LAPACKE_TEST_HIGH_LEVEL (LAPACKE_TEST_LAYER == LAPACKE_TEST_LAYER_HIGH)
/** Nonzero when the wrappers should use row-major layout. */
#define LAPACKE_TEST_ROW_MAJOR (LAPACKE_TEST_LAYOUT == LAPACK_ROW_MAJOR)

/** Map a LAPACKE info return value back to Fortran numbering and report
 *  argument errors through the testing XERBLA. */
lapack_int lapacke_test_info(const char *srname, lapack_int ret);

/** lapacke_test_info for the routines LAPACKE declares layout-independent,
 *  whose info carries no matrix_layout shift to undo. */
lapack_int lapacke_test_info_unshifted(const char *srname, lapack_int ret);

#if LAPACKE_TEST_LAYOUT == LAPACK_ROW_MAJOR

/* One cm_to_rm/rm_to_cm pair per storage scheme and precision, mirroring
 * LAPACKE's own <x><type>_trans usage: cm_to_rm allocates a row-major copy of
 * the caller's column-major matrix and returns NULL when the allocation
 * fails, rm_to_cm copies the result back. The row-major leading dimension
 * cm_to_rm reports is always MAX(1, n), the smallest value LAPACKE accepts.
 */

/******************************************************************************/
/*                        ge: an m-by-n general matrix.                       */
/******************************************************************************/
float *lapacke_test_sge_cm_to_rm(lapack_int m, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr);
void lapacke_test_sge_rm_to_cm(lapack_int m, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda);
double *lapacke_test_dge_cm_to_rm(lapack_int m, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr);
void lapacke_test_dge_rm_to_cm(lapack_int m, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda);
lapack_complex_float *lapacke_test_cge_cm_to_rm(lapack_int m, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);
void lapacke_test_cge_rm_to_cm(lapack_int m, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda);
lapack_complex_double *lapacke_test_zge_cm_to_rm(lapack_int m, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);
void lapacke_test_zge_rm_to_cm(lapack_int m, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda);

/******************************************************************************/
/*                     gb: an m-by-n general band matrix.                     */
/******************************************************************************/
float *lapacke_test_sgb_cm_to_rm(lapack_int m, lapack_int n, lapack_int kl,
                                 lapack_int ku, const float *a, lapack_int lda,
                                 lapack_int *ldr);
void lapacke_test_sgb_rm_to_cm(lapack_int m, lapack_int n, lapack_int kl,
                               lapack_int ku, const float *r, lapack_int ldr,
                               float *a, lapack_int lda);
double *lapacke_test_dgb_cm_to_rm(lapack_int m, lapack_int n, lapack_int kl,
                                  lapack_int ku, const double *a,
                                  lapack_int lda, lapack_int *ldr);
void lapacke_test_dgb_rm_to_cm(lapack_int m, lapack_int n, lapack_int kl,
                               lapack_int ku, const double *r, lapack_int ldr,
                               double *a, lapack_int lda);
lapack_complex_float *lapacke_test_cgb_cm_to_rm(lapack_int m, lapack_int n,
                                                lapack_int kl, lapack_int ku,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);
void lapacke_test_cgb_rm_to_cm(lapack_int m, lapack_int n, lapack_int kl,
                               lapack_int ku, const lapack_complex_float *r,
                               lapack_int ldr, lapack_complex_float *a,
                               lapack_int lda);
lapack_complex_double *lapacke_test_zgb_cm_to_rm(lapack_int m, lapack_int n,
                                                 lapack_int kl, lapack_int ku,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);
void lapacke_test_zgb_rm_to_cm(lapack_int m, lapack_int n, lapack_int kl,
                               lapack_int ku, const lapack_complex_double *r,
                               lapack_int ldr, lapack_complex_double *a,
                               lapack_int lda);

/******************************************************************************/
/*       po: the uplo triangle of a symmetric positive definite matrix.       */
/******************************************************************************/
float *lapacke_test_spo_cm_to_rm(char uplo, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr);
void lapacke_test_spo_rm_to_cm(char uplo, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda);
double *lapacke_test_dpo_cm_to_rm(char uplo, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr);
void lapacke_test_dpo_rm_to_cm(char uplo, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda);
lapack_complex_float *lapacke_test_cpo_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);
void lapacke_test_cpo_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda);
lapack_complex_double *lapacke_test_zpo_cm_to_rm(char uplo, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);
void lapacke_test_zpo_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda);

/******************************************************************************/
/*                sy: the uplo triangle of a symmetric matrix.                */
/******************************************************************************/
float *lapacke_test_ssy_cm_to_rm(char uplo, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr);
void lapacke_test_ssy_rm_to_cm(char uplo, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda);
double *lapacke_test_dsy_cm_to_rm(char uplo, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr);
void lapacke_test_dsy_rm_to_cm(char uplo, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda);
lapack_complex_float *lapacke_test_csy_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);
void lapacke_test_csy_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda);
lapack_complex_double *lapacke_test_zsy_cm_to_rm(char uplo, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);
void lapacke_test_zsy_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda);

/******************************************************************************/
/*                he: the uplo triangle of a Hermitian matrix.                */
/******************************************************************************/
lapack_complex_float *lapacke_test_che_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);
void lapacke_test_che_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda);
lapack_complex_double *lapacke_test_zhe_cm_to_rm(char uplo, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);
void lapacke_test_zhe_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda);

/******************************************************************************/
/*             tr: the referenced triangle of a triangular matrix.            */
/******************************************************************************/
float *lapacke_test_str_cm_to_rm(char uplo, char diag, lapack_int n,
                                 const float *a, lapack_int lda,
                                 lapack_int *ldr);
void lapacke_test_str_rm_to_cm(char uplo, char diag, lapack_int n,
                               const float *r, lapack_int ldr, float *a,
                               lapack_int lda);
double *lapacke_test_dtr_cm_to_rm(char uplo, char diag, lapack_int n,
                                  const double *a, lapack_int lda,
                                  lapack_int *ldr);
void lapacke_test_dtr_rm_to_cm(char uplo, char diag, lapack_int n,
                               const double *r, lapack_int ldr, double *a,
                               lapack_int lda);
lapack_complex_float *lapacke_test_ctr_cm_to_rm(char uplo, char diag,
                                                lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);
void lapacke_test_ctr_rm_to_cm(char uplo, char diag, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda);
lapack_complex_double *lapacke_test_ztr_cm_to_rm(char uplo, char diag,
                                                 lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);
void lapacke_test_ztr_rm_to_cm(char uplo, char diag, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda);

/******************************************************************************/
/*               pb: a symmetric positive definite band matrix.               */
/******************************************************************************/
float *lapacke_test_spb_cm_to_rm(char uplo, lapack_int n, lapack_int kd,
                                 const float *a, lapack_int lda,
                                 lapack_int *ldr);
void lapacke_test_spb_rm_to_cm(char uplo, lapack_int n, lapack_int kd,
                               const float *r, lapack_int ldr, float *a,
                               lapack_int lda);
double *lapacke_test_dpb_cm_to_rm(char uplo, lapack_int n, lapack_int kd,
                                  const double *a, lapack_int lda,
                                  lapack_int *ldr);
void lapacke_test_dpb_rm_to_cm(char uplo, lapack_int n, lapack_int kd,
                               const double *r, lapack_int ldr, double *a,
                               lapack_int lda);
lapack_complex_float *lapacke_test_cpb_cm_to_rm(char uplo, lapack_int n,
                                                lapack_int kd,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);
void lapacke_test_cpb_rm_to_cm(char uplo, lapack_int n, lapack_int kd,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda);
lapack_complex_double *lapacke_test_zpb_cm_to_rm(char uplo, lapack_int n,
                                                 lapack_int kd,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);
void lapacke_test_zpb_rm_to_cm(char uplo, lapack_int n, lapack_int kd,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda);

/******************************************************************************/
/*                        tb: a triangular band matrix.                       */
/******************************************************************************/
float *lapacke_test_stb_cm_to_rm(char uplo, char diag, lapack_int n,
                                 lapack_int kd, const float *a, lapack_int lda,
                                 lapack_int *ldr);
void lapacke_test_stb_rm_to_cm(char uplo, char diag, lapack_int n,
                               lapack_int kd, const float *r, lapack_int ldr,
                               float *a, lapack_int lda);
double *lapacke_test_dtb_cm_to_rm(char uplo, char diag, lapack_int n,
                                  lapack_int kd, const double *a,
                                  lapack_int lda, lapack_int *ldr);
void lapacke_test_dtb_rm_to_cm(char uplo, char diag, lapack_int n,
                               lapack_int kd, const double *r, lapack_int ldr,
                               double *a, lapack_int lda);
lapack_complex_float *lapacke_test_ctb_cm_to_rm(char uplo, char diag,
                                                lapack_int n, lapack_int kd,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);
void lapacke_test_ctb_rm_to_cm(char uplo, char diag, lapack_int n,
                               lapack_int kd, const lapack_complex_float *r,
                               lapack_int ldr, lapack_complex_float *a,
                               lapack_int lda);
lapack_complex_double *lapacke_test_ztb_cm_to_rm(char uplo, char diag,
                                                 lapack_int n, lapack_int kd,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);
void lapacke_test_ztb_rm_to_cm(char uplo, char diag, lapack_int n,
                               lapack_int kd, const lapack_complex_double *r,
                               lapack_int ldr, lapack_complex_double *a,
                               lapack_int lda);

/******************************************************************************/
/*              pp: a packed symmetric positive definite matrix.              */
/******************************************************************************/
float *lapacke_test_spp_cm_to_rm(char uplo, lapack_int n, const float *ap);
void lapacke_test_spp_rm_to_cm(char uplo, lapack_int n, const float *r,
                               float *ap);
double *lapacke_test_dpp_cm_to_rm(char uplo, lapack_int n, const double *ap);
void lapacke_test_dpp_rm_to_cm(char uplo, lapack_int n, const double *r,
                               double *ap);
lapack_complex_float *lapacke_test_cpp_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *ap);
void lapacke_test_cpp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r,
                               lapack_complex_float *ap);
lapack_complex_double *
lapacke_test_zpp_cm_to_rm(char uplo, lapack_int n,
                          const lapack_complex_double *ap);
void lapacke_test_zpp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r,
                               lapack_complex_double *ap);

/******************************************************************************/
/*                       sp: a packed symmetric matrix.                       */
/******************************************************************************/
float *lapacke_test_ssp_cm_to_rm(char uplo, lapack_int n, const float *ap);
void lapacke_test_ssp_rm_to_cm(char uplo, lapack_int n, const float *r,
                               float *ap);
double *lapacke_test_dsp_cm_to_rm(char uplo, lapack_int n, const double *ap);
void lapacke_test_dsp_rm_to_cm(char uplo, lapack_int n, const double *r,
                               double *ap);
lapack_complex_float *lapacke_test_csp_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *ap);
void lapacke_test_csp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r,
                               lapack_complex_float *ap);
lapack_complex_double *
lapacke_test_zsp_cm_to_rm(char uplo, lapack_int n,
                          const lapack_complex_double *ap);
void lapacke_test_zsp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r,
                               lapack_complex_double *ap);

/******************************************************************************/
/*                       hp: a packed Hermitian matrix.                       */
/******************************************************************************/
lapack_complex_float *lapacke_test_chp_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *ap);
void lapacke_test_chp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r,
                               lapack_complex_float *ap);
lapack_complex_double *
lapacke_test_zhp_cm_to_rm(char uplo, lapack_int n,
                          const lapack_complex_double *ap);
void lapacke_test_zhp_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r,
                               lapack_complex_double *ap);

/******************************************************************************/
/*                       tp: a packed triangular matrix.                      */
/******************************************************************************/
float *lapacke_test_stp_cm_to_rm(char uplo, char diag, lapack_int n,
                                 const float *ap);
void lapacke_test_stp_rm_to_cm(char uplo, char diag, lapack_int n,
                               const float *r, float *ap);
double *lapacke_test_dtp_cm_to_rm(char uplo, char diag, lapack_int n,
                                  const double *ap);
void lapacke_test_dtp_rm_to_cm(char uplo, char diag, lapack_int n,
                               const double *r, double *ap);
lapack_complex_float *lapacke_test_ctp_cm_to_rm(char uplo, char diag,
                                                lapack_int n,
                                                const lapack_complex_float *ap);
void lapacke_test_ctp_rm_to_cm(char uplo, char diag, lapack_int n,
                               const lapack_complex_float *r,
                               lapack_complex_float *ap);
lapack_complex_double *
lapacke_test_ztp_cm_to_rm(char uplo, char diag, lapack_int n,
                          const lapack_complex_double *ap);
void lapacke_test_ztp_rm_to_cm(char uplo, char diag, lapack_int n,
                               const lapack_complex_double *r,
                               lapack_complex_double *ap);

/** Number of rows xLASWP reaches, the way LAPACKE counts them. */
lapack_int lapacke_test_laswp_rows(lapack_int k1, lapack_int k2,
                                   const lapack_int *ipiv, lapack_int incx);

/** Report a failed shadow buffer allocation and set info to
 *  LAPACK_TRANSPOSE_MEMORY_ERROR. */
void lapacke_test_report_alloc_failure(const char *srname, lapack_int *info);

#endif /* LAPACK_ROW_MAJOR */
#endif /* LAPACKE_TEST_WRAPPERS_H */
