/*****************************************************************************
 * LAPACKE test wrappers (double precision, spike subset)
 *
 * Each wrapper exports a Fortran-callable symbol <NAME>_TEST with the exact
 * argument list of the corresponding LAPACK routine and forwards to
 * LAPACKE; see lapacke_test_wrappers.h for the wrapper scheme and the
 * (layout, layer) build combinations.
 *****************************************************************************/

#include <stdio.h>

#include "lapacke_test_wrappers.h"

/*****************************************************************************
 * DGETRF( M, N, A, LDA, IPIV, INFO )
 *****************************************************************************/
#define DGETRF_TEST LAPACK_GLOBAL(dgetrf_test, DGETRF_TEST)
void DGETRF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret;
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR, *m, *n, a, *lda, ipiv);
#else
    ret = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, *m, *n, a, *lda, ipiv);
#endif
#else
    lapack_int lda_r;
    double *a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGETRF", info);
        return;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, *m, *n, a_r, lda_r, ipiv);
#else
    ret = LAPACKE_dgetrf_work(LAPACK_ROW_MAJOR, *m, *n, a_r, lda_r, ipiv);
#endif
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGETRF", ret);
}

/*****************************************************************************
 * DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 *****************************************************************************/
#define DGETRS_TEST LAPACK_GLOBAL(dgetrs_test, DGETRS_TEST)
void DGETRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, const lapack_int *ipiv,
                 double *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret;
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgetrs(LAPACK_COL_MAJOR, *trans, *n, *nrhs, a, *lda, ipiv, b,
                         *ldb);
#else
    ret = LAPACKE_dgetrs_work(LAPACK_COL_MAJOR, *trans, *n, *nrhs, a, *lda,
                              ipiv, b, *ldb);
#endif
#else
    lapack_int lda_r, ldb_r;
    double *a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    double *b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGETRS", info);
        return;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, *trans, *n, *nrhs, a_r, lda_r, ipiv,
                         b_r, ldb_r);
#else
    ret = LAPACKE_dgetrs_work(LAPACK_ROW_MAJOR, *trans, *n, *nrhs, a_r, lda_r,
                              ipiv, b_r, ldb_r);
#endif
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGETRS", ret);
}

/*****************************************************************************
 * DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 *****************************************************************************/
#define DGETRI_TEST LAPACK_GLOBAL(dgetri_test, DGETRI_TEST)
void DGETRI_TEST(const lapack_int *n, double *a, const lapack_int *lda,
                 const lapack_int *ipiv, double *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret;
#if LAPACKE_TEST_HIGH_LEVEL_API
    if (*lwork == -1) {
        ret = LAPACKE_dgetri_work(LAPACK_COL_MAJOR, *n, a, *lda, ipiv, work,
                                  *lwork);
        *info = lapacke_test_info("DGETRI", ret);
        return;
    }
#endif
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgetri(LAPACK_COL_MAJOR, *n, a, *lda, ipiv);
#else
    ret =
        LAPACKE_dgetri_work(LAPACK_COL_MAJOR, *n, a, *lda, ipiv, work, *lwork);
#endif
#else
    lapack_int lda_r;
    double *a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGETRI", info);
        return;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR, *n, a_r, lda_r, ipiv);
#else
    ret = LAPACKE_dgetri_work(LAPACK_ROW_MAJOR, *n, a_r, lda_r, ipiv, work,
                              *lwork);
#endif
    lapacke_test_dge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGETRI", ret);
}

/*****************************************************************************
 * DPOTRF( UPLO, N, A, LDA, INFO )
 *****************************************************************************/
#define DPOTRF_TEST LAPACK_GLOBAL(dpotrf_test, DPOTRF_TEST)
void DPOTRF_TEST(const char *uplo, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret;
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dpotrf(LAPACK_COL_MAJOR, *uplo, *n, a, *lda);
#else
    ret = LAPACKE_dpotrf_work(LAPACK_COL_MAJOR, *uplo, *n, a, *lda);
#endif
#else
    lapack_int lda_r;
    double *a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DPOTRF", info);
        return;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, *uplo, *n, a_r, lda_r);
#else
    ret = LAPACKE_dpotrf_work(LAPACK_ROW_MAJOR, *uplo, *n, a_r, lda_r);
#endif
    lapacke_test_dpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DPOTRF", ret);
}

/*****************************************************************************
 * DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 *****************************************************************************/
#define DPOTRS_TEST LAPACK_GLOBAL(dpotrs_test, DPOTRS_TEST)
void DPOTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret;
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dpotrs(LAPACK_COL_MAJOR, *uplo, *n, *nrhs, a, *lda, b, *ldb);
#else
    ret = LAPACKE_dpotrs_work(LAPACK_COL_MAJOR, *uplo, *n, *nrhs, a, *lda, b,
                              *ldb);
#endif
#else
    lapack_int lda_r, ldb_r;
    double *a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    double *b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DPOTRS", info);
        return;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dpotrs(LAPACK_ROW_MAJOR, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                         ldb_r);
#else
    ret = LAPACKE_dpotrs_work(LAPACK_ROW_MAJOR, *uplo, *n, *nrhs, a_r, lda_r,
                              b_r, ldb_r);
#endif
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DPOTRS", ret);
}

/*****************************************************************************
 * DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
#define DGEQRF_TEST LAPACK_GLOBAL(dgeqrf_test, DGEQRF_TEST)
void DGEQRF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, double *tau, double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret;
#if LAPACKE_TEST_HIGH_LEVEL_API
    if (*lwork == -1) {
        ret = LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, *m, *n, a, *lda, tau, work,
                                  *lwork);
        *info = lapacke_test_info("DGEQRF", ret);
        return;
    }
#endif
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, *m, *n, a, *lda, tau);
#else
    ret = LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, *m, *n, a, *lda, tau, work,
                              *lwork);
#endif
#else
    lapack_int lda_r;
    double *a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGEQRF", info);
        return;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, *m, *n, a_r, lda_r, tau);
#else
    ret = LAPACKE_dgeqrf_work(LAPACK_ROW_MAJOR, *m, *n, a_r, lda_r, tau, work,
                              *lwork);
#endif
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGEQRF", ret);
}

/*****************************************************************************
 * DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
#define DORGQR_TEST LAPACK_GLOBAL(dorgqr_test, DORGQR_TEST)
void DORGQR_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 double *a, const lapack_int *lda, const double *tau,
                 double *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret;
#if LAPACKE_TEST_HIGH_LEVEL_API
    if (*lwork == -1) {
        ret = LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, *m, *n, *k, a, *lda, tau,
                                  work, *lwork);
        *info = lapacke_test_info("DORGQR", ret);
        return;
    }
#endif
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dorgqr(LAPACK_COL_MAJOR, *m, *n, *k, a, *lda, tau);
#else
    ret = LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, *m, *n, *k, a, *lda, tau, work,
                              *lwork);
#endif
#else
    lapack_int lda_r;
    double *a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DORGQR", info);
        return;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = LAPACKE_dorgqr_work(LAPACK_ROW_MAJOR, *m, *n, *k, a_r, lda_r, tau,
                              work, *lwork);
#endif
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DORGQR", ret);
}

/*****************************************************************************
 * DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 *****************************************************************************/
#define DORMQR_TEST LAPACK_GLOBAL(dormqr_test, DORMQR_TEST)
void DORMQR_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k, const double *a,
                 const lapack_int *lda, const double *tau, double *c,
                 const lapack_int *ldc, double *work, const lapack_int *lwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret;
#if LAPACKE_TEST_HIGH_LEVEL_API
    if (*lwork == -1) {
        ret = LAPACKE_dormqr_work(LAPACK_COL_MAJOR, *side, *trans, *m, *n, *k,
                                  a, *lda, tau, c, *ldc, work, *lwork);
        *info = lapacke_test_info("DORMQR", ret);
        return;
    }
#endif
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dormqr(LAPACK_COL_MAJOR, *side, *trans, *m, *n, *k, a, *lda,
                         tau, c, *ldc);
#else
    ret = LAPACKE_dormqr_work(LAPACK_COL_MAJOR, *side, *trans, *m, *n, *k, a,
                              *lda, tau, c, *ldc, work, *lwork);
#endif
#else
    lapack_int lda_r, ldc_r;
    lapack_int r = LAPACKE_lsame(*side, 'l') ? *m : *n;
    double *a_r = lapacke_test_dge_cm_to_rm(r, *k, a, *lda, &lda_r);
    double *c_r = lapacke_test_dge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("DORMQR", info);
        return;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dormqr(LAPACK_ROW_MAJOR, *side, *trans, *m, *n, *k, a_r,
                         lda_r, tau, c_r, ldc_r);
#else
    ret = LAPACKE_dormqr_work(LAPACK_ROW_MAJOR, *side, *trans, *m, *n, *k, a_r,
                              lda_r, tau, c_r, ldc_r, work, *lwork);
#endif
    lapacke_test_dge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free(a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("DORMQR", ret);
}

/*****************************************************************************
 * DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 *****************************************************************************/
#define DGELS_TEST LAPACK_GLOBAL(dgels_test, DGELS_TEST)
void DGELS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                const lapack_int *nrhs, double *a, const lapack_int *lda,
                double *b, const lapack_int *ldb, double *work,
                const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret;
#if LAPACKE_TEST_HIGH_LEVEL_API
    if (*lwork == -1) {
        ret = LAPACKE_dgels_work(LAPACK_COL_MAJOR, *trans, *m, *n, *nrhs, a,
                                 *lda, b, *ldb, work, *lwork);
        *info = lapacke_test_info("DGELS", ret);
        return;
    }
#endif
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgels(LAPACK_COL_MAJOR, *trans, *m, *n, *nrhs, a, *lda, b,
                        *ldb);
#else
    ret = LAPACKE_dgels_work(LAPACK_COL_MAJOR, *trans, *m, *n, *nrhs, a, *lda,
                             b, *ldb, work, *lwork);
#endif
#else
    lapack_int lda_r, ldb_r;
    lapack_int brows = MAX(*m, *n);
    double *a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    double *b_r = lapacke_test_dge_cm_to_rm(brows, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGELS", info);
        return;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    ret = LAPACKE_dgels(LAPACK_ROW_MAJOR, *trans, *m, *n, *nrhs, a_r, lda_r,
                        b_r, ldb_r);
#else
    ret = LAPACKE_dgels_work(LAPACK_ROW_MAJOR, *trans, *m, *n, *nrhs, a_r,
                             lda_r, b_r, ldb_r, work, *lwork);
#endif
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(brows, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGELS", ret);
}

/*****************************************************************************
 * DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
 *****************************************************************************/
#define DLANGE_TEST LAPACK_GLOBAL(dlange_test, DLANGE_TEST)
double DLANGE_TEST(const char *norm, const lapack_int *m, const lapack_int *n,
                   const double *a, const lapack_int *lda, double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len
#endif
)
{
#if LAPACKE_TEST_LAYOUT == LAPACK_COL_MAJOR
#if LAPACKE_TEST_HIGH_LEVEL_API
    return LAPACKE_dlange(LAPACK_COL_MAJOR, *norm, *m, *n, a, *lda);
#else
    return LAPACKE_dlange_work(LAPACK_COL_MAJOR, *norm, *m, *n, a, *lda, work);
#endif
#else
    lapack_int lda_r;
    double res;
    double *a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        fprintf(stderr, "LAPACKE test wrapper DLANGE: shadow buffer"
                        " allocation failed\n");
        return 0.0;
    }
#if LAPACKE_TEST_HIGH_LEVEL_API
    res = LAPACKE_dlange(LAPACK_ROW_MAJOR, *norm, *m, *n, a_r, lda_r);
#else
    res =
        LAPACKE_dlange_work(LAPACK_ROW_MAJOR, *norm, *m, *n, a_r, lda_r, work);
#endif
    LAPACKE_free(a_r);
    return res;
#endif
}
