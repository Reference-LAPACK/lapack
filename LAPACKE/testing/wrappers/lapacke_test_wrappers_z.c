/******************************************************************************
 * LAPACKE test wrappers (double precision complex)
 *
 * Each wrapper exports a Fortran-callable symbol <NAME>_TEST with the exact
 * argument list of the corresponding LAPACK routine and forwards to
 * LAPACKE; see lapacke_test_wrappers.h for the wrapper scheme and the
 * (layout, layer) build combinations.
 *
 * The wrappers are grouped by the LIN test path that drives them.
 ******************************************************************************/

#include <stdio.h>

#include "lapacke_test_wrappers.h"

static const int layout = LAPACKE_TEST_LAYOUT;

// ========================================================================== //
//                            General matrices (GE)                           //
// ========================================================================== //

/******************************************************************************
 * ZGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZGECON_TEST LAPACK_GLOBAL_SUFFIX(zgecon_test, ZGECON_TEST)
void ZGECON_TEST(const char *norm, const lapack_int *n,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const double *anorm, double *rcond,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGECON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgecon)(layout, *norm, *n, a_r, lda_r, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_zgecon_work)(layout, *norm, *n, a_r, lda_r, *anorm,
                                          rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("ZGECON", ret);
}

/******************************************************************************
 * ZGECXX( FACT, USESD, M, N, SESEL_ROWS, SEL_DESEL_COLS, KMAXFREE, ABSTOL,
 * RELTOL, A, LDA, K, MAXC2NRMK, RELMAXC2NRMK, FNRMK, IPIV, JPIV, TAU, C, LDC,
 * QRC, LDQRC, X, LDX, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
 ******************************************************************************/
#define ZGECXX_TEST LAPACK_GLOBAL_SUFFIX(zgecxx_test, ZGECXX_TEST)
void ZGECXX_TEST(const char *fact, const char *usesd, const lapack_int *m,
                 const lapack_int *n, const lapack_int *sesel_rows,
                 const lapack_int *sel_desel_cols, const lapack_int *kmaxfree,
                 const double *abstol, const double *reltol,
                 lapack_complex_double *a, const lapack_int *lda, lapack_int *k,
                 double *maxc2nrmk, double *relmaxc2nrmk, double *fnrmk,
                 lapack_int *ipiv, lapack_int *jpiv, lapack_complex_double *tau,
                 lapack_complex_double *c, const lapack_int *ldc,
                 lapack_complex_double *qrc, const lapack_int *ldqrc,
                 lapack_complex_double *x, const lapack_int *ldx,
                 lapack_complex_double *work, const lapack_int *lwork,
                 double *rwork, const lapack_int *lrwork, lapack_int *iwork,
                 const lapack_int *liwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN usesd_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgecxx_work)(
            LAPACK_COL_MAJOR, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
            (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a, *lda,
            k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c, *ldc, qrc,
            *ldqrc, x, *ldx, work, *lwork, rwork, *lrwork, iwork, *liwork);
        *info = lapacke_test_info("ZGECXX", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *c_r = c;
    lapack_int ldc_r = *ldc;
    lapack_complex_double *qrc_r = qrc;
    lapack_int ldqrc_r = *ldqrc;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    c_r = lapacke_test_zge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    qrc_r = lapacke_test_zge_cm_to_rm(*m, MIN(*m, *n), qrc, *ldqrc, &ldqrc_r);
    x_r = lapacke_test_zge_cm_to_rm(*m, *n, x, *ldx, &ldx_r);
    if (a_r == NULL || c_r == NULL || qrc_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(c_r);
        LAPACKE_free(qrc_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZGECXX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgecxx)(
        layout, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
        (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a_r, lda_r,
        k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c_r, ldc_r, qrc_r,
        ldqrc_r, x_r, ldx_r);
#else
    ret = API_SUFFIX(LAPACKE_zgecxx_work)(
        layout, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
        (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a_r, lda_r,
        k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c_r, ldc_r, qrc_r,
        ldqrc_r, x_r, ldx_r, work, *lwork, rwork, *lrwork, iwork, *liwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    lapacke_test_zge_rm_to_cm(*m, MIN(*m, *n), qrc_r, ldqrc_r, qrc, *ldqrc);
    lapacke_test_zge_rm_to_cm(*m, *n, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(c_r);
    LAPACKE_free(qrc_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZGECXX", ret);
}

/******************************************************************************
 * ZGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )
 ******************************************************************************/
#define ZGEEQU_TEST LAPACK_GLOBAL_SUFFIX(zgeequ_test, ZGEEQU_TEST)
void ZGEEQU_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_complex_double *a, const lapack_int *lda,
                 double *r, double *c, double *rowcnd, double *colcnd,
                 double *amax, lapack_int *info)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGEEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeequ)(layout, *m, *n, a_r, lda_r, r, c, rowcnd,
                                     colcnd, amax);
#else
    ret = API_SUFFIX(LAPACKE_zgeequ_work)(layout, *m, *n, a_r, lda_r, r, c,
                                          rowcnd, colcnd, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("ZGEEQU", ret);
}

/******************************************************************************
 * ZGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define ZGERFS_TEST LAPACK_GLOBAL_SUFFIX(zgerfs_test, ZGERFS_TEST)
void ZGERFS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *af, const lapack_int *ldaf,
                 const lapack_int *ipiv, const lapack_complex_double *b,
                 const lapack_int *ldb, lapack_complex_double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    af_r = lapacke_test_zge_cm_to_rm(*n, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZGERFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgerfs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     af_r, ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zgerfs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZGERFS", ret);
}

/******************************************************************************
 * ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZGESV_TEST LAPACK_GLOBAL_SUFFIX(zgesv_test, ZGESV_TEST)
void ZGESV_TEST(const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_double *a, const lapack_int *lda,
                lapack_int *ipiv, lapack_complex_double *b,
                const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGESV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgesv)(layout, *n, *nrhs, a_r, lda_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgesv_work)(layout, *n, *nrhs, a_r, lda_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGESV", ret);
}

/******************************************************************************
 * ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZGESVX_TEST LAPACK_GLOBAL_SUFFIX(zgesvx_test, ZGESVX_TEST)
void ZGESVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_double *a,
                 const lapack_int *lda, lapack_complex_double *af,
                 const lapack_int *ldaf, lapack_int *ipiv, char *equed,
                 double *r, double *c, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_complex_double *x,
                 const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    af_r = lapacke_test_zge_cm_to_rm(*n, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(af_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZGESVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgesvx)(
        layout, *fact, *trans, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, equed,
        r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, rwork);
#else
    ret = API_SUFFIX(LAPACKE_zgesvx_work)(
        layout, *fact, *trans, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, equed,
        r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(af_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZGESVX", ret);
}

/******************************************************************************
 * ZGETF2( M, N, A, LDA, IPIV, INFO )
 ******************************************************************************/
#define ZGETF2_TEST LAPACK_GLOBAL_SUFFIX(zgetf2_test, ZGETF2_TEST)
void ZGETF2_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGETF2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgetf2)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zgetf2_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGETF2", ret);
}

/******************************************************************************
 * ZGETRF( M, N, A, LDA, IPIV, INFO )
 ******************************************************************************/
#define ZGETRF_TEST LAPACK_GLOBAL_SUFFIX(zgetrf_test, ZGETRF_TEST)
void ZGETRF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgetrf)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zgetrf_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGETRF", ret);
}

/******************************************************************************
 * ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGETRI_TEST LAPACK_GLOBAL_SUFFIX(zgetri_test, ZGETRI_TEST)
void ZGETRI_TEST(const lapack_int *n, lapack_complex_double *a,
                 const lapack_int *lda, const lapack_int *ipiv,
                 lapack_complex_double *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgetri_work)(LAPACK_COL_MAJOR, *n, a, *lda,
                                              ipiv, work, *lwork);
        *info = lapacke_test_info("ZGETRI", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgetri)(layout, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zgetri_work)(layout, *n, a_r, lda_r, ipiv, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGETRI", ret);
}

/******************************************************************************
 * ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZGETRS_TEST LAPACK_GLOBAL_SUFFIX(zgetrs_test, ZGETRS_TEST)
void ZGETRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgetrs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgetrs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGETRS", ret);
}

// ========================================================================== //
//                         General band matrices (GB)                         //
// ========================================================================== //

/******************************************************************************
 * ZGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZGBCON_TEST LAPACK_GLOBAL_SUFFIX(zgbcon_test, ZGBCON_TEST)
void ZGBCON_TEST(const char *norm, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_complex_double *ab,
                 const lapack_int *ldab, const lapack_int *ipiv,
                 const double *anorm, double *rcond,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("ZGBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgbcon)(layout, *norm, *n, *kl, *ku, ab_r, ldab_r,
                                     ipiv, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_zgbcon_work)(layout, *norm, *n, *kl, *ku, ab_r,
                                          ldab_r, ipiv, *anorm, rcond, work,
                                          rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("ZGBCON", ret);
}

/******************************************************************************
 * ZGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFO )
 ******************************************************************************/
#define ZGBEQU_TEST LAPACK_GLOBAL_SUFFIX(zgbequ_test, ZGBEQU_TEST)
void ZGBEQU_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_complex_double *ab,
                 const lapack_int *ldab, double *r, double *c, double *rowcnd,
                 double *colcnd, double *amax, lapack_int *info)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zgb_cm_to_rm(*m, *n, *kl, *ku, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("ZGBEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgbequ)(layout, *m, *n, *kl, *ku, ab_r, ldab_r, r,
                                     c, rowcnd, colcnd, amax);
#else
    ret = API_SUFFIX(LAPACKE_zgbequ_work)(layout, *m, *n, *kl, *ku, ab_r,
                                          ldab_r, r, c, rowcnd, colcnd, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("ZGBEQU", ret);
}

/******************************************************************************
 * ZGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX,
 * FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZGBRFS_TEST LAPACK_GLOBAL_SUFFIX(zgbrfs_test, ZGBRFS_TEST)
void ZGBRFS_TEST(const char *trans, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_int *nrhs,
                 const lapack_complex_double *ab, const lapack_int *ldab,
                 const lapack_complex_double *afb, const lapack_int *ldafb,
                 const lapack_int *ipiv, const lapack_complex_double *b,
                 const lapack_int *ldb, lapack_complex_double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const lapack_complex_double *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zgb_cm_to_rm(*n, *n, *kl, *ku, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_zgb_cm_to_rm(*n, *n, *kl, *kl + *ku, afb, *ldafb,
                                      &ldafb_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)afb_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZGBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgbrfs)(layout, *trans, *n, *kl, *ku, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, ipiv, b_r, ldb_r,
                                     x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zgbrfs_work)(
        layout, *trans, *n, *kl, *ku, *nrhs, ab_r, ldab_r, afb_r, ldafb_r, ipiv,
        b_r, ldb_r, x_r, ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)afb_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZGBRFS", ret);
}

/******************************************************************************
 * ZGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZGBSV_TEST LAPACK_GLOBAL_SUFFIX(zgbsv_test, ZGBSV_TEST)
void ZGBSV_TEST(const lapack_int *n, const lapack_int *kl, const lapack_int *ku,
                const lapack_int *nrhs, lapack_complex_double *ab,
                const lapack_int *ldab, lapack_int *ipiv,
                lapack_complex_double *b, const lapack_int *ldb,
                lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGBSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgbsv)(layout, *n, *kl, *ku, *nrhs, ab_r, ldab_r,
                                    ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgbsv_work)(layout, *n, *kl, *ku, *nrhs, ab_r,
                                         ldab_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zgb_rm_to_cm(*n, *n, *kl, *kl + *ku, ab_r, ldab_r, ab, *ldab);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGBSV", ret);
}

/******************************************************************************
 * ZGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, EQUED, R,
 * C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZGBSVX_TEST LAPACK_GLOBAL_SUFFIX(zgbsvx_test, ZGBSVX_TEST)
void ZGBSVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *kl, const lapack_int *ku,
                 const lapack_int *nrhs, lapack_complex_double *ab,
                 const lapack_int *ldab, lapack_complex_double *afb,
                 const lapack_int *ldafb, lapack_int *ipiv, char *equed,
                 double *r, double *c, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_complex_double *x,
                 const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_double *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zgb_cm_to_rm(*n, *n, *kl, *ku, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_zgb_cm_to_rm(*n, *n, *kl, *kl + *ku, afb, *ldafb,
                                      &ldafb_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(afb_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZGBSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgbsvx)(layout, *fact, *trans, *n, *kl, *ku, *nrhs,
                                     ab_r, ldab_r, afb_r, ldafb_r, ipiv, equed,
                                     r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr,
                                     berr, rwork);
#else
    ret = API_SUFFIX(LAPACKE_zgbsvx_work)(
        layout, *fact, *trans, *n, *kl, *ku, *nrhs, ab_r, ldab_r, afb_r,
        ldafb_r, ipiv, equed, r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr,
        work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zgb_rm_to_cm(*n, *n, *kl, *ku, ab_r, ldab_r, ab, *ldab);
    lapacke_test_zgb_rm_to_cm(*n, *n, *kl, *kl + *ku, afb_r, ldafb_r, afb,
                              *ldafb);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(ab_r);
    LAPACKE_free(afb_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZGBSVX", ret);
}

/******************************************************************************
 * ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
 ******************************************************************************/
#define ZGBTRF_TEST LAPACK_GLOBAL_SUFFIX(zgbtrf_test, ZGBTRF_TEST)
void ZGBTRF_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, lapack_complex_double *ab,
                 const lapack_int *ldab, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zgb_cm_to_rm(*m, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("ZGBTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgbtrf)(layout, *m, *n, *kl, *ku, ab_r, ldab_r,
                                     ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zgbtrf_work)(layout, *m, *n, *kl, *ku, ab_r,
                                          ldab_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zgb_rm_to_cm(*m, *n, *kl, *kl + *ku, ab_r, ldab_r, ab, *ldab);
    LAPACKE_free(ab_r);
#endif
    *info = lapacke_test_info("ZGBTRF", ret);
}

/******************************************************************************
 * ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZGBTRS_TEST LAPACK_GLOBAL_SUFFIX(zgbtrs_test, ZGBTRS_TEST)
void ZGBTRS_TEST(const char *trans, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_int *nrhs,
                 const lapack_complex_double *ab, const lapack_int *ldab,
                 const lapack_int *ipiv, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgbtrs)(layout, *trans, *n, *kl, *ku, *nrhs, ab_r,
                                     ldab_r, ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgbtrs_work)(layout, *trans, *n, *kl, *ku, *nrhs,
                                          ab_r, ldab_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGBTRS", ret);
}

// ========================================================================== //
//                      General tridiagonal matrices (GT)                     //
// ========================================================================== //

/******************************************************************************
 * ZGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define ZGTCON_TEST LAPACK_GLOBAL_SUFFIX(zgtcon_test, ZGTCON_TEST)
void ZGTCON_TEST(const char *norm, const lapack_int *n,
                 const lapack_complex_double *dl,
                 const lapack_complex_double *d,
                 const lapack_complex_double *du,
                 const lapack_complex_double *du2, const lapack_int *ipiv,
                 const double *anorm, double *rcond,
                 lapack_complex_double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgtcon)(*norm, *n, dl, d, du, du2, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_zgtcon_work)(*norm, *n, dl, d, du, du2, ipiv,
                                          *anorm, rcond, work);
#endif

    *info = lapacke_test_info_unshifted("ZGTCON", ret);
}

/******************************************************************************
 * ZGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX,
 * FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZGTRFS_TEST LAPACK_GLOBAL_SUFFIX(zgtrfs_test, ZGTRFS_TEST)
void
ZGTRFS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
            const lapack_complex_double *dl, const lapack_complex_double *d,
            const lapack_complex_double *du, const lapack_complex_double *dlf,
            const lapack_complex_double *df, const lapack_complex_double *duf,
            const lapack_complex_double *du2, const lapack_int *ipiv,
            const lapack_complex_double *b, const lapack_int *ldb,
            lapack_complex_double *x, const lapack_int *ldx, double *ferr,
            double *berr, lapack_complex_double *work, double *rwork,
            lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
            ,
            FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZGTRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgtrfs)(layout, *trans, *n, *nrhs, dl, d, du, dlf,
                                     df, duf, du2, ipiv, b_r, ldb_r, x_r, ldx_r,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zgtrfs_work)(layout, *trans, *n, *nrhs, dl, d, du,
                                          dlf, df, duf, du2, ipiv, b_r, ldb_r,
                                          x_r, ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZGTRFS", ret);
}

/******************************************************************************
 * ZGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
 ******************************************************************************/
#define ZGTSV_TEST LAPACK_GLOBAL_SUFFIX(zgtsv_test, ZGTSV_TEST)
void ZGTSV_TEST(const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_double *dl, lapack_complex_double *d,
                lapack_complex_double *du, lapack_complex_double *b,
                const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("ZGTSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgtsv)(layout, *n, *nrhs, dl, d, du, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgtsv_work)(layout, *n, *nrhs, dl, d, du, b_r,
                                         ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGTSV", ret);
}

/******************************************************************************
 * ZGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZGTSVX_TEST LAPACK_GLOBAL_SUFFIX(zgtsvx_test, ZGTSVX_TEST)
void ZGTSVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_double *dl,
                 const lapack_complex_double *d,
                 const lapack_complex_double *du, lapack_complex_double *dlf,
                 lapack_complex_double *df, lapack_complex_double *duf,
                 lapack_complex_double *du2, lapack_int *ipiv,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *rcond,
                 double *ferr, double *berr, lapack_complex_double *work,
                 double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZGTSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgtsvx)(layout, *fact, *trans, *n, *nrhs, dl, d,
                                     du, dlf, df, duf, du2, ipiv, b_r, ldb_r,
                                     x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zgtsvx_work)(
        layout, *fact, *trans, *n, *nrhs, dl, d, du, dlf, df, duf, du2, ipiv,
        b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZGTSVX", ret);
}

/******************************************************************************
 * ZGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
 ******************************************************************************/
#define ZGTTRF_TEST LAPACK_GLOBAL_SUFFIX(zgttrf_test, ZGTTRF_TEST)
void ZGTTRF_TEST(const lapack_int *n, lapack_complex_double *dl,
                 lapack_complex_double *d, lapack_complex_double *du,
                 lapack_complex_double *du2, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgttrf)(*n, dl, d, du, du2, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zgttrf_work)(*n, dl, d, du, du2, ipiv);
#endif

    *info = lapacke_test_info_unshifted("ZGTTRF", ret);
}

/******************************************************************************
 * ZGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZGTTRS_TEST LAPACK_GLOBAL_SUFFIX(zgttrs_test, ZGTTRS_TEST)
void ZGTTRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *dl,
                 const lapack_complex_double *d,
                 const lapack_complex_double *du,
                 const lapack_complex_double *du2, const lapack_int *ipiv,
                 lapack_complex_double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("ZGTTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgttrs)(layout, *trans, *n, *nrhs, dl, d, du, du2,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgttrs_work)(layout, *trans, *n, *nrhs, dl, d, du,
                                          du2, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGTTRS", ret);
}

// ========================================================================== //
//             Symmetric/Hermitian positive definite matrices (PO)            //
// ========================================================================== //

/******************************************************************************
 * ZPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZPOCON_TEST LAPACK_GLOBAL_SUFFIX(zpocon_test, ZPOCON_TEST)
void ZPOCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const double *anorm, double *rcond,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZPOCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpocon)(layout, *uplo, *n, a_r, lda_r, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_zpocon_work)(layout, *uplo, *n, a_r, lda_r, *anorm,
                                          rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("ZPOCON", ret);
}

/******************************************************************************
 * ZPOEQU( N, A, LDA, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define ZPOEQU_TEST LAPACK_GLOBAL_SUFFIX(zpoequ_test, ZPOEQU_TEST)
void ZPOEQU_TEST(const lapack_int *n, const lapack_complex_double *a,
                 const lapack_int *lda, double *s, double *scond, double *amax,
                 lapack_int *info)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZPOEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpoequ)(layout, *n, a_r, lda_r, s, scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_zpoequ_work)(layout, *n, a_r, lda_r, s, scond,
                                          amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("ZPOEQU", ret);
}

/******************************************************************************
 * ZPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define ZPORFS_TEST LAPACK_GLOBAL_SUFFIX(zporfs_test, ZPORFS_TEST)
void ZPORFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *af, const lapack_int *ldaf,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *ferr,
                 double *berr, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZPORFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zporfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_zporfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, b_r, ldb_r, x_r, ldx_r,
                                          ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZPORFS", ret);
}

/******************************************************************************
 * ZPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define ZPOSV_TEST LAPACK_GLOBAL_SUFFIX(zposv_test, ZPOSV_TEST)
void ZPOSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_double *a, const lapack_int *lda,
                lapack_complex_double *b, const lapack_int *ldb,
                lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZPOSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zposv)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zposv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZPOSV", ret);
}

/******************************************************************************
 * ZPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S, B, LDB, X, LDX,
 * RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZPOSVX_TEST LAPACK_GLOBAL_SUFFIX(zposvx_test, ZPOSVX_TEST)
void ZPOSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_double *a,
                 const lapack_int *lda, lapack_complex_double *af,
                 const lapack_int *ldaf, char *equed, double *s,
                 lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *rcond,
                 double *ferr, double *berr, lapack_complex_double *work,
                 double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(af_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZPOSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zposvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, equed, s, b_r, ldb_r,
                                     x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zposvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, equed, s,
        b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_zpo_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(af_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZPOSVX", ret);
}

/******************************************************************************
 * ZPOTRF( UPLO, N, A, LDA, INFO )
 ******************************************************************************/
#define ZPOTRF_TEST LAPACK_GLOBAL_SUFFIX(zpotrf_test, ZPOTRF_TEST)
void ZPOTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZPOTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpotrf)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_zpotrf_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZPOTRF", ret);
}

/******************************************************************************
 * ZPOTRI( UPLO, N, A, LDA, INFO )
 ******************************************************************************/
#define ZPOTRI_TEST LAPACK_GLOBAL_SUFFIX(zpotri_test, ZPOTRI_TEST)
void ZPOTRI_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZPOTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpotri)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_zpotri_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZPOTRI", ret);
}

/******************************************************************************
 * ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define ZPOTRS_TEST LAPACK_GLOBAL_SUFFIX(zpotrs_test, ZPOTRS_TEST)
void ZPOTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZPOTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpotrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zpotrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZPOTRS", ret);
}

/******************************************************************************
 * ZPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
 ******************************************************************************/
#define ZPSTRF_TEST LAPACK_GLOBAL_SUFFIX(zpstrf_test, ZPSTRF_TEST)
void ZPSTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *piv, lapack_int *rank, const double *tol,
                 double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZPSTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpstrf)(layout, *uplo, *n, a_r, lda_r, piv, rank,
                                     *tol);
#else
    ret = API_SUFFIX(LAPACKE_zpstrf_work)(layout, *uplo, *n, a_r, lda_r, piv,
                                          rank, *tol, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZPSTRF", ret);
}

// ========================================================================== //
//                   Packed positive definite matrices (PP)                   //
// ========================================================================== //

/******************************************************************************
 * ZPPCON( UPLO, N, AP, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZPPCON_TEST LAPACK_GLOBAL_SUFFIX(zppcon_test, ZPPCON_TEST)
void ZPPCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_double *ap, const double *anorm,
                 double *rcond, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZPPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zppcon)(layout, *uplo, *n, ap_r, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_zppcon_work)(layout, *uplo, *n, ap_r, *anorm,
                                          rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("ZPPCON", ret);
}

/******************************************************************************
 * ZPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define ZPPEQU_TEST LAPACK_GLOBAL_SUFFIX(zppequ_test, ZPPEQU_TEST)
void ZPPEQU_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_double *ap, double *s, double *scond,
                 double *amax, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZPPEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zppequ)(layout, *uplo, *n, ap_r, s, scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_zppequ_work)(layout, *uplo, *n, ap_r, s, scond,
                                          amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("ZPPEQU", ret);
}

/******************************************************************************
 * ZPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO
 * )
 ******************************************************************************/
#define ZPPRFS_TEST LAPACK_GLOBAL_SUFFIX(zpprfs_test, ZPPRFS_TEST)
void ZPPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *ap,
                 const lapack_complex_double *afp,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *ferr,
                 double *berr, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
    const lapack_complex_double *ap_r = ap;
    const lapack_complex_double *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("ZPPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r, b_r,
                                     ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zpprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          b_r, ldb_r, x_r, ldx_r, ferr, berr,
                                          work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("ZPPRFS", ret);
}

/******************************************************************************
 * ZPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
 ******************************************************************************/
#define ZPPSV_TEST LAPACK_GLOBAL_SUFFIX(zppsv_test, ZPPSV_TEST)
void ZPPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_double *ap, lapack_complex_double *b,
                const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("ZPPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zppsv)(layout, *uplo, *n, *nrhs, ap_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zppsv_work)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                         ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_zpp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("ZPPSV", ret);
}

/******************************************************************************
 * ZPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZPPSVX_TEST LAPACK_GLOBAL_SUFFIX(zppsvx_test, ZPPSVX_TEST)
void ZPPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_double *ap,
                 lapack_complex_double *afp, char *equed, double *s,
                 lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *rcond,
                 double *ferr, double *berr, lapack_complex_double *work,
                 double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *ap_r = ap;
    lapack_complex_double *afp_r = afp;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZPPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zppsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, equed, s, b_r, ldb_r, x_r, ldx_r,
                                     rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zppsvx_work)(
        layout, *fact, *uplo, *n, *nrhs, ap_r, afp_r, equed, s, b_r, ldb_r, x_r,
        ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_zpp_rm_to_cm(*uplo, *n, ap_r, ap);
    lapacke_test_zpp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZPPSVX", ret);
}

/******************************************************************************
 * ZPPTRF( UPLO, N, AP, INFO )
 ******************************************************************************/
#define ZPPTRF_TEST LAPACK_GLOBAL_SUFFIX(zpptrf_test, ZPPTRF_TEST)
void ZPPTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *ap, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZPPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpptrf)(layout, *uplo, *n, ap_r);
#else
    ret = API_SUFFIX(LAPACKE_zpptrf_work)(layout, *uplo, *n, ap_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("ZPPTRF", ret);
}

/******************************************************************************
 * ZPPTRI( UPLO, N, AP, INFO )
 ******************************************************************************/
#define ZPPTRI_TEST LAPACK_GLOBAL_SUFFIX(zpptri_test, ZPPTRI_TEST)
void ZPPTRI_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *ap, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZPPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpptri)(layout, *uplo, *n, ap_r);
#else
    ret = API_SUFFIX(LAPACKE_zpptri_work)(layout, *uplo, *n, ap_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("ZPPTRI", ret);
}

/******************************************************************************
 * ZPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
 ******************************************************************************/
#define ZPPTRS_TEST LAPACK_GLOBAL_SUFFIX(zpptrs_test, ZPPTRS_TEST)
void ZPPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *ap, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_zpp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("ZPPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpptrs)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zpptrs_work)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                          ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("ZPPTRS", ret);
}

// ========================================================================== //
//                    Positive definite band matrices (PB)                    //
// ========================================================================== //

/******************************************************************************
 * ZPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZPBCON_TEST LAPACK_GLOBAL_SUFFIX(zpbcon_test, ZPBCON_TEST)
void ZPBCON_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_complex_double *ab, const lapack_int *ldab,
                 const double *anorm, double *rcond,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("ZPBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpbcon)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_zpbcon_work)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                          *anorm, rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("ZPBCON", ret);
}

/******************************************************************************
 * ZPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define ZPBEQU_TEST LAPACK_GLOBAL_SUFFIX(zpbequ_test, ZPBEQU_TEST)
void ZPBEQU_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_complex_double *ab, const lapack_int *ldab,
                 double *s, double *scond, double *amax, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("ZPBEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpbequ)(layout, *uplo, *n, *kd, ab_r, ldab_r, s,
                                     scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_zpbequ_work)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                          s, scond, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("ZPBEQU", ret);
}

/******************************************************************************
 * ZPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define ZPBRFS_TEST LAPACK_GLOBAL_SUFFIX(zpbrfs_test, ZPBRFS_TEST)
void ZPBRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const lapack_complex_double *ab,
                 const lapack_int *ldab, const lapack_complex_double *afb,
                 const lapack_int *ldafb, const lapack_complex_double *b,
                 const lapack_int *ldb, lapack_complex_double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const lapack_complex_double *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_zpb_cm_to_rm(*uplo, *n, *kd, afb, *ldafb, &ldafb_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)afb_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZPBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpbrfs)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, b_r, ldb_r, x_r,
                                     ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zpbrfs_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                          ldab_r, afb_r, ldafb_r, b_r, ldb_r,
                                          x_r, ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)afb_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZPBRFS", ret);
}

/******************************************************************************
 * ZPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define ZPBSV_TEST LAPACK_GLOBAL_SUFFIX(zpbsv_test, ZPBSV_TEST)
void ZPBSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                const lapack_int *nrhs, lapack_complex_double *ab,
                const lapack_int *ldab, lapack_complex_double *b,
                const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZPBSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpbsv)(layout, *uplo, *n, *kd, *nrhs, ab_r, ldab_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zpbsv_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                         ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZPBSV", ret);
}

/******************************************************************************
 * ZPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZPBSVX_TEST LAPACK_GLOBAL_SUFFIX(zpbsvx_test, ZPBSVX_TEST)
void ZPBSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *kd, const lapack_int *nrhs,
                 lapack_complex_double *ab, const lapack_int *ldab,
                 lapack_complex_double *afb, const lapack_int *ldafb,
                 char *equed, double *s, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_complex_double *x,
                 const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_double *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_zpb_cm_to_rm(*uplo, *n, *kd, afb, *ldafb, &ldafb_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(afb_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZPBSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpbsvx)(layout, *fact, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, equed, s, b_r,
                                     ldb_r, x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zpbsvx_work)(
        layout, *fact, *uplo, *n, *kd, *nrhs, ab_r, ldab_r, afb_r, ldafb_r,
        equed, s, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    lapacke_test_zpb_rm_to_cm(*uplo, *n, *kd, afb_r, ldafb_r, afb, *ldafb);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(ab_r);
    LAPACKE_free(afb_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZPBSVX", ret);
}

/******************************************************************************
 * ZPBTRF( UPLO, N, KD, AB, LDAB, INFO )
 ******************************************************************************/
#define ZPBTRF_TEST LAPACK_GLOBAL_SUFFIX(zpbtrf_test, ZPBTRF_TEST)
void ZPBTRF_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 lapack_complex_double *ab, const lapack_int *ldab,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("ZPBTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpbtrf)(layout, *uplo, *n, *kd, ab_r, ldab_r);
#else
    ret = API_SUFFIX(LAPACKE_zpbtrf_work)(layout, *uplo, *n, *kd, ab_r, ldab_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    LAPACKE_free(ab_r);
#endif
    *info = lapacke_test_info("ZPBTRF", ret);
}

/******************************************************************************
 * ZPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define ZPBTRS_TEST LAPACK_GLOBAL_SUFFIX(zpbtrs_test, ZPBTRS_TEST)
void ZPBTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const lapack_complex_double *ab,
                 const lapack_int *ldab, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_zpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZPBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpbtrs)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zpbtrs_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                          ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZPBTRS", ret);
}

// ========================================================================== //
//                 Positive definite tridiagonal matrices (PT)                //
// ========================================================================== //

/******************************************************************************
 * ZPTCON( N, D, E, ANORM, RCOND, RWORK, INFO )
 ******************************************************************************/
#define ZPTCON_TEST LAPACK_GLOBAL_SUFFIX(zptcon_test, ZPTCON_TEST)
void ZPTCON_TEST(const lapack_int *n, const double *d,
                 const lapack_complex_double *e, const double *anorm,
                 double *rcond, double *rwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zptcon)(*n, d, e, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_zptcon_work)(*n, d, e, *anorm, rcond, rwork);
#endif

    *info = lapacke_test_info_unshifted("ZPTCON", ret);
}

/******************************************************************************
 * ZPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK,
 * INFO )
 ******************************************************************************/
#define ZPTRFS_TEST LAPACK_GLOBAL_SUFFIX(zptrfs_test, ZPTRFS_TEST)
void ZPTRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *d, const lapack_complex_double *e,
                 const double *df, const lapack_complex_double *ef,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *ferr,
                 double *berr, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZPTRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zptrfs)(layout, *uplo, *n, *nrhs, d, e, df, ef,
                                     b_r, ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zptrfs_work)(layout, *uplo, *n, *nrhs, d, e, df,
                                          ef, b_r, ldb_r, x_r, ldx_r, ferr,
                                          berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZPTRFS", ret);
}

/******************************************************************************
 * ZPTSV( N, NRHS, D, E, B, LDB, INFO )
 ******************************************************************************/
#define ZPTSV_TEST LAPACK_GLOBAL_SUFFIX(zptsv_test, ZPTSV_TEST)
void ZPTSV_TEST(const lapack_int *n, const lapack_int *nrhs, double *d,
                lapack_complex_double *e, lapack_complex_double *b,
                const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("ZPTSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zptsv)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zptsv_work)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZPTSV", ret);
}

/******************************************************************************
 * ZPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define ZPTSVX_TEST LAPACK_GLOBAL_SUFFIX(zptsvx_test, ZPTSVX_TEST)
void ZPTSVX_TEST(const char *fact, const lapack_int *n, const lapack_int *nrhs,
                 const double *d, const lapack_complex_double *e, double *df,
                 lapack_complex_double *ef, const lapack_complex_double *b,
                 const lapack_int *ldb, lapack_complex_double *x,
                 const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZPTSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zptsvx)(layout, *fact, *n, *nrhs, d, e, df, ef,
                                     b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zptsvx_work)(layout, *fact, *n, *nrhs, d, e, df,
                                          ef, b_r, ldb_r, x_r, ldx_r, rcond,
                                          ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZPTSVX", ret);
}

/******************************************************************************
 * ZPTTRF( N, D, E, INFO )
 ******************************************************************************/
#define ZPTTRF_TEST LAPACK_GLOBAL_SUFFIX(zpttrf_test, ZPTTRF_TEST)
void ZPTTRF_TEST(const lapack_int *n, double *d, lapack_complex_double *e,
                 lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpttrf)(*n, d, e);
#else
    ret = API_SUFFIX(LAPACKE_zpttrf_work)(*n, d, e);
#endif

    *info = lapacke_test_info_unshifted("ZPTTRF", ret);
}

/******************************************************************************
 * ZPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )
 ******************************************************************************/
#define ZPTTRS_TEST LAPACK_GLOBAL_SUFFIX(zpttrs_test, ZPTTRS_TEST)
void ZPTTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *d, const lapack_complex_double *e,
                 lapack_complex_double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("ZPTTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpttrs)(layout, *uplo, *n, *nrhs, d, e, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zpttrs_work)(layout, *uplo, *n, *nrhs, d, e, b_r,
                                          ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZPTTRS", ret);
}

// ========================================================================== //
//                     Symmetric indefinite matrices (SY)                     //
// ========================================================================== //

/******************************************************************************
 * ZSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define ZSYCON_TEST LAPACK_GLOBAL_SUFFIX(zsycon_test, ZSYCON_TEST)
void ZSYCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_int *ipiv, const double *anorm, double *rcond,
                 lapack_complex_double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZSYCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsycon)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_zsycon_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          *anorm, rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("ZSYCON", ret);
}

/******************************************************************************
 * ZSYCON_3( UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define ZSYCON_3_TEST LAPACK_GLOBAL_SUFFIX(zsycon_3_test, ZSYCON_3_TEST)
void ZSYCON_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_complex_double *a, const lapack_int *lda,
                   const lapack_complex_double *e, const lapack_int *ipiv,
                   const double *anorm, double *rcond,
                   lapack_complex_double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZSYCON_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsycon_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv,
                                       *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_zsycon_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, *anorm, rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("ZSYCON_3", ret);
}

/******************************************************************************
 * ZSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define ZSYRFS_TEST LAPACK_GLOBAL_SUFFIX(zsyrfs_test, ZSYRFS_TEST)
void ZSYRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *af, const lapack_int *ldaf,
                 const lapack_int *ipiv, const lapack_complex_double *b,
                 const lapack_int *ldb, lapack_complex_double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZSYRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsyrfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_zsyrfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZSYRFS", ret);
}

/******************************************************************************
 * ZSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZSYSV_TEST LAPACK_GLOBAL_SUFFIX(zsysv_test, ZSYSV_TEST)
void ZSYSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_double *a, const lapack_int *lda,
                lapack_int *ipiv, lapack_complex_double *b,
                const lapack_int *ldb, lapack_complex_double *work,
                const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zsysv_work)(LAPACK_COL_MAJOR, *uplo, *n, *nrhs,
                                             a, *lda, ipiv, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("ZSYSV", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZSYSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsysv)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zsysv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZSYSV", ret);
}

/******************************************************************************
 * ZSYSV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZSYSV_RK_TEST LAPACK_GLOBAL_SUFFIX(zsysv_rk_test, ZSYSV_RK_TEST)
void ZSYSV_RK_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, lapack_complex_double *a,
                   const lapack_int *lda, lapack_complex_double *e,
                   lapack_int *ipiv, lapack_complex_double *b,
                   const lapack_int *ldb, lapack_complex_double *work,
                   const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zsysv_rk_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, e, ipiv, b,
                                                *ldb, work, *lwork);
        *info = lapacke_test_info("ZSYSV_RK", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZSYSV_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsysv_rk)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zsysv_rk_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r, work,
                                            *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZSYSV_RK", ret);
}

/******************************************************************************
 * ZSYSV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZSYSV_ROOK_TEST LAPACK_GLOBAL_SUFFIX(zsysv_rook_test, ZSYSV_ROOK_TEST)
void ZSYSV_ROOK_TEST(const char *uplo, const lapack_int *n,
                     const lapack_int *nrhs, lapack_complex_double *a,
                     const lapack_int *lda, lapack_int *ipiv,
                     lapack_complex_double *b, const lapack_int *ldb,
                     lapack_complex_double *work, const lapack_int *lwork,
                     lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                     ,
                     FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zsysv_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                  *nrhs, a, *lda, ipiv, b, *ldb,
                                                  work, *lwork);
        *info = lapacke_test_info("ZSYSV_ROOK", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZSYSV_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsysv_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zsysv_rook_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZSYSV_ROOK", ret);
}

/******************************************************************************
 * ZSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND,
 * FERR, BERR, WORK, LWORK, RWORK, INFO )
 ******************************************************************************/
#define ZSYSVX_TEST LAPACK_GLOBAL_SUFFIX(zsysvx_test, ZSYSVX_TEST)
void ZSYSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_double *a,
                 const lapack_int *lda, lapack_complex_double *af,
                 const lapack_int *ldaf, lapack_int *ipiv,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *rcond,
                 double *ferr, double *berr, lapack_complex_double *work,
                 const lapack_int *lwork, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zsysvx_work)(
            LAPACK_COL_MAJOR, *fact, *uplo, *n, *nrhs, a, *lda, af, *ldaf, ipiv,
            b, *ldb, x, *ldx, rcond, ferr, berr, work, *lwork, rwork);
        *info = lapacke_test_info("ZSYSVX", ret);
        return;
    }

    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZSYSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsysvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                     ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zsysvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, b_r,
        ldb_r, x_r, ldx_r, rcond, ferr, berr, work, *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZSYSVX", ret);
}

/******************************************************************************
 * ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZSYTRF_TEST LAPACK_GLOBAL_SUFFIX(zsytrf_test, ZSYTRF_TEST)
void ZSYTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *ipiv, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zsytrf_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                              *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("ZSYTRF", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZSYTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytrf)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zsytrf_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZSYTRF", ret);
}

/******************************************************************************
 * ZSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZSYTRF_RK_TEST LAPACK_GLOBAL_SUFFIX(zsytrf_rk_test, ZSYTRF_RK_TEST)
void ZSYTRF_RK_TEST(const char *uplo, const lapack_int *n,
                    lapack_complex_double *a, const lapack_int *lda,
                    lapack_complex_double *e, lapack_int *ipiv,
                    lapack_complex_double *work, const lapack_int *lwork,
                    lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zsytrf_rk_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("ZSYTRF_RK", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZSYTRF_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytrf_rk)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zsytrf_rk_work)(layout, *uplo, *n, a_r, lda_r, e,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZSYTRF_RK", ret);
}

/******************************************************************************
 * ZSYTRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZSYTRF_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(zsytrf_rook_test, ZSYTRF_ROOK_TEST)
void ZSYTRF_ROOK_TEST(const char *uplo, const lapack_int *n,
                      lapack_complex_double *a, const lapack_int *lda,
                      lapack_int *ipiv, lapack_complex_double *work,
                      const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zsytrf_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                   a, *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("ZSYTRF_ROOK", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZSYTRF_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytrf_rook)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zsytrf_rook_work)(layout, *uplo, *n, a_r, lda_r,
                                               ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZSYTRF_ROOK", ret);
}

/******************************************************************************
 * ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
 ******************************************************************************/
#define ZSYTRI_TEST LAPACK_GLOBAL_SUFFIX(zsytri_test, ZSYTRI_TEST)
void ZSYTRI_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_double *work,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZSYTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytri)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zsytri_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZSYTRI", ret);
}

/******************************************************************************
 * ZSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZSYTRI2_TEST LAPACK_GLOBAL_SUFFIX(zsytri2_test, ZSYTRI2_TEST)
void ZSYTRI2_TEST(const char *uplo, const lapack_int *n,
                  lapack_complex_double *a, const lapack_int *lda,
                  const lapack_int *ipiv, lapack_complex_double *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zsytri2_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                               *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("ZSYTRI2", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZSYTRI2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytri2)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zsytri2_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZSYTRI2", ret);
}

/******************************************************************************
 * ZSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
 ******************************************************************************/
#define ZSYTRI2X_TEST LAPACK_GLOBAL_SUFFIX(zsytri2x_test, ZSYTRI2X_TEST)
void ZSYTRI2X_TEST(const char *uplo, const lapack_int *n,
                   lapack_complex_double *a, const lapack_int *lda,
                   const lapack_int *ipiv, lapack_complex_double *work,
                   const lapack_int *nb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZSYTRI2X", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytri2x)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                       *nb);
#else
    ret = API_SUFFIX(LAPACKE_zsytri2x_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                            work, *nb);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZSYTRI2X", ret);
}

/******************************************************************************
 * ZSYTRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZSYTRI_3_TEST LAPACK_GLOBAL_SUFFIX(zsytri_3_test, ZSYTRI_3_TEST)
void ZSYTRI_3_TEST(const char *uplo, const lapack_int *n,
                   lapack_complex_double *a, const lapack_int *lda,
                   const lapack_complex_double *e, const lapack_int *ipiv,
                   lapack_complex_double *work, const lapack_int *lwork,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zsytri_3_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("ZSYTRI_3", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZSYTRI_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytri_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zsytri_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZSYTRI_3", ret);
}

/******************************************************************************
 * ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZSYTRS_TEST LAPACK_GLOBAL_SUFFIX(zsytrs_test, ZSYTRS_TEST)
void ZSYTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZSYTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                     b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zsytrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZSYTRS", ret);
}

/******************************************************************************
 * ZSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )
 ******************************************************************************/
#define ZSYTRS2_TEST LAPACK_GLOBAL_SUFFIX(zsytrs2_test, ZSYTRS2_TEST)
void ZSYTRS2_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                  const lapack_complex_double *a, const lapack_int *lda,
                  const lapack_int *ipiv, lapack_complex_double *b,
                  const lapack_int *ldb, lapack_complex_double *work,
                  lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZSYTRS2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytrs2)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                      ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zsytrs2_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                           ipiv, b_r, ldb_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZSYTRS2", ret);
}

/******************************************************************************
 * ZSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZSYTRS_3_TEST LAPACK_GLOBAL_SUFFIX(zsytrs_3_test, ZSYTRS_3_TEST)
void ZSYTRS_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, const lapack_complex_double *a,
                   const lapack_int *lda, const lapack_complex_double *e,
                   const lapack_int *ipiv, lapack_complex_double *b,
                   const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZSYTRS_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytrs_3)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zsytrs_3_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZSYTRS_3", ret);
}

/******************************************************************************
 * ZSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZSYTRS_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(zsytrs_rook_test, ZSYTRS_ROOK_TEST)
void ZSYTRS_ROOK_TEST(const char *uplo, const lapack_int *n,
                      const lapack_int *nrhs, const lapack_complex_double *a,
                      const lapack_int *lda, const lapack_int *ipiv,
                      lapack_complex_double *b, const lapack_int *ldb,
                      lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZSYTRS_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsytrs_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zsytrs_rook_work)(layout, *uplo, *n, *nrhs, a_r,
                                               lda_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZSYTRS_ROOK", ret);
}

// ========================================================================== //
//                  Packed symmetric indefinite matrices (SP)                 //
// ========================================================================== //

/******************************************************************************
 * ZSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define ZSPCON_TEST LAPACK_GLOBAL_SUFFIX(zspcon_test, ZSPCON_TEST)
void ZSPCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_double *ap, const lapack_int *ipiv,
                 const double *anorm, double *rcond,
                 lapack_complex_double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zsp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZSPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zspcon)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_zspcon_work)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                          rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("ZSPCON", ret);
}

/******************************************************************************
 * ZSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define ZSPRFS_TEST LAPACK_GLOBAL_SUFFIX(zsprfs_test, ZSPRFS_TEST)
void ZSPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *ap,
                 const lapack_complex_double *afp, const lapack_int *ipiv,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *ferr,
                 double *berr, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
    const lapack_complex_double *ap_r = ap;
    const lapack_complex_double *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_zsp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_zsp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("ZSPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                     ipiv, b_r, ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zsprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                          berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("ZSPRFS", ret);
}

/******************************************************************************
 * ZSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZSPSV_TEST LAPACK_GLOBAL_SUFFIX(zspsv_test, ZSPSV_TEST)
void ZSPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_double *ap, lapack_int *ipiv,
                lapack_complex_double *b, const lapack_int *ldb,
                lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_zsp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("ZSPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zspsv)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zspsv_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_zsp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("ZSPSV", ret);
}

/******************************************************************************
 * ZSPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZSPSVX_TEST LAPACK_GLOBAL_SUFFIX(zspsvx_test, ZSPSVX_TEST)
void ZSPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_double *ap,
                 lapack_complex_double *afp, lapack_int *ipiv,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *rcond,
                 double *ferr, double *berr, lapack_complex_double *work,
                 double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_double *ap_r = ap;
    lapack_complex_double *afp_r = afp;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_zsp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_zsp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZSPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zspsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, ipiv, b_r, ldb_r, x_r, ldx_r, rcond,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zspsvx_work)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                          afp_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                          rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZSPSVX", ret);
}

/******************************************************************************
 * ZSPTRF( UPLO, N, AP, IPIV, INFO )
 ******************************************************************************/
#define ZSPTRF_TEST LAPACK_GLOBAL_SUFFIX(zsptrf_test, ZSPTRF_TEST)
void ZSPTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *ap, lapack_int *ipiv, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zsp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZSPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsptrf)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zsptrf_work)(layout, *uplo, *n, ap_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("ZSPTRF", ret);
}

/******************************************************************************
 * ZSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
 ******************************************************************************/
#define ZSPTRI_TEST LAPACK_GLOBAL_SUFFIX(zsptri_test, ZSPTRI_TEST)
void ZSPTRI_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *ap, const lapack_int *ipiv,
                 lapack_complex_double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zsp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZSPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsptri)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zsptri_work)(layout, *uplo, *n, ap_r, ipiv, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zsp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("ZSPTRI", ret);
}

/******************************************************************************
 * ZSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZSPTRS_TEST LAPACK_GLOBAL_SUFFIX(zsptrs_test, ZSPTRS_TEST)
void ZSPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *ap, const lapack_int *ipiv,
                 lapack_complex_double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_zsp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("ZSPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zsptrs)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zsptrs_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("ZSPTRS", ret);
}

// ========================================================================== //
//                     Hermitian indefinite matrices (HE)                     //
// ========================================================================== //

/******************************************************************************
 * ZHECON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define ZHECON_TEST LAPACK_GLOBAL_SUFFIX(zhecon_test, ZHECON_TEST)
void ZHECON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_int *ipiv, const double *anorm, double *rcond,
                 lapack_complex_double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHECON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhecon)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_zhecon_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          *anorm, rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("ZHECON", ret);
}

/******************************************************************************
 * ZHECON_3( UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define ZHECON_3_TEST LAPACK_GLOBAL_SUFFIX(zhecon_3_test, ZHECON_3_TEST)
void ZHECON_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_complex_double *a, const lapack_int *lda,
                   const lapack_complex_double *e, const lapack_int *ipiv,
                   const double *anorm, double *rcond,
                   lapack_complex_double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHECON_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhecon_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv,
                                       *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_zhecon_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, *anorm, rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("ZHECON_3", ret);
}

/******************************************************************************
 * ZHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define ZHERFS_TEST LAPACK_GLOBAL_SUFFIX(zherfs_test, ZHERFS_TEST)
void ZHERFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *af, const lapack_int *ldaf,
                 const lapack_int *ipiv, const lapack_complex_double *b,
                 const lapack_int *ldb, lapack_complex_double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZHERFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zherfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_zherfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZHERFS", ret);
}

/******************************************************************************
 * ZHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHESV_TEST LAPACK_GLOBAL_SUFFIX(zhesv_test, ZHESV_TEST)
void ZHESV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_double *a, const lapack_int *lda,
                lapack_int *ipiv, lapack_complex_double *b,
                const lapack_int *ldb, lapack_complex_double *work,
                const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhesv_work)(LAPACK_COL_MAJOR, *uplo, *n, *nrhs,
                                             a, *lda, ipiv, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("ZHESV", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZHESV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhesv)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhesv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZHESV", ret);
}

/******************************************************************************
 * ZHESV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHESV_AA_TEST LAPACK_GLOBAL_SUFFIX(zhesv_aa_test, ZHESV_AA_TEST)
void ZHESV_AA_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, lapack_complex_double *a,
                   const lapack_int *lda, lapack_int *ipiv,
                   lapack_complex_double *b, const lapack_int *ldb,
                   lapack_complex_double *work, const lapack_int *lwork,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhesv_aa_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, ipiv, b, *ldb,
                                                work, *lwork);
        *info = lapacke_test_info("ZHESV_AA", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZHESV_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhesv_aa)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhesv_aa_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZHESV_AA", ret);
}

/******************************************************************************
 * ZHESV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHESV_RK_TEST LAPACK_GLOBAL_SUFFIX(zhesv_rk_test, ZHESV_RK_TEST)
void ZHESV_RK_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, lapack_complex_double *a,
                   const lapack_int *lda, lapack_complex_double *e,
                   lapack_int *ipiv, lapack_complex_double *b,
                   const lapack_int *ldb, lapack_complex_double *work,
                   const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhesv_rk_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, e, ipiv, b,
                                                *ldb, work, *lwork);
        *info = lapacke_test_info("ZHESV_RK", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZHESV_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhesv_rk)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhesv_rk_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r, work,
                                            *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZHESV_RK", ret);
}

/******************************************************************************
 * ZHESVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND,
 * FERR, BERR, WORK, LWORK, RWORK, INFO )
 ******************************************************************************/
#define ZHESVX_TEST LAPACK_GLOBAL_SUFFIX(zhesvx_test, ZHESVX_TEST)
void ZHESVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_double *a,
                 const lapack_int *lda, lapack_complex_double *af,
                 const lapack_int *ldaf, lapack_int *ipiv,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *rcond,
                 double *ferr, double *berr, lapack_complex_double *work,
                 const lapack_int *lwork, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhesvx_work)(
            LAPACK_COL_MAJOR, *fact, *uplo, *n, *nrhs, a, *lda, af, *ldaf, ipiv,
            b, *ldb, x, *ldx, rcond, ferr, berr, work, *lwork, rwork);
        *info = lapacke_test_info("ZHESVX", ret);
        return;
    }

    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZHESVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhesvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                     ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zhesvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, b_r,
        ldb_r, x_r, ldx_r, rcond, ferr, berr, work, *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZHESVX", ret);
}

/******************************************************************************
 * ZHETRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHETRF_TEST LAPACK_GLOBAL_SUFFIX(zhetrf_test, ZHETRF_TEST)
void ZHETRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *ipiv, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhetrf_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                              *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("ZHETRF", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetrf)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zhetrf_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZHETRF", ret);
}

/******************************************************************************
 * ZHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHETRF_AA_TEST LAPACK_GLOBAL_SUFFIX(zhetrf_aa_test, ZHETRF_AA_TEST)
void ZHETRF_AA_TEST(const char *uplo, const lapack_int *n,
                    lapack_complex_double *a, const lapack_int *lda,
                    lapack_int *ipiv, lapack_complex_double *work,
                    const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhetrf_aa_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("ZHETRF_AA", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHETRF_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetrf_aa)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zhetrf_aa_work)(layout, *uplo, *n, a_r, lda_r,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZHETRF_AA", ret);
}

/******************************************************************************
 * ZHETRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHETRF_RK_TEST LAPACK_GLOBAL_SUFFIX(zhetrf_rk_test, ZHETRF_RK_TEST)
void ZHETRF_RK_TEST(const char *uplo, const lapack_int *n,
                    lapack_complex_double *a, const lapack_int *lda,
                    lapack_complex_double *e, lapack_int *ipiv,
                    lapack_complex_double *work, const lapack_int *lwork,
                    lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhetrf_rk_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("ZHETRF_RK", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHETRF_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetrf_rk)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zhetrf_rk_work)(layout, *uplo, *n, a_r, lda_r, e,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZHETRF_RK", ret);
}

/******************************************************************************
 * ZHETRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHETRF_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(zhetrf_rook_test, ZHETRF_ROOK_TEST)
void ZHETRF_ROOK_TEST(const char *uplo, const lapack_int *n,
                      lapack_complex_double *a, const lapack_int *lda,
                      lapack_int *ipiv, lapack_complex_double *work,
                      const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhetrf_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                   a, *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("ZHETRF_ROOK", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHETRF_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetrf_rook)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zhetrf_rook_work)(layout, *uplo, *n, a_r, lda_r,
                                               ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZHETRF_ROOK", ret);
}

/******************************************************************************
 * ZHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
 ******************************************************************************/
#define ZHETRI_TEST LAPACK_GLOBAL_SUFFIX(zhetri_test, ZHETRI_TEST)
void ZHETRI_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_double *work,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetri)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zhetri_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZHETRI", ret);
}

/******************************************************************************
 * ZHETRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHETRI2_TEST LAPACK_GLOBAL_SUFFIX(zhetri2_test, ZHETRI2_TEST)
void ZHETRI2_TEST(const char *uplo, const lapack_int *n,
                  lapack_complex_double *a, const lapack_int *lda,
                  const lapack_int *ipiv, lapack_complex_double *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhetri2_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                               *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("ZHETRI2", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHETRI2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetri2)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zhetri2_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZHETRI2", ret);
}

/******************************************************************************
 * ZHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
 ******************************************************************************/
#define ZHETRI2X_TEST LAPACK_GLOBAL_SUFFIX(zhetri2x_test, ZHETRI2X_TEST)
void ZHETRI2X_TEST(const char *uplo, const lapack_int *n,
                   lapack_complex_double *a, const lapack_int *lda,
                   const lapack_int *ipiv, lapack_complex_double *work,
                   const lapack_int *nb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*lda, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHETRI2X", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetri2x)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                       *nb);
#else
    ret = API_SUFFIX(LAPACKE_zhetri2x_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                            work, *nb);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*lda, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZHETRI2X", ret);
}

/******************************************************************************
 * ZHETRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHETRI_3_TEST LAPACK_GLOBAL_SUFFIX(zhetri_3_test, ZHETRI_3_TEST)
void ZHETRI_3_TEST(const char *uplo, const lapack_int *n,
                   lapack_complex_double *a, const lapack_int *lda,
                   const lapack_complex_double *e, const lapack_int *ipiv,
                   lapack_complex_double *work, const lapack_int *lwork,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhetri_3_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("ZHETRI_3", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZHETRI_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetri_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zhetri_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhe_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZHETRI_3", ret);
}

/******************************************************************************
 * ZHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZHETRS_TEST LAPACK_GLOBAL_SUFFIX(zhetrs_test, ZHETRS_TEST)
void ZHETRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZHETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                     b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhetrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZHETRS", ret);
}

/******************************************************************************
 * ZHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )
 ******************************************************************************/
#define ZHETRS2_TEST LAPACK_GLOBAL_SUFFIX(zhetrs2_test, ZHETRS2_TEST)
void ZHETRS2_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                  const lapack_complex_double *a, const lapack_int *lda,
                  const lapack_int *ipiv, lapack_complex_double *b,
                  const lapack_int *ldb, lapack_complex_double *work,
                  lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZHETRS2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetrs2)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                      ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhetrs2_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                           ipiv, b_r, ldb_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZHETRS2", ret);
}

/******************************************************************************
 * ZHETRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZHETRS_3_TEST LAPACK_GLOBAL_SUFFIX(zhetrs_3_test, ZHETRS_3_TEST)
void ZHETRS_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, const lapack_complex_double *a,
                   const lapack_int *lda, const lapack_complex_double *e,
                   const lapack_int *ipiv, lapack_complex_double *b,
                   const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZHETRS_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetrs_3)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhetrs_3_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZHETRS_3", ret);
}

/******************************************************************************
 * ZHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZHETRS_AA_TEST LAPACK_GLOBAL_SUFFIX(zhetrs_aa_test, ZHETRS_AA_TEST)
void ZHETRS_AA_TEST(const char *uplo, const lapack_int *n,
                    const lapack_int *nrhs, const lapack_complex_double *a,
                    const lapack_int *lda, const lapack_int *ipiv,
                    lapack_complex_double *b, const lapack_int *ldb,
                    lapack_complex_double *work, const lapack_int *lwork,
                    lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zhetrs_aa_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                 *nrhs, a, *lda, ipiv, b, *ldb,
                                                 work, *lwork);
        *info = lapacke_test_info("ZHETRS_AA", ret);
        return;
    }

    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZHETRS_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetrs_aa)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                        ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhetrs_aa_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZHETRS_AA", ret);
}

/******************************************************************************
 * ZHETRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZHETRS_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(zhetrs_rook_test, ZHETRS_ROOK_TEST)
void ZHETRS_ROOK_TEST(const char *uplo, const lapack_int *n,
                      const lapack_int *nrhs, const lapack_complex_double *a,
                      const lapack_int *lda, const lapack_int *ipiv,
                      lapack_complex_double *b, const lapack_int *ldb,
                      lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZHETRS_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhetrs_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhetrs_rook_work)(layout, *uplo, *n, *nrhs, a_r,
                                               lda_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZHETRS_ROOK", ret);
}

// ========================================================================== //
//                  Packed Hermitian indefinite matrices (HP)                 //
// ========================================================================== //

/******************************************************************************
 * ZHPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define ZHPCON_TEST LAPACK_GLOBAL_SUFFIX(zhpcon_test, ZHPCON_TEST)
void ZHPCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_double *ap, const lapack_int *ipiv,
                 const double *anorm, double *rcond,
                 lapack_complex_double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zhp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZHPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhpcon)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_zhpcon_work)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                          rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("ZHPCON", ret);
}

/******************************************************************************
 * ZHPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define ZHPRFS_TEST LAPACK_GLOBAL_SUFFIX(zhprfs_test, ZHPRFS_TEST)
void ZHPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *ap,
                 const lapack_complex_double *afp, const lapack_int *ipiv,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *ferr,
                 double *berr, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
    const lapack_complex_double *ap_r = ap;
    const lapack_complex_double *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_zhp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_zhp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("ZHPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                     ipiv, b_r, ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zhprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                          berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("ZHPRFS", ret);
}

/******************************************************************************
 * ZHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZHPSV_TEST LAPACK_GLOBAL_SUFFIX(zhpsv_test, ZHPSV_TEST)
void ZHPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_double *ap, lapack_int *ipiv,
                lapack_complex_double *b, const lapack_int *ldb,
                lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_zhp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("ZHPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhpsv)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhpsv_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_zhp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("ZHPSV", ret);
}

/******************************************************************************
 * ZHPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZHPSVX_TEST LAPACK_GLOBAL_SUFFIX(zhpsvx_test, ZHPSVX_TEST)
void ZHPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_double *ap,
                 lapack_complex_double *afp, lapack_int *ipiv,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 lapack_complex_double *x, const lapack_int *ldx, double *rcond,
                 double *ferr, double *berr, lapack_complex_double *work,
                 double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_double *ap_r = ap;
    lapack_complex_double *afp_r = afp;
    lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_zhp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_zhp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("ZHPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhpsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, ipiv, b_r, ldb_r, x_r, ldx_r, rcond,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_zhpsvx_work)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                          afp_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                          rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_zge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("ZHPSVX", ret);
}

/******************************************************************************
 * ZHPTRF( UPLO, N, AP, IPIV, INFO )
 ******************************************************************************/
#define ZHPTRF_TEST LAPACK_GLOBAL_SUFFIX(zhptrf_test, ZHPTRF_TEST)
void ZHPTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *ap, lapack_int *ipiv, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zhp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZHPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhptrf)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zhptrf_work)(layout, *uplo, *n, ap_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("ZHPTRF", ret);
}

/******************************************************************************
 * ZHPTRI( UPLO, N, AP, IPIV, WORK, INFO )
 ******************************************************************************/
#define ZHPTRI_TEST LAPACK_GLOBAL_SUFFIX(zhptri_test, ZHPTRI_TEST)
void ZHPTRI_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *ap, const lapack_int *ipiv,
                 lapack_complex_double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_zhp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("ZHPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhptri)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zhptri_work)(layout, *uplo, *n, ap_r, ipiv, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zhp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("ZHPTRI", ret);
}

/******************************************************************************
 * ZHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define ZHPTRS_TEST LAPACK_GLOBAL_SUFFIX(zhptrs_test, ZHPTRS_TEST)
void ZHPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *ap, const lapack_int *ipiv,
                 lapack_complex_double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_zhp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("ZHPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zhptrs)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zhptrs_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("ZHPTRS", ret);
}

// ========================================================================== //
//                          Triangular matrices (TR)                          //
// ========================================================================== //

/******************************************************************************
 * ZTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZTRCON_TEST LAPACK_GLOBAL_SUFFIX(ztrcon_test, ZTRCON_TEST)
void ZTRCON_TEST(const char *norm, const char *uplo, const char *diag,
                 const lapack_int *n, const lapack_complex_double *a,
                 const lapack_int *lda, double *rcond,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ztr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZTRCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ztrcon)(layout, *norm, *uplo, *diag, *n, a_r,
                                     lda_r, rcond);
#else
    ret = API_SUFFIX(LAPACKE_ztrcon_work)(layout, *norm, *uplo, *diag, *n, a_r,
                                          lda_r, rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("ZTRCON", ret);
}

/******************************************************************************
 * ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, LDX, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define ZTRRFS_TEST LAPACK_GLOBAL_SUFFIX(ztrrfs_test, ZTRRFS_TEST)
void ZTRRFS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *b, const lapack_int *ldb,
                 const lapack_complex_double *x, const lapack_int *ldx,
                 double *ferr, double *berr, lapack_complex_double *work,
                 double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ztr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)x_r);
        lapacke_test_report_alloc_failure("ZTRRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ztrrfs)(layout, *uplo, *trans, *diag, *n, *nrhs,
                                     a_r, lda_r, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_ztrrfs_work)(layout, *uplo, *trans, *diag, *n,
                                          *nrhs, a_r, lda_r, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)x_r);
#endif
    *info = lapacke_test_info("ZTRRFS", ret);
}

/******************************************************************************
 * ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
 ******************************************************************************/
#define ZTRTRI_TEST LAPACK_GLOBAL_SUFFIX(ztrtri_test, ZTRTRI_TEST)
void ZTRTRI_TEST(const char *uplo, const char *diag, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ztr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZTRTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ztrtri)(layout, *uplo, *diag, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_ztrtri_work)(layout, *uplo, *diag, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ztr_rm_to_cm(*uplo, *diag, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZTRTRI", ret);
}

/******************************************************************************
 * ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define ZTRTRS_TEST LAPACK_GLOBAL_SUFFIX(ztrtrs_test, ZTRTRS_TEST)
void ZTRTRS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ztr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZTRTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ztrtrs)(layout, *uplo, *trans, *diag, *n, *nrhs,
                                     a_r, lda_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ztrtrs_work)(layout, *uplo, *trans, *diag, *n,
                                          *nrhs, a_r, lda_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZTRTRS", ret);
}

// ========================================================================== //
//                        Triangular band matrices (TB)                       //
// ========================================================================== //

/******************************************************************************
 * ZTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define ZTBCON_TEST LAPACK_GLOBAL_SUFFIX(ztbcon_test, ZTBCON_TEST)
void ZTBCON_TEST(const char *norm, const char *uplo, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_complex_double *ab, const lapack_int *ldab,
                 double *rcond, lapack_complex_double *work, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_ztb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("ZTBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ztbcon)(layout, *norm, *uplo, *diag, *n, *kd, ab_r,
                                     ldab_r, rcond);
#else
    ret = API_SUFFIX(LAPACKE_ztbcon_work)(layout, *norm, *uplo, *diag, *n, *kd,
                                          ab_r, ldab_r, rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("ZTBCON", ret);
}

/******************************************************************************
 * ZTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define ZTBRFS_TEST LAPACK_GLOBAL_SUFFIX(ztbrfs_test, ZTBRFS_TEST)
void ZTBRFS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const lapack_complex_double *ab,
                 const lapack_int *ldab, const lapack_complex_double *b,
                 const lapack_int *ldb, const lapack_complex_double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 lapack_complex_double *work, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_ztb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)x_r);
        lapacke_test_report_alloc_failure("ZTBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ztbrfs)(layout, *uplo, *trans, *diag, *n, *kd,
                                     *nrhs, ab_r, ldab_r, b_r, ldb_r, x_r,
                                     ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_ztbrfs_work)(layout, *uplo, *trans, *diag, *n, *kd,
                                          *nrhs, ab_r, ldab_r, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)x_r);
#endif
    *info = lapacke_test_info("ZTBRFS", ret);
}

/******************************************************************************
 * ZTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define ZTBTRS_TEST LAPACK_GLOBAL_SUFFIX(ztbtrs_test, ZTBTRS_TEST)
void ZTBTRS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const lapack_complex_double *ab,
                 const lapack_int *ldab, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_ztb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZTBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ztbtrs)(layout, *uplo, *trans, *diag, *n, *kd,
                                     *nrhs, ab_r, ldab_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ztbtrs_work)(layout, *uplo, *trans, *diag, *n, *kd,
                                          *nrhs, ab_r, ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZTBTRS", ret);
}

// ========================================================================== //
//                   Orthogonal factorizations (QR/LQ/QL/RQ)                  //
// ========================================================================== //

/******************************************************************************
 * ZGELQ2( M, N, A, LDA, TAU, WORK, INFO )
 ******************************************************************************/
#define ZGELQ2_TEST LAPACK_GLOBAL_SUFFIX(zgelq2_test, ZGELQ2_TEST)
void ZGELQ2_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *tau, lapack_complex_double *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGELQ2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgelq2)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zgelq2_work)(layout, *m, *n, a_r, lda_r, tau,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGELQ2", ret);
}

/******************************************************************************
 * ZGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGELQF_TEST LAPACK_GLOBAL_SUFFIX(zgelqf_test, ZGELQF_TEST)
void ZGELQF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgelqf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("ZGELQF", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGELQF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgelqf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zgelqf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGELQF", ret);
}

/******************************************************************************
 * ZGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGEQLF_TEST LAPACK_GLOBAL_SUFFIX(zgeqlf_test, ZGEQLF_TEST)
void ZGEQLF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgeqlf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("ZGEQLF", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGEQLF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqlf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zgeqlf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGEQLF", ret);
}

/******************************************************************************
 * ZGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGEQR_TEST LAPACK_GLOBAL_SUFFIX(zgeqr_test, ZGEQR_TEST)
void ZGEQR_TEST(const lapack_int *m, const lapack_int *n,
                lapack_complex_double *a, const lapack_int *lda,
                lapack_complex_double *t, const lapack_int *tsize,
                lapack_complex_double *work, const lapack_int *lwork,
                lapack_int *info)
{
    lapack_int ret = 0;
    if (*tsize == -1 || *lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgeqr_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                             t, *tsize, work, *lwork);
        *info = lapacke_test_info("ZGEQR", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGEQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqr)(layout, *m, *n, a_r, lda_r, t, *tsize);
#else
    ret = API_SUFFIX(LAPACKE_zgeqr_work)(layout, *m, *n, a_r, lda_r, t, *tsize,
                                         work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGEQR", ret);
}

/******************************************************************************
 * ZGEQR2( M, N, A, LDA, TAU, WORK, INFO )
 ******************************************************************************/
#define ZGEQR2_TEST LAPACK_GLOBAL_SUFFIX(zgeqr2_test, ZGEQR2_TEST)
void ZGEQR2_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *tau, lapack_complex_double *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGEQR2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqr2)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zgeqr2_work)(layout, *m, *n, a_r, lda_r, tau,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGEQR2", ret);
}

/******************************************************************************
 * ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGEQRF_TEST LAPACK_GLOBAL_SUFFIX(zgeqrf_test, ZGEQRF_TEST)
void ZGEQRF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgeqrf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("ZGEQRF", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGEQRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqrf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zgeqrf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGEQRF", ret);
}

/******************************************************************************
 * ZGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGEQRFP_TEST LAPACK_GLOBAL_SUFFIX(zgeqrfp_test, ZGEQRFP_TEST)
void ZGEQRFP_TEST(const lapack_int *m, const lapack_int *n,
                  lapack_complex_double *a, const lapack_int *lda,
                  lapack_complex_double *tau, lapack_complex_double *work,
                  const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgeqrfp_work)(LAPACK_COL_MAJOR, *m, *n, a,
                                               *lda, tau, work, *lwork);
        *info = lapacke_test_info("ZGEQRFP", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGEQRFP", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqrfp)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zgeqrfp_work)(layout, *m, *n, a_r, lda_r, tau,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGEQRFP", ret);
}

/******************************************************************************
 * ZGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGERQF_TEST LAPACK_GLOBAL_SUFFIX(zgerqf_test, ZGERQF_TEST)
void ZGERQF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgerqf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("ZGERQF", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGERQF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgerqf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zgerqf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGERQF", ret);
}

/******************************************************************************
 * ZUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZUNGLQ_TEST LAPACK_GLOBAL_SUFFIX(zunglq_test, ZUNGLQ_TEST)
void ZUNGLQ_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zunglq_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("ZUNGLQ", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZUNGLQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zunglq)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zunglq_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZUNGLQ", ret);
}

/******************************************************************************
 * ZUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZUNGQL_TEST LAPACK_GLOBAL_SUFFIX(zungql_test, ZUNGQL_TEST)
void ZUNGQL_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zungql_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("ZUNGQL", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZUNGQL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zungql)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zungql_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZUNGQL", ret);
}

/******************************************************************************
 * ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZUNGQR_TEST LAPACK_GLOBAL_SUFFIX(zungqr_test, ZUNGQR_TEST)
void ZUNGQR_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zungqr_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("ZUNGQR", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZUNGQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zungqr)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zungqr_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZUNGQR", ret);
}

/******************************************************************************
 * ZUNGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZUNGRQ_TEST LAPACK_GLOBAL_SUFFIX(zungrq_test, ZUNGRQ_TEST)
void ZUNGRQ_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zungrq_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("ZUNGRQ", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZUNGRQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zungrq)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zungrq_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZUNGRQ", ret);
}

/******************************************************************************
 * ZUNMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZUNMLQ_TEST LAPACK_GLOBAL_SUFFIX(zunmlq_test, ZUNMLQ_TEST)
void ZUNMLQ_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *tau, lapack_complex_double *c,
                 const lapack_int *ldc, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zunmlq_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("ZUNMLQ", ret);
        return;
    }

    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_zge_cm_to_rm(*k, nq, a, *lda, &lda_r);
    c_r = lapacke_test_zge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("ZUNMLQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zunmlq)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_zunmlq_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("ZUNMLQ", ret);
}

/******************************************************************************
 * ZUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZUNMQL_TEST LAPACK_GLOBAL_SUFFIX(zunmql_test, ZUNMQL_TEST)
void ZUNMQL_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *tau, lapack_complex_double *c,
                 const lapack_int *ldc, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zunmql_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("ZUNMQL", ret);
        return;
    }

    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_zge_cm_to_rm(nq, *k, a, *lda, &lda_r);
    c_r = lapacke_test_zge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("ZUNMQL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zunmql)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_zunmql_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("ZUNMQL", ret);
}

/******************************************************************************
 * ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZUNMQR_TEST LAPACK_GLOBAL_SUFFIX(zunmqr_test, ZUNMQR_TEST)
void ZUNMQR_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *tau, lapack_complex_double *c,
                 const lapack_int *ldc, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zunmqr_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("ZUNMQR", ret);
        return;
    }

    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_zge_cm_to_rm(nq, *k, a, *lda, &lda_r);
    c_r = lapacke_test_zge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("ZUNMQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zunmqr)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_zunmqr_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("ZUNMQR", ret);
}

// ========================================================================== //
//            Blocked and tall-skinny factorizations (QRT/TSQR/HR)            //
// ========================================================================== //

/******************************************************************************
 * ZGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
 ******************************************************************************/
#define ZGEQRT_TEST LAPACK_GLOBAL_SUFFIX(zgeqrt_test, ZGEQRT_TEST)
void ZGEQRT_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *nb,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *t, const lapack_int *ldt,
                 lapack_complex_double *work, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_zge_cm_to_rm(*nb, MIN(*m, *n), t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("ZGEQRT", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqrt)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                     ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_zgeqrt_work)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                          ldt_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*nb, MIN(*m, *n), t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("ZGEQRT", ret);
}

/******************************************************************************
 * ZGEQRT2( M, N, A, LDA, T, LDT, INFO )
 ******************************************************************************/
#define ZGEQRT2_TEST LAPACK_GLOBAL_SUFFIX(zgeqrt2_test, ZGEQRT2_TEST)
void ZGEQRT2_TEST(const lapack_int *m, const lapack_int *n,
                  lapack_complex_double *a, const lapack_int *lda,
                  lapack_complex_double *t, const lapack_int *ldt,
                  lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_zge_cm_to_rm(*n, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("ZGEQRT2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqrt2)(layout, *m, *n, a_r, lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_zgeqrt2_work)(layout, *m, *n, a_r, lda_r, t_r,
                                           ldt_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("ZGEQRT2", ret);
}

/******************************************************************************
 * ZGEQRT3( M, N, A, LDA, T, LDT, INFO )
 ******************************************************************************/
#define ZGEQRT3_TEST LAPACK_GLOBAL_SUFFIX(zgeqrt3_test, ZGEQRT3_TEST)
void ZGEQRT3_TEST(const lapack_int *m, const lapack_int *n,
                  lapack_complex_double *a, const lapack_int *lda,
                  lapack_complex_double *t, const lapack_int *ldt,
                  lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_zge_cm_to_rm(*n, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("ZGEQRT3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqrt3)(layout, *m, *n, a_r, lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_zgeqrt3_work)(layout, *m, *n, a_r, lda_r, t_r,
                                           ldt_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*n, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("ZGEQRT3", ret);
}

/******************************************************************************
 * ZGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGETSQRHRT_TEST LAPACK_GLOBAL_SUFFIX(zgetsqrhrt_test, ZGETSQRHRT_TEST)
void ZGETSQRHRT_TEST(const lapack_int *m, const lapack_int *n,
                     const lapack_int *mb1, const lapack_int *nb1,
                     const lapack_int *nb2, lapack_complex_double *a,
                     const lapack_int *lda, lapack_complex_double *t,
                     const lapack_int *ldt, lapack_complex_double *work,
                     const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgetsqrhrt_work)(LAPACK_COL_MAJOR, *m, *n,
                                                  *mb1, *nb1, *nb2, a, *lda, t,
                                                  *ldt, work, *lwork);
        *info = lapacke_test_info("ZGETSQRHRT", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_zge_cm_to_rm(*nb2, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("ZGETSQRHRT", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgetsqrhrt)(layout, *m, *n, *mb1, *nb1, *nb2, a_r,
                                         lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_zgetsqrhrt_work)(
        layout, *m, *n, *mb1, *nb1, *nb2, a_r, lda_r, t_r, ldt_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*nb2, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("ZGETSQRHRT", ret);
}

/******************************************************************************
 * ZUNHR_COL( M, N, NB, A, LDA, T, LDT, D, INFO )
 ******************************************************************************/
#define ZUNHR_COL_TEST LAPACK_GLOBAL_SUFFIX(zunhr_col_test, ZUNHR_COL_TEST)
void ZUNHR_COL_TEST(const lapack_int *m, const lapack_int *n,
                    const lapack_int *nb, lapack_complex_double *a,
                    const lapack_int *lda, lapack_complex_double *t,
                    const lapack_int *ldt, lapack_complex_double *d,
                    lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_zge_cm_to_rm(*ldt, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("ZUNHR_COL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zunhr_col)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                        ldt_r, d);
#else
    ret = API_SUFFIX(LAPACKE_zunhr_col_work)(layout, *m, *n, *nb, a_r, lda_r,
                                             t_r, ldt_r, d);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(*ldt, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("ZUNHR_COL", ret);
}

// ========================================================================== //
//                    Rank-revealing factorizations (QP/TZ)                   //
// ========================================================================== //

/******************************************************************************
 * ZGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, INFO )
 ******************************************************************************/
#define ZGEQP3_TEST LAPACK_GLOBAL_SUFFIX(zgeqp3_test, ZGEQP3_TEST)
void ZGEQP3_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *jpvt, lapack_complex_double *tau,
                 lapack_complex_double *work, const lapack_int *lwork,
                 double *rwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgeqp3_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              jpvt, tau, work, *lwork, rwork);
        *info = lapacke_test_info("ZGEQP3", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGEQP3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqp3)(layout, *m, *n, a_r, lda_r, jpvt, tau);
#else
    ret = API_SUFFIX(LAPACKE_zgeqp3_work)(layout, *m, *n, a_r, lda_r, jpvt, tau,
                                          work, *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGEQP3", ret);
}

/******************************************************************************
 * ZTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZTZRZF_TEST LAPACK_GLOBAL_SUFFIX(ztzrzf_test, ZTZRZF_TEST)
void ZTZRZF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ztzrzf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("ZTZRZF", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZTZRZF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ztzrzf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_ztzrzf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZTZRZF", ret);
}

// ========================================================================== //
//                             Least squares (LS)                             //
// ========================================================================== //

/******************************************************************************
 * ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGELS_TEST LAPACK_GLOBAL_SUFFIX(zgels_test, ZGELS_TEST)
void ZGELS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                const lapack_int *nrhs, lapack_complex_double *a,
                const lapack_int *lda, lapack_complex_double *b,
                const lapack_int *ldb, lapack_complex_double *work,
                const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgels_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                             *nrhs, a, *lda, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("ZGELS", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGELS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgels)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgels_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                         lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGELS", ret);
}

/******************************************************************************
 * ZGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, RWORK,
 * IWORK, INFO )
 ******************************************************************************/
#define ZGELSD_TEST LAPACK_GLOBAL_SUFFIX(zgelsd_test, ZGELSD_TEST)
void ZGELSD_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_double *a,
                 const lapack_int *lda, lapack_complex_double *b,
                 const lapack_int *ldb, double *s, const double *rcond,
                 lapack_int *rank, lapack_complex_double *work,
                 const lapack_int *lwork, double *rwork, lapack_int *iwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgelsd_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, s, *rcond, rank,
                                              work, *lwork, rwork, iwork);
        *info = lapacke_test_info("ZGELSD", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGELSD", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgelsd)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, s, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_zgelsd_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, s, *rcond, rank, work,
                                          *lwork, rwork, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGELSD", ret);
}

/******************************************************************************
 * ZGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, RWORK, INFO
 * )
 ******************************************************************************/
#define ZGELSS_TEST LAPACK_GLOBAL_SUFFIX(zgelss_test, ZGELSS_TEST)
void ZGELSS_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_double *a,
                 const lapack_int *lda, lapack_complex_double *b,
                 const lapack_int *ldb, double *s, const double *rcond,
                 lapack_int *rank, lapack_complex_double *work,
                 const lapack_int *lwork, double *rwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgelss_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, s, *rcond, rank,
                                              work, *lwork, rwork);
        *info = lapacke_test_info("ZGELSS", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGELSS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgelss)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, s, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_zgelss_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, s, *rcond, rank, work,
                                          *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGELSS", ret);
}

/******************************************************************************
 * ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, RWORK,
 * INFO )
 ******************************************************************************/
#define ZGELSY_TEST LAPACK_GLOBAL_SUFFIX(zgelsy_test, ZGELSY_TEST)
void ZGELSY_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_double *a,
                 const lapack_int *lda, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_int *jpvt, const double *rcond,
                 lapack_int *rank, lapack_complex_double *work,
                 const lapack_int *lwork, double *rwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgelsy_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, jpvt, *rcond,
                                              rank, work, *lwork, rwork);
        *info = lapacke_test_info("ZGELSY", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGELSY", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgelsy)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, jpvt, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_zgelsy_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, jpvt, *rcond, rank, work,
                                          *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGELSY", ret);
}

/******************************************************************************
 * ZGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define ZGETSLS_TEST LAPACK_GLOBAL_SUFFIX(zgetsls_test, ZGETSLS_TEST)
void ZGETSLS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                  const lapack_int *nrhs, lapack_complex_double *a,
                  const lapack_int *lda, lapack_complex_double *b,
                  const lapack_int *ldb, lapack_complex_double *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgetsls_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                               *nrhs, a, *lda, b, *ldb, work,
                                               *lwork);
        *info = lapacke_test_info("ZGETSLS", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGETSLS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgetsls)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                      b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgetsls_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                           lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGETSLS", ret);
}

// ========================================================================== //
//                             SVD and bidiagonal                             //
// ========================================================================== //

/******************************************************************************
 * ZBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, RWORK, INFO
 * )
 ******************************************************************************/
#define ZBDSQR_TEST LAPACK_GLOBAL_SUFFIX(zbdsqr_test, ZBDSQR_TEST)
void ZBDSQR_TEST(const char *uplo, const lapack_int *n, const lapack_int *ncvt,
                 const lapack_int *nru, const lapack_int *ncc, double *d,
                 double *e, lapack_complex_double *vt, const lapack_int *ldvt,
                 lapack_complex_double *u, const lapack_int *ldu,
                 lapack_complex_double *c, const lapack_int *ldc, double *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *vt_r = vt;
    lapack_int ldvt_r = *ldvt;
    lapack_complex_double *u_r = u;
    lapack_int ldu_r = *ldu;
    lapack_complex_double *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    vt_r = lapacke_test_zge_cm_to_rm(*n, *ncvt, vt, *ldvt, &ldvt_r);
    u_r = lapacke_test_zge_cm_to_rm(*nru, *n, u, *ldu, &ldu_r);
    c_r = lapacke_test_zge_cm_to_rm(*n, *ncc, c, *ldc, &ldc_r);
    if (vt_r == NULL || u_r == NULL || c_r == NULL) {
        LAPACKE_free(vt_r);
        LAPACKE_free(u_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("ZBDSQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zbdsqr)(layout, *uplo, *n, *ncvt, *nru, *ncc, d, e,
                                     vt_r, ldvt_r, u_r, ldu_r, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_zbdsqr_work)(layout, *uplo, *n, *ncvt, *nru, *ncc,
                                          d, e, vt_r, ldvt_r, u_r, ldu_r, c_r,
                                          ldc_r, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *ncvt, vt_r, ldvt_r, vt, *ldvt);
    lapacke_test_zge_rm_to_cm(*nru, *n, u_r, ldu_r, u, *ldu);
    lapacke_test_zge_rm_to_cm(*n, *ncc, c_r, ldc_r, c, *ldc);
    LAPACKE_free(vt_r);
    LAPACKE_free(u_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("ZBDSQR", ret);
}

/******************************************************************************
 * ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK,
 * INFO )
 ******************************************************************************/
#define ZGESVD_TEST LAPACK_GLOBAL_SUFFIX(zgesvd_test, ZGESVD_TEST)
void ZGESVD_TEST(const char *jobu, const char *jobvt, const lapack_int *m,
                 const lapack_int *n, lapack_complex_double *a,
                 const lapack_int *lda, double *s, lapack_complex_double *u,
                 const lapack_int *ldu, lapack_complex_double *vt,
                 const lapack_int *ldvt, lapack_complex_double *work,
                 const lapack_int *lwork, double *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN jobu_len, FORTRAN_STRLEN jobvt_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_zgesvd_work)(LAPACK_COL_MAJOR, *jobu, *jobvt,
                                              *m, *n, a, *lda, s, u, *ldu, vt,
                                              *ldvt, work, *lwork, rwork);
        *info = lapacke_test_info("ZGESVD", ret);
        return;
    }

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_double *u_r = u;
    lapack_int ldu_r = *ldu;
    lapack_complex_double *vt_r = vt;
    lapack_int ldvt_r = *ldvt;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int urows = (API_SUFFIX(LAPACKE_lsame)(*jobu, 'a') ||
                              API_SUFFIX(LAPACKE_lsame)(*jobu, 's'))
                                 ? *m
                                 : 1;
    const lapack_int ucols = API_SUFFIX(LAPACKE_lsame)(*jobu, 'a') ? *m
                             : API_SUFFIX(LAPACKE_lsame)(*jobu, 's')
                                 ? MIN(*m, *n)
                                 : 1;
    const lapack_int vtrows = API_SUFFIX(LAPACKE_lsame)(*jobvt, 'a') ? *n
                              : API_SUFFIX(LAPACKE_lsame)(*jobvt, 's')
                                  ? MIN(*m, *n)
                                  : 1;
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    u_r = lapacke_test_zge_cm_to_rm(urows, ucols, u, *ldu, &ldu_r);
    vt_r = lapacke_test_zge_cm_to_rm(vtrows, *n, vt, *ldvt, &ldvt_r);
    if (a_r == NULL || u_r == NULL || vt_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(u_r);
        LAPACKE_free(vt_r);
        lapacke_test_report_alloc_failure("ZGESVD", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgesvd)(layout, *jobu, *jobvt, *m, *n, a_r, lda_r,
                                     s, u_r, ldu_r, vt_r, ldvt_r, rwork);
#else
    ret = API_SUFFIX(LAPACKE_zgesvd_work)(layout, *jobu, *jobvt, *m, *n, a_r,
                                          lda_r, s, u_r, ldu_r, vt_r, ldvt_r,
                                          work, *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(urows, ucols, u_r, ldu_r, u, *ldu);
    lapacke_test_zge_rm_to_cm(vtrows, *n, vt_r, ldvt_r, vt, *ldvt);
    LAPACKE_free(a_r);
    LAPACKE_free(u_r);
    LAPACKE_free(vt_r);
#endif
    *info = lapacke_test_info("ZGESVD", ret);
}

// ========================================================================== //
//                             Auxiliary routines                             //
// ========================================================================== //

/******************************************************************************
 * ZLACGV( N, X, INCX )
 ******************************************************************************/
#define ZLACGV_TEST LAPACK_GLOBAL_SUFFIX(zlacgv_test, ZLACGV_TEST)
void ZLACGV_TEST(const lapack_int *n, lapack_complex_double *x,
                 const lapack_int *incx)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zlacgv)(*n, x, *incx);
#else
    ret = API_SUFFIX(LAPACKE_zlacgv_work)(*n, x, *incx);
#endif

    lapacke_test_info_unshifted("ZLACGV", ret);
}

/******************************************************************************
 * DOUBLE PRECISION FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
 ******************************************************************************/
#define ZLANGE_TEST LAPACK_GLOBAL_SUFFIX(zlange_test, ZLANGE_TEST)
double ZLANGE_TEST(const char *norm, const lapack_int *m, const lapack_int *n,
                   const lapack_complex_double *a, const lapack_int *lda,
                   double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len
#endif
)
{
    double res = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("ZLANGE", &info);
        return (double)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_zlange)(layout, *norm, *m, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_zlange_work)(layout, *norm, *m, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * DOUBLE PRECISION FUNCTION ZLANHE( NORM, UPLO, N, A, LDA, WORK )
 ******************************************************************************/
#define ZLANHE_TEST LAPACK_GLOBAL_SUFFIX(zlanhe_test, ZLANHE_TEST)
double ZLANHE_TEST(const char *norm, const char *uplo, const lapack_int *n,
                   const lapack_complex_double *a, const lapack_int *lda,
                   double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    double res = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zhe_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("ZLANHE", &info);
        return (double)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_zlanhe)(layout, *norm, *uplo, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_zlanhe_work)(layout, *norm, *uplo, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * DOUBLE PRECISION FUNCTION ZLANSY( NORM, UPLO, N, A, LDA, WORK )
 ******************************************************************************/
#define ZLANSY_TEST LAPACK_GLOBAL_SUFFIX(zlansy_test, ZLANSY_TEST)
double ZLANSY_TEST(const char *norm, const char *uplo, const lapack_int *n,
                   const lapack_complex_double *a, const lapack_int *lda,
                   double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    double res = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("ZLANSY", &info);
        return (double)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_zlansy)(layout, *norm, *uplo, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_zlansy_work)(layout, *norm, *uplo, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * DOUBLE PRECISION FUNCTION ZLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
 ******************************************************************************/
#define ZLANTR_TEST LAPACK_GLOBAL_SUFFIX(zlantr_test, ZLANTR_TEST)
double ZLANTR_TEST(const char *norm, const char *uplo, const char *diag,
                   const lapack_int *m, const lapack_int *n,
                   const lapack_complex_double *a, const lapack_int *lda,
                   double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                   FORTRAN_STRLEN diag_len
#endif
)
{
    double res = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("ZLANTR", &info);
        return (double)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_zlantr)(layout, *norm, *uplo, *diag, *m, *n, a_r,
                                     lda_r);
#else
    res = API_SUFFIX(LAPACKE_zlantr_work)(layout, *norm, *uplo, *diag, *m, *n,
                                          a_r, lda_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * ZLARNV( IDIST, ISEED, N, X )
 ******************************************************************************/
#define ZLARNV_TEST LAPACK_GLOBAL_SUFFIX(zlarnv_test, ZLARNV_TEST)
void ZLARNV_TEST(const lapack_int *idist, lapack_int *iseed,
                 const lapack_int *n, lapack_complex_double *x)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zlarnv)(*idist, iseed, *n, x);
#else
    ret = API_SUFFIX(LAPACKE_zlarnv_work)(*idist, iseed, *n, x);
#endif

    lapacke_test_info_unshifted("ZLARNV", ret);
}

/******************************************************************************
 * ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
 ******************************************************************************/
#define ZLASCL_TEST LAPACK_GLOBAL_SUFFIX(zlascl_test, ZLASCL_TEST)
void ZLASCL_TEST(const char *type, const lapack_int *kl, const lapack_int *ku,
                 const double *cfrom, const double *cto, const lapack_int *m,
                 const lapack_int *n, lapack_complex_double *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN type_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int arows = API_SUFFIX(LAPACKE_lsame)(*type, 'b')   ? *kl + 1
                             : API_SUFFIX(LAPACKE_lsame)(*type, 'q') ? *ku + 1
                             : API_SUFFIX(LAPACKE_lsame)(*type, 'z')
                                 ? 2 * *kl + *ku + 1
                                 : *m;
    a_r = lapacke_test_zge_cm_to_rm(arows, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZLASCL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zlascl)(layout, *type, *kl, *ku, *cfrom, *cto, *m,
                                     *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_zlascl_work)(layout, *type, *kl, *ku, *cfrom, *cto,
                                          *m, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(arows, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZLASCL", ret);
}

/******************************************************************************
 * ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
 ******************************************************************************/
#define ZLASET_TEST LAPACK_GLOBAL_SUFFIX(zlaset_test, ZLASET_TEST)
void ZLASET_TEST(const char *uplo, const lapack_int *m, const lapack_int *n,
                 const lapack_complex_double *alpha,
                 const lapack_complex_double *beta, lapack_complex_double *a,
                 const lapack_int *lda
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("ZLASET", &info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zlaset)(layout, *uplo, *m, *n, *alpha, *beta, a_r,
                                     lda_r);
#else
    ret = API_SUFFIX(LAPACKE_zlaset_work)(layout, *uplo, *m, *n, *alpha, *beta,
                                          a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    lapacke_test_info("ZLASET", ret);
}

/******************************************************************************
 * ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
 ******************************************************************************/
#define ZLASWP_TEST LAPACK_GLOBAL_SUFFIX(zlaswp_test, ZLASWP_TEST)
void ZLASWP_TEST(const lapack_int *n, lapack_complex_double *a,
                 const lapack_int *lda, const lapack_int *k1,
                 const lapack_int *k2, const lapack_int *ipiv,
                 const lapack_int *incx)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(
        lapacke_test_laswp_rows(*k1, *k2, ipiv, *incx), *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("ZLASWP", &info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zlaswp)(layout, *n, a_r, lda_r, *k1, *k2, ipiv,
                                     *incx);
#else
    ret = API_SUFFIX(LAPACKE_zlaswp_work)(layout, *n, a_r, lda_r, *k1, *k2,
                                          ipiv, *incx);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(lapacke_test_laswp_rows(*k1, *k2, ipiv, *incx),
                              *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    lapacke_test_info("ZLASWP", ret);
}

/******************************************************************************
 * ZROT( N, CX, INCX, CY, INCY, C, S )
 ******************************************************************************/
#define ZROT_TEST LAPACK_GLOBAL_SUFFIX(zrot_test, ZROT_TEST)
void ZROT_TEST(const lapack_int *n, lapack_complex_double *cx,
               const lapack_int *incx, lapack_complex_double *cy,
               const lapack_int *incy, const double *c,
               const lapack_complex_double *s)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zrot)(*n, cx, *incx, cy, *incy, *c, *s);
#else
    ret = API_SUFFIX(LAPACKE_zrot_work)(*n, cx, *incx, cy, *incy, *c, *s);
#endif

    lapacke_test_info_unshifted("ZROT", ret);
}
