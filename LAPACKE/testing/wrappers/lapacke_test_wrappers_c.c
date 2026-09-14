/******************************************************************************
 * LAPACKE test wrappers (single precision complex)
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
 * CGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define CGECON_TEST LAPACK_GLOBAL_SUFFIX(cgecon_test, CGECON_TEST)
void CGECON_TEST(const char *norm, const lapack_int *n,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const float *anorm, float *rcond, lapack_complex_float *work,
                 float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGECON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgecon)(layout, *norm, *n, a_r, lda_r, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_cgecon_work)(layout, *norm, *n, a_r, lda_r, *anorm,
                                          rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("CGECON", ret);
}

/******************************************************************************
 * CGECXX( FACT, USESD, M, N, SESEL_ROWS, SEL_DESEL_COLS, KMAXFREE, ABSTOL,
 * RELTOL, A, LDA, K, MAXC2NRMK, RELMAXC2NRMK, FNRMK, IPIV, JPIV, TAU, C, LDC,
 * QRC, LDQRC, X, LDX, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
 ******************************************************************************/
#define CGECXX_TEST LAPACK_GLOBAL_SUFFIX(cgecxx_test, CGECXX_TEST)
void CGECXX_TEST(const char *fact, const char *usesd, const lapack_int *m,
                 const lapack_int *n, const lapack_int *sesel_rows,
                 const lapack_int *sel_desel_cols, const lapack_int *kmaxfree,
                 const float *abstol, const float *reltol,
                 lapack_complex_float *a, const lapack_int *lda, lapack_int *k,
                 float *maxc2nrmk, float *relmaxc2nrmk, float *fnrmk,
                 lapack_int *ipiv, lapack_int *jpiv, lapack_complex_float *tau,
                 lapack_complex_float *c, const lapack_int *ldc,
                 lapack_complex_float *qrc, const lapack_int *ldqrc,
                 lapack_complex_float *x, const lapack_int *ldx,
                 lapack_complex_float *work, const lapack_int *lwork,
                 float *rwork, const lapack_int *lrwork, lapack_int *iwork,
                 const lapack_int *liwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN usesd_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgecxx_work)(
            LAPACK_COL_MAJOR, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
            (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a, *lda,
            k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c, *ldc, qrc,
            *ldqrc, x, *ldx, work, *lwork, rwork, *lrwork, iwork, *liwork);
        *info = lapacke_test_info("CGECXX", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *c_r = c;
    lapack_int ldc_r = *ldc;
    lapack_complex_float *qrc_r = qrc;
    lapack_int ldqrc_r = *ldqrc;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    c_r = lapacke_test_cge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    qrc_r = lapacke_test_cge_cm_to_rm(*m, MIN(*m, *n), qrc, *ldqrc, &ldqrc_r);
    x_r = lapacke_test_cge_cm_to_rm(*m, *n, x, *ldx, &ldx_r);
    if (a_r == NULL || c_r == NULL || qrc_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(c_r);
        LAPACKE_free(qrc_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CGECXX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgecxx)(
        layout, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
        (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a_r, lda_r,
        k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c_r, ldc_r, qrc_r,
        ldqrc_r, x_r, ldx_r);
#else
    ret = API_SUFFIX(LAPACKE_cgecxx_work)(
        layout, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
        (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a_r, lda_r,
        k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c_r, ldc_r, qrc_r,
        ldqrc_r, x_r, ldx_r, work, *lwork, rwork, *lrwork, iwork, *liwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    lapacke_test_cge_rm_to_cm(*m, MIN(*m, *n), qrc_r, ldqrc_r, qrc, *ldqrc);
    lapacke_test_cge_rm_to_cm(*m, *n, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(c_r);
    LAPACKE_free(qrc_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CGECXX", ret);
}

/******************************************************************************
 * CGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )
 ******************************************************************************/
#define CGEEQU_TEST LAPACK_GLOBAL_SUFFIX(cgeequ_test, CGEEQU_TEST)
void CGEEQU_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_complex_float *a, const lapack_int *lda, float *r,
                 float *c, float *rowcnd, float *colcnd, float *amax,
                 lapack_int *info)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGEEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeequ)(layout, *m, *n, a_r, lda_r, r, c, rowcnd,
                                     colcnd, amax);
#else
    ret = API_SUFFIX(LAPACKE_cgeequ_work)(layout, *m, *n, a_r, lda_r, r, c,
                                          rowcnd, colcnd, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("CGEEQU", ret);
}

/******************************************************************************
 * CGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define CGERFS_TEST LAPACK_GLOBAL_SUFFIX(cgerfs_test, CGERFS_TEST)
void CGERFS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *af, const lapack_int *ldaf,
                 const lapack_int *ipiv, const lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    af_r = lapacke_test_cge_cm_to_rm(*n, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CGERFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgerfs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     af_r, ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cgerfs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CGERFS", ret);
}

/******************************************************************************
 * CGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CGESV_TEST LAPACK_GLOBAL_SUFFIX(cgesv_test, CGESV_TEST)
void CGESV_TEST(const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_float *a, const lapack_int *lda,
                lapack_int *ipiv, lapack_complex_float *b,
                const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGESV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgesv)(layout, *n, *nrhs, a_r, lda_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgesv_work)(layout, *n, *nrhs, a_r, lda_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGESV", ret);
}

/******************************************************************************
 * CGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CGESVX_TEST LAPACK_GLOBAL_SUFFIX(cgesvx_test, CGESVX_TEST)
void CGESVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_float *a,
                 const lapack_int *lda, lapack_complex_float *af,
                 const lapack_int *ldaf, lapack_int *ipiv, char *equed,
                 float *r, float *c, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *rcond, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    af_r = lapacke_test_cge_cm_to_rm(*n, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(af_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CGESVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgesvx)(
        layout, *fact, *trans, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, equed,
        r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, rwork);
#else
    ret = API_SUFFIX(LAPACKE_cgesvx_work)(
        layout, *fact, *trans, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, equed,
        r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(af_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CGESVX", ret);
}

/******************************************************************************
 * CGETF2( M, N, A, LDA, IPIV, INFO )
 ******************************************************************************/
#define CGETF2_TEST LAPACK_GLOBAL_SUFFIX(cgetf2_test, CGETF2_TEST)
void CGETF2_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGETF2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgetf2)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_cgetf2_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGETF2", ret);
}

/******************************************************************************
 * CGETRF( M, N, A, LDA, IPIV, INFO )
 ******************************************************************************/
#define CGETRF_TEST LAPACK_GLOBAL_SUFFIX(cgetrf_test, CGETRF_TEST)
void CGETRF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgetrf)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_cgetrf_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGETRF", ret);
}

/******************************************************************************
 * CGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGETRI_TEST LAPACK_GLOBAL_SUFFIX(cgetri_test, CGETRI_TEST)
void CGETRI_TEST(const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, const lapack_int *ipiv,
                 lapack_complex_float *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgetri_work)(LAPACK_COL_MAJOR, *n, a, *lda,
                                              ipiv, work, *lwork);
        *info = lapacke_test_info("CGETRI", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgetri)(layout, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_cgetri_work)(layout, *n, a_r, lda_r, ipiv, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGETRI", ret);
}

/******************************************************************************
 * CGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CGETRS_TEST LAPACK_GLOBAL_SUFFIX(cgetrs_test, CGETRS_TEST)
void CGETRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgetrs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgetrs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGETRS", ret);
}

// ========================================================================== //
//                         General band matrices (GB)                         //
// ========================================================================== //

/******************************************************************************
 * CGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define CGBCON_TEST LAPACK_GLOBAL_SUFFIX(cgbcon_test, CGBCON_TEST)
void CGBCON_TEST(const char *norm, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_complex_float *ab,
                 const lapack_int *ldab, const lapack_int *ipiv,
                 const float *anorm, float *rcond, lapack_complex_float *work,
                 float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("CGBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgbcon)(layout, *norm, *n, *kl, *ku, ab_r, ldab_r,
                                     ipiv, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_cgbcon_work)(layout, *norm, *n, *kl, *ku, ab_r,
                                          ldab_r, ipiv, *anorm, rcond, work,
                                          rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("CGBCON", ret);
}

/******************************************************************************
 * CGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFO )
 ******************************************************************************/
#define CGBEQU_TEST LAPACK_GLOBAL_SUFFIX(cgbequ_test, CGBEQU_TEST)
void CGBEQU_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_complex_float *ab,
                 const lapack_int *ldab, float *r, float *c, float *rowcnd,
                 float *colcnd, float *amax, lapack_int *info)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cgb_cm_to_rm(*m, *n, *kl, *ku, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("CGBEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgbequ)(layout, *m, *n, *kl, *ku, ab_r, ldab_r, r,
                                     c, rowcnd, colcnd, amax);
#else
    ret = API_SUFFIX(LAPACKE_cgbequ_work)(layout, *m, *n, *kl, *ku, ab_r,
                                          ldab_r, r, c, rowcnd, colcnd, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("CGBEQU", ret);
}

/******************************************************************************
 * CGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX,
 * FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CGBRFS_TEST LAPACK_GLOBAL_SUFFIX(cgbrfs_test, CGBRFS_TEST)
void CGBRFS_TEST(const char *trans, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_int *nrhs,
                 const lapack_complex_float *ab, const lapack_int *ldab,
                 const lapack_complex_float *afb, const lapack_int *ldafb,
                 const lapack_int *ipiv, const lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const lapack_complex_float *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cgb_cm_to_rm(*n, *n, *kl, *ku, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_cgb_cm_to_rm(*n, *n, *kl, *kl + *ku, afb, *ldafb,
                                      &ldafb_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)afb_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CGBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgbrfs)(layout, *trans, *n, *kl, *ku, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, ipiv, b_r, ldb_r,
                                     x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cgbrfs_work)(
        layout, *trans, *n, *kl, *ku, *nrhs, ab_r, ldab_r, afb_r, ldafb_r, ipiv,
        b_r, ldb_r, x_r, ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)afb_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CGBRFS", ret);
}

/******************************************************************************
 * CGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CGBSV_TEST LAPACK_GLOBAL_SUFFIX(cgbsv_test, CGBSV_TEST)
void CGBSV_TEST(const lapack_int *n, const lapack_int *kl, const lapack_int *ku,
                const lapack_int *nrhs, lapack_complex_float *ab,
                const lapack_int *ldab, lapack_int *ipiv,
                lapack_complex_float *b, const lapack_int *ldb,
                lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGBSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgbsv)(layout, *n, *kl, *ku, *nrhs, ab_r, ldab_r,
                                    ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgbsv_work)(layout, *n, *kl, *ku, *nrhs, ab_r,
                                         ldab_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cgb_rm_to_cm(*n, *n, *kl, *kl + *ku, ab_r, ldab_r, ab, *ldab);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGBSV", ret);
}

/******************************************************************************
 * CGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, EQUED, R,
 * C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CGBSVX_TEST LAPACK_GLOBAL_SUFFIX(cgbsvx_test, CGBSVX_TEST)
void CGBSVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *kl, const lapack_int *ku,
                 const lapack_int *nrhs, lapack_complex_float *ab,
                 const lapack_int *ldab, lapack_complex_float *afb,
                 const lapack_int *ldafb, lapack_int *ipiv, char *equed,
                 float *r, float *c, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *rcond, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_float *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cgb_cm_to_rm(*n, *n, *kl, *ku, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_cgb_cm_to_rm(*n, *n, *kl, *kl + *ku, afb, *ldafb,
                                      &ldafb_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(afb_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CGBSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgbsvx)(layout, *fact, *trans, *n, *kl, *ku, *nrhs,
                                     ab_r, ldab_r, afb_r, ldafb_r, ipiv, equed,
                                     r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr,
                                     berr, rwork);
#else
    ret = API_SUFFIX(LAPACKE_cgbsvx_work)(
        layout, *fact, *trans, *n, *kl, *ku, *nrhs, ab_r, ldab_r, afb_r,
        ldafb_r, ipiv, equed, r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr,
        work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cgb_rm_to_cm(*n, *n, *kl, *ku, ab_r, ldab_r, ab, *ldab);
    lapacke_test_cgb_rm_to_cm(*n, *n, *kl, *kl + *ku, afb_r, ldafb_r, afb,
                              *ldafb);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(ab_r);
    LAPACKE_free(afb_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CGBSVX", ret);
}

/******************************************************************************
 * CGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
 ******************************************************************************/
#define CGBTRF_TEST LAPACK_GLOBAL_SUFFIX(cgbtrf_test, CGBTRF_TEST)
void CGBTRF_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, lapack_complex_float *ab,
                 const lapack_int *ldab, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cgb_cm_to_rm(*m, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("CGBTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgbtrf)(layout, *m, *n, *kl, *ku, ab_r, ldab_r,
                                     ipiv);
#else
    ret = API_SUFFIX(LAPACKE_cgbtrf_work)(layout, *m, *n, *kl, *ku, ab_r,
                                          ldab_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cgb_rm_to_cm(*m, *n, *kl, *kl + *ku, ab_r, ldab_r, ab, *ldab);
    LAPACKE_free(ab_r);
#endif
    *info = lapacke_test_info("CGBTRF", ret);
}

/******************************************************************************
 * CGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CGBTRS_TEST LAPACK_GLOBAL_SUFFIX(cgbtrs_test, CGBTRS_TEST)
void CGBTRS_TEST(const char *trans, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_int *nrhs,
                 const lapack_complex_float *ab, const lapack_int *ldab,
                 const lapack_int *ipiv, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgbtrs)(layout, *trans, *n, *kl, *ku, *nrhs, ab_r,
                                     ldab_r, ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgbtrs_work)(layout, *trans, *n, *kl, *ku, *nrhs,
                                          ab_r, ldab_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGBTRS", ret);
}

// ========================================================================== //
//                      General tridiagonal matrices (GT)                     //
// ========================================================================== //

/******************************************************************************
 * CGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define CGTCON_TEST LAPACK_GLOBAL_SUFFIX(cgtcon_test, CGTCON_TEST)
void CGTCON_TEST(const char *norm, const lapack_int *n,
                 const lapack_complex_float *dl, const lapack_complex_float *d,
                 const lapack_complex_float *du,
                 const lapack_complex_float *du2, const lapack_int *ipiv,
                 const float *anorm, float *rcond, lapack_complex_float *work,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgtcon)(*norm, *n, dl, d, du, du2, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_cgtcon_work)(*norm, *n, dl, d, du, du2, ipiv,
                                          *anorm, rcond, work);
#endif

    *info = lapacke_test_info_unshifted("CGTCON", ret);
}

/******************************************************************************
 * CGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX,
 * FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CGTRFS_TEST LAPACK_GLOBAL_SUFFIX(cgtrfs_test, CGTRFS_TEST)
void CGTRFS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *dl, const lapack_complex_float *d,
                 const lapack_complex_float *du,
                 const lapack_complex_float *dlf,
                 const lapack_complex_float *df,
                 const lapack_complex_float *duf,
                 const lapack_complex_float *du2, const lapack_int *ipiv,
                 const lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *ferr,
                 float *berr, lapack_complex_float *work, float *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CGTRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgtrfs)(layout, *trans, *n, *nrhs, dl, d, du, dlf,
                                     df, duf, du2, ipiv, b_r, ldb_r, x_r, ldx_r,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cgtrfs_work)(layout, *trans, *n, *nrhs, dl, d, du,
                                          dlf, df, duf, du2, ipiv, b_r, ldb_r,
                                          x_r, ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CGTRFS", ret);
}

/******************************************************************************
 * CGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
 ******************************************************************************/
#define CGTSV_TEST LAPACK_GLOBAL_SUFFIX(cgtsv_test, CGTSV_TEST)
void CGTSV_TEST(const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_float *dl, lapack_complex_float *d,
                lapack_complex_float *du, lapack_complex_float *b,
                const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("CGTSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgtsv)(layout, *n, *nrhs, dl, d, du, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgtsv_work)(layout, *n, *nrhs, dl, d, du, b_r,
                                         ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGTSV", ret);
}

/******************************************************************************
 * CGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CGTSVX_TEST LAPACK_GLOBAL_SUFFIX(cgtsvx_test, CGTSVX_TEST)
void CGTSVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_float *dl,
                 const lapack_complex_float *d, const lapack_complex_float *du,
                 lapack_complex_float *dlf, lapack_complex_float *df,
                 lapack_complex_float *duf, lapack_complex_float *du2,
                 lapack_int *ipiv, const lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *rcond, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CGTSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgtsvx)(layout, *fact, *trans, *n, *nrhs, dl, d,
                                     du, dlf, df, duf, du2, ipiv, b_r, ldb_r,
                                     x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cgtsvx_work)(
        layout, *fact, *trans, *n, *nrhs, dl, d, du, dlf, df, duf, du2, ipiv,
        b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CGTSVX", ret);
}

/******************************************************************************
 * CGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
 ******************************************************************************/
#define CGTTRF_TEST LAPACK_GLOBAL_SUFFIX(cgttrf_test, CGTTRF_TEST)
void CGTTRF_TEST(const lapack_int *n, lapack_complex_float *dl,
                 lapack_complex_float *d, lapack_complex_float *du,
                 lapack_complex_float *du2, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgttrf)(*n, dl, d, du, du2, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_cgttrf_work)(*n, dl, d, du, du2, ipiv);
#endif

    *info = lapacke_test_info_unshifted("CGTTRF", ret);
}

/******************************************************************************
 * CGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CGTTRS_TEST LAPACK_GLOBAL_SUFFIX(cgttrs_test, CGTTRS_TEST)
void CGTTRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *dl, const lapack_complex_float *d,
                 const lapack_complex_float *du,
                 const lapack_complex_float *du2, const lapack_int *ipiv,
                 lapack_complex_float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("CGTTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgttrs)(layout, *trans, *n, *nrhs, dl, d, du, du2,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgttrs_work)(layout, *trans, *n, *nrhs, dl, d, du,
                                          du2, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGTTRS", ret);
}

// ========================================================================== //
//             Symmetric/Hermitian positive definite matrices (PO)            //
// ========================================================================== //

/******************************************************************************
 * CPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define CPOCON_TEST LAPACK_GLOBAL_SUFFIX(cpocon_test, CPOCON_TEST)
void CPOCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const float *anorm, float *rcond, lapack_complex_float *work,
                 float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CPOCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpocon)(layout, *uplo, *n, a_r, lda_r, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_cpocon_work)(layout, *uplo, *n, a_r, lda_r, *anorm,
                                          rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("CPOCON", ret);
}

/******************************************************************************
 * CPOEQU( N, A, LDA, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define CPOEQU_TEST LAPACK_GLOBAL_SUFFIX(cpoequ_test, CPOEQU_TEST)
void CPOEQU_TEST(const lapack_int *n, const lapack_complex_float *a,
                 const lapack_int *lda, float *s, float *scond, float *amax,
                 lapack_int *info)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CPOEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpoequ)(layout, *n, a_r, lda_r, s, scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_cpoequ_work)(layout, *n, a_r, lda_r, s, scond,
                                          amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("CPOEQU", ret);
}

/******************************************************************************
 * CPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define CPORFS_TEST LAPACK_GLOBAL_SUFFIX(cporfs_test, CPORFS_TEST)
void CPORFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *af, const lapack_int *ldaf,
                 const lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *ferr,
                 float *berr, lapack_complex_float *work, float *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CPORFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cporfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_cporfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, b_r, ldb_r, x_r, ldx_r,
                                          ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CPORFS", ret);
}

/******************************************************************************
 * CPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define CPOSV_TEST LAPACK_GLOBAL_SUFFIX(cposv_test, CPOSV_TEST)
void CPOSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_float *a, const lapack_int *lda,
                lapack_complex_float *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CPOSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cposv)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cposv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CPOSV", ret);
}

/******************************************************************************
 * CPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S, B, LDB, X, LDX,
 * RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CPOSVX_TEST LAPACK_GLOBAL_SUFFIX(cposvx_test, CPOSVX_TEST)
void CPOSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_float *a,
                 const lapack_int *lda, lapack_complex_float *af,
                 const lapack_int *ldaf, char *equed, float *s,
                 lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *rcond,
                 float *ferr, float *berr, lapack_complex_float *work,
                 float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(af_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CPOSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cposvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, equed, s, b_r, ldb_r,
                                     x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cposvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, equed, s,
        b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_cpo_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(af_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CPOSVX", ret);
}

/******************************************************************************
 * CPOTRF( UPLO, N, A, LDA, INFO )
 ******************************************************************************/
#define CPOTRF_TEST LAPACK_GLOBAL_SUFFIX(cpotrf_test, CPOTRF_TEST)
void CPOTRF_TEST(const char *uplo, const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CPOTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpotrf)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_cpotrf_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CPOTRF", ret);
}

/******************************************************************************
 * CPOTRI( UPLO, N, A, LDA, INFO )
 ******************************************************************************/
#define CPOTRI_TEST LAPACK_GLOBAL_SUFFIX(cpotri_test, CPOTRI_TEST)
void CPOTRI_TEST(const char *uplo, const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CPOTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpotri)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_cpotri_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CPOTRI", ret);
}

/******************************************************************************
 * CPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define CPOTRS_TEST LAPACK_GLOBAL_SUFFIX(cpotrs_test, CPOTRS_TEST)
void CPOTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CPOTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpotrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cpotrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CPOTRS", ret);
}

/******************************************************************************
 * CPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
 ******************************************************************************/
#define CPSTRF_TEST LAPACK_GLOBAL_SUFFIX(cpstrf_test, CPSTRF_TEST)
void CPSTRF_TEST(const char *uplo, const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, lapack_int *piv, lapack_int *rank,
                 const float *tol, float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CPSTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpstrf)(layout, *uplo, *n, a_r, lda_r, piv, rank,
                                     *tol);
#else
    ret = API_SUFFIX(LAPACKE_cpstrf_work)(layout, *uplo, *n, a_r, lda_r, piv,
                                          rank, *tol, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CPSTRF", ret);
}

// ========================================================================== //
//                   Packed positive definite matrices (PP)                   //
// ========================================================================== //

/******************************************************************************
 * CPPCON( UPLO, N, AP, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define CPPCON_TEST LAPACK_GLOBAL_SUFFIX(cppcon_test, CPPCON_TEST)
void CPPCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_float *ap, const float *anorm,
                 float *rcond, lapack_complex_float *work, float *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CPPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cppcon)(layout, *uplo, *n, ap_r, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_cppcon_work)(layout, *uplo, *n, ap_r, *anorm,
                                          rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("CPPCON", ret);
}

/******************************************************************************
 * CPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define CPPEQU_TEST LAPACK_GLOBAL_SUFFIX(cppequ_test, CPPEQU_TEST)
void CPPEQU_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_float *ap, float *s, float *scond,
                 float *amax, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CPPEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cppequ)(layout, *uplo, *n, ap_r, s, scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_cppequ_work)(layout, *uplo, *n, ap_r, s, scond,
                                          amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("CPPEQU", ret);
}

/******************************************************************************
 * CPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO
 * )
 ******************************************************************************/
#define CPPRFS_TEST LAPACK_GLOBAL_SUFFIX(cpprfs_test, CPPRFS_TEST)
void CPPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *ap,
                 const lapack_complex_float *afp, const lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
    const lapack_complex_float *ap_r = ap;
    const lapack_complex_float *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("CPPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r, b_r,
                                     ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cpprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          b_r, ldb_r, x_r, ldx_r, ferr, berr,
                                          work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("CPPRFS", ret);
}

/******************************************************************************
 * CPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
 ******************************************************************************/
#define CPPSV_TEST LAPACK_GLOBAL_SUFFIX(cppsv_test, CPPSV_TEST)
void CPPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_float *ap, lapack_complex_float *b,
                const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("CPPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cppsv)(layout, *uplo, *n, *nrhs, ap_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cppsv_work)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                         ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_cpp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("CPPSV", ret);
}

/******************************************************************************
 * CPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CPPSVX_TEST LAPACK_GLOBAL_SUFFIX(cppsvx_test, CPPSVX_TEST)
void CPPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_float *ap,
                 lapack_complex_float *afp, char *equed, float *s,
                 lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *rcond,
                 float *ferr, float *berr, lapack_complex_float *work,
                 float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *ap_r = ap;
    lapack_complex_float *afp_r = afp;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CPPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cppsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, equed, s, b_r, ldb_r, x_r, ldx_r,
                                     rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cppsvx_work)(
        layout, *fact, *uplo, *n, *nrhs, ap_r, afp_r, equed, s, b_r, ldb_r, x_r,
        ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_cpp_rm_to_cm(*uplo, *n, ap_r, ap);
    lapacke_test_cpp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CPPSVX", ret);
}

/******************************************************************************
 * CPPTRF( UPLO, N, AP, INFO )
 ******************************************************************************/
#define CPPTRF_TEST LAPACK_GLOBAL_SUFFIX(cpptrf_test, CPPTRF_TEST)
void CPPTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_float *ap, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CPPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpptrf)(layout, *uplo, *n, ap_r);
#else
    ret = API_SUFFIX(LAPACKE_cpptrf_work)(layout, *uplo, *n, ap_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("CPPTRF", ret);
}

/******************************************************************************
 * CPPTRI( UPLO, N, AP, INFO )
 ******************************************************************************/
#define CPPTRI_TEST LAPACK_GLOBAL_SUFFIX(cpptri_test, CPPTRI_TEST)
void CPPTRI_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_float *ap, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CPPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpptri)(layout, *uplo, *n, ap_r);
#else
    ret = API_SUFFIX(LAPACKE_cpptri_work)(layout, *uplo, *n, ap_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("CPPTRI", ret);
}

/******************************************************************************
 * CPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
 ******************************************************************************/
#define CPPTRS_TEST LAPACK_GLOBAL_SUFFIX(cpptrs_test, CPPTRS_TEST)
void CPPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *ap, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_cpp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("CPPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpptrs)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cpptrs_work)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                          ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("CPPTRS", ret);
}

// ========================================================================== //
//                    Positive definite band matrices (PB)                    //
// ========================================================================== //

/******************************************************************************
 * CPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define CPBCON_TEST LAPACK_GLOBAL_SUFFIX(cpbcon_test, CPBCON_TEST)
void CPBCON_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_complex_float *ab, const lapack_int *ldab,
                 const float *anorm, float *rcond, lapack_complex_float *work,
                 float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("CPBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpbcon)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_cpbcon_work)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                          *anorm, rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("CPBCON", ret);
}

/******************************************************************************
 * CPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define CPBEQU_TEST LAPACK_GLOBAL_SUFFIX(cpbequ_test, CPBEQU_TEST)
void CPBEQU_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_complex_float *ab, const lapack_int *ldab,
                 float *s, float *scond, float *amax, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("CPBEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpbequ)(layout, *uplo, *n, *kd, ab_r, ldab_r, s,
                                     scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_cpbequ_work)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                          s, scond, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("CPBEQU", ret);
}

/******************************************************************************
 * CPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define CPBRFS_TEST LAPACK_GLOBAL_SUFFIX(cpbrfs_test, CPBRFS_TEST)
void CPBRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const lapack_complex_float *ab,
                 const lapack_int *ldab, const lapack_complex_float *afb,
                 const lapack_int *ldafb, const lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const lapack_complex_float *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_cpb_cm_to_rm(*uplo, *n, *kd, afb, *ldafb, &ldafb_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)afb_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CPBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpbrfs)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, b_r, ldb_r, x_r,
                                     ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cpbrfs_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                          ldab_r, afb_r, ldafb_r, b_r, ldb_r,
                                          x_r, ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)afb_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CPBRFS", ret);
}

/******************************************************************************
 * CPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define CPBSV_TEST LAPACK_GLOBAL_SUFFIX(cpbsv_test, CPBSV_TEST)
void CPBSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                const lapack_int *nrhs, lapack_complex_float *ab,
                const lapack_int *ldab, lapack_complex_float *b,
                const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CPBSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpbsv)(layout, *uplo, *n, *kd, *nrhs, ab_r, ldab_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cpbsv_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                         ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CPBSV", ret);
}

/******************************************************************************
 * CPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CPBSVX_TEST LAPACK_GLOBAL_SUFFIX(cpbsvx_test, CPBSVX_TEST)
void CPBSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *kd, const lapack_int *nrhs,
                 lapack_complex_float *ab, const lapack_int *ldab,
                 lapack_complex_float *afb, const lapack_int *ldafb,
                 char *equed, float *s, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *rcond, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_float *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_cpb_cm_to_rm(*uplo, *n, *kd, afb, *ldafb, &ldafb_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(afb_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CPBSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpbsvx)(layout, *fact, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, equed, s, b_r,
                                     ldb_r, x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cpbsvx_work)(
        layout, *fact, *uplo, *n, *kd, *nrhs, ab_r, ldab_r, afb_r, ldafb_r,
        equed, s, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    lapacke_test_cpb_rm_to_cm(*uplo, *n, *kd, afb_r, ldafb_r, afb, *ldafb);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(ab_r);
    LAPACKE_free(afb_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CPBSVX", ret);
}

/******************************************************************************
 * CPBTRF( UPLO, N, KD, AB, LDAB, INFO )
 ******************************************************************************/
#define CPBTRF_TEST LAPACK_GLOBAL_SUFFIX(cpbtrf_test, CPBTRF_TEST)
void CPBTRF_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 lapack_complex_float *ab, const lapack_int *ldab,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("CPBTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpbtrf)(layout, *uplo, *n, *kd, ab_r, ldab_r);
#else
    ret = API_SUFFIX(LAPACKE_cpbtrf_work)(layout, *uplo, *n, *kd, ab_r, ldab_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    LAPACKE_free(ab_r);
#endif
    *info = lapacke_test_info("CPBTRF", ret);
}

/******************************************************************************
 * CPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define CPBTRS_TEST LAPACK_GLOBAL_SUFFIX(cpbtrs_test, CPBTRS_TEST)
void CPBTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const lapack_complex_float *ab,
                 const lapack_int *ldab, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_cpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CPBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpbtrs)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cpbtrs_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                          ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CPBTRS", ret);
}

// ========================================================================== //
//                 Positive definite tridiagonal matrices (PT)                //
// ========================================================================== //

/******************************************************************************
 * CPTCON( N, D, E, ANORM, RCOND, RWORK, INFO )
 ******************************************************************************/
#define CPTCON_TEST LAPACK_GLOBAL_SUFFIX(cptcon_test, CPTCON_TEST)
void CPTCON_TEST(const lapack_int *n, const float *d,
                 const lapack_complex_float *e, const float *anorm,
                 float *rcond, float *rwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cptcon)(*n, d, e, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_cptcon_work)(*n, d, e, *anorm, rcond, rwork);
#endif

    *info = lapacke_test_info_unshifted("CPTCON", ret);
}

/******************************************************************************
 * CPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK,
 * INFO )
 ******************************************************************************/
#define CPTRFS_TEST LAPACK_GLOBAL_SUFFIX(cptrfs_test, CPTRFS_TEST)
void CPTRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *d, const lapack_complex_float *e, const float *df,
                 const lapack_complex_float *ef, const lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CPTRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cptrfs)(layout, *uplo, *n, *nrhs, d, e, df, ef,
                                     b_r, ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cptrfs_work)(layout, *uplo, *n, *nrhs, d, e, df,
                                          ef, b_r, ldb_r, x_r, ldx_r, ferr,
                                          berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CPTRFS", ret);
}

/******************************************************************************
 * CPTSV( N, NRHS, D, E, B, LDB, INFO )
 ******************************************************************************/
#define CPTSV_TEST LAPACK_GLOBAL_SUFFIX(cptsv_test, CPTSV_TEST)
void CPTSV_TEST(const lapack_int *n, const lapack_int *nrhs, float *d,
                lapack_complex_float *e, lapack_complex_float *b,
                const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("CPTSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cptsv)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cptsv_work)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CPTSV", ret);
}

/******************************************************************************
 * CPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define CPTSVX_TEST LAPACK_GLOBAL_SUFFIX(cptsvx_test, CPTSVX_TEST)
void CPTSVX_TEST(const char *fact, const lapack_int *n, const lapack_int *nrhs,
                 const float *d, const lapack_complex_float *e, float *df,
                 lapack_complex_float *ef, const lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *rcond, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CPTSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cptsvx)(layout, *fact, *n, *nrhs, d, e, df, ef,
                                     b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cptsvx_work)(layout, *fact, *n, *nrhs, d, e, df,
                                          ef, b_r, ldb_r, x_r, ldx_r, rcond,
                                          ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CPTSVX", ret);
}

/******************************************************************************
 * CPTTRF( N, D, E, INFO )
 ******************************************************************************/
#define CPTTRF_TEST LAPACK_GLOBAL_SUFFIX(cpttrf_test, CPTTRF_TEST)
void CPTTRF_TEST(const lapack_int *n, float *d, lapack_complex_float *e,
                 lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpttrf)(*n, d, e);
#else
    ret = API_SUFFIX(LAPACKE_cpttrf_work)(*n, d, e);
#endif

    *info = lapacke_test_info_unshifted("CPTTRF", ret);
}

/******************************************************************************
 * CPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )
 ******************************************************************************/
#define CPTTRS_TEST LAPACK_GLOBAL_SUFFIX(cpttrs_test, CPTTRS_TEST)
void CPTTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *d, const lapack_complex_float *e,
                 lapack_complex_float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("CPTTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpttrs)(layout, *uplo, *n, *nrhs, d, e, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cpttrs_work)(layout, *uplo, *n, *nrhs, d, e, b_r,
                                          ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CPTTRS", ret);
}

// ========================================================================== //
//                     Symmetric indefinite matrices (SY)                     //
// ========================================================================== //

/******************************************************************************
 * CSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define CSYCON_TEST LAPACK_GLOBAL_SUFFIX(csycon_test, CSYCON_TEST)
void CSYCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_int *ipiv, const float *anorm, float *rcond,
                 lapack_complex_float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CSYCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csycon)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_csycon_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          *anorm, rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("CSYCON", ret);
}

/******************************************************************************
 * CSYCON_3( UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define CSYCON_3_TEST LAPACK_GLOBAL_SUFFIX(csycon_3_test, CSYCON_3_TEST)
void CSYCON_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_complex_float *a, const lapack_int *lda,
                   const lapack_complex_float *e, const lapack_int *ipiv,
                   const float *anorm, float *rcond, lapack_complex_float *work,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CSYCON_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csycon_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv,
                                       *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_csycon_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, *anorm, rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("CSYCON_3", ret);
}

/******************************************************************************
 * CSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define CSYRFS_TEST LAPACK_GLOBAL_SUFFIX(csyrfs_test, CSYRFS_TEST)
void CSYRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *af, const lapack_int *ldaf,
                 const lapack_int *ipiv, const lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_csy_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CSYRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csyrfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_csyrfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CSYRFS", ret);
}

/******************************************************************************
 * CSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define CSYSV_TEST LAPACK_GLOBAL_SUFFIX(csysv_test, CSYSV_TEST)
void CSYSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_float *a, const lapack_int *lda,
                lapack_int *ipiv, lapack_complex_float *b,
                const lapack_int *ldb, lapack_complex_float *work,
                const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_csysv_work)(LAPACK_COL_MAJOR, *uplo, *n, *nrhs,
                                             a, *lda, ipiv, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("CSYSV", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CSYSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csysv)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_csysv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CSYSV", ret);
}

/******************************************************************************
 * CSYSV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define CSYSV_RK_TEST LAPACK_GLOBAL_SUFFIX(csysv_rk_test, CSYSV_RK_TEST)
void CSYSV_RK_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, lapack_complex_float *a,
                   const lapack_int *lda, lapack_complex_float *e,
                   lapack_int *ipiv, lapack_complex_float *b,
                   const lapack_int *ldb, lapack_complex_float *work,
                   const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_csysv_rk_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, e, ipiv, b,
                                                *ldb, work, *lwork);
        *info = lapacke_test_info("CSYSV_RK", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CSYSV_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csysv_rk)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_csysv_rk_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r, work,
                                            *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CSYSV_RK", ret);
}

/******************************************************************************
 * CSYSV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define CSYSV_ROOK_TEST LAPACK_GLOBAL_SUFFIX(csysv_rook_test, CSYSV_ROOK_TEST)
void CSYSV_ROOK_TEST(const char *uplo, const lapack_int *n,
                     const lapack_int *nrhs, lapack_complex_float *a,
                     const lapack_int *lda, lapack_int *ipiv,
                     lapack_complex_float *b, const lapack_int *ldb,
                     lapack_complex_float *work, const lapack_int *lwork,
                     lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                     ,
                     FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_csysv_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                  *nrhs, a, *lda, ipiv, b, *ldb,
                                                  work, *lwork);
        *info = lapacke_test_info("CSYSV_ROOK", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CSYSV_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csysv_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_csysv_rook_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CSYSV_ROOK", ret);
}

/******************************************************************************
 * CSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND,
 * FERR, BERR, WORK, LWORK, RWORK, INFO )
 ******************************************************************************/
#define CSYSVX_TEST LAPACK_GLOBAL_SUFFIX(csysvx_test, CSYSVX_TEST)
void CSYSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_float *a,
                 const lapack_int *lda, lapack_complex_float *af,
                 const lapack_int *ldaf, lapack_int *ipiv,
                 const lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *rcond,
                 float *ferr, float *berr, lapack_complex_float *work,
                 const lapack_int *lwork, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_csysvx_work)(
            LAPACK_COL_MAJOR, *fact, *uplo, *n, *nrhs, a, *lda, af, *ldaf, ipiv,
            b, *ldb, x, *ldx, rcond, ferr, berr, work, *lwork, rwork);
        *info = lapacke_test_info("CSYSVX", ret);
        return;
    }

    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_csy_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CSYSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csysvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                     ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_csysvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, b_r,
        ldb_r, x_r, ldx_r, rcond, ferr, berr, work, *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CSYSVX", ret);
}

/******************************************************************************
 * CSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CSYTRF_TEST LAPACK_GLOBAL_SUFFIX(csytrf_test, CSYTRF_TEST)
void CSYTRF_TEST(const char *uplo, const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, lapack_int *ipiv,
                 lapack_complex_float *work, const lapack_int *lwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_csytrf_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                              *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("CSYTRF", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CSYTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytrf)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_csytrf_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CSYTRF", ret);
}

/******************************************************************************
 * CSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CSYTRF_RK_TEST LAPACK_GLOBAL_SUFFIX(csytrf_rk_test, CSYTRF_RK_TEST)
void CSYTRF_RK_TEST(const char *uplo, const lapack_int *n,
                    lapack_complex_float *a, const lapack_int *lda,
                    lapack_complex_float *e, lapack_int *ipiv,
                    lapack_complex_float *work, const lapack_int *lwork,
                    lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_csytrf_rk_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("CSYTRF_RK", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CSYTRF_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytrf_rk)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_csytrf_rk_work)(layout, *uplo, *n, a_r, lda_r, e,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CSYTRF_RK", ret);
}

/******************************************************************************
 * CSYTRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CSYTRF_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(csytrf_rook_test, CSYTRF_ROOK_TEST)
void CSYTRF_ROOK_TEST(const char *uplo, const lapack_int *n,
                      lapack_complex_float *a, const lapack_int *lda,
                      lapack_int *ipiv, lapack_complex_float *work,
                      const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_csytrf_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                   a, *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("CSYTRF_ROOK", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CSYTRF_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytrf_rook)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_csytrf_rook_work)(layout, *uplo, *n, a_r, lda_r,
                                               ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CSYTRF_ROOK", ret);
}

/******************************************************************************
 * CSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
 ******************************************************************************/
#define CSYTRI_TEST LAPACK_GLOBAL_SUFFIX(csytri_test, CSYTRI_TEST)
void CSYTRI_TEST(const char *uplo, const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, const lapack_int *ipiv,
                 lapack_complex_float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CSYTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytri)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_csytri_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CSYTRI", ret);
}

/******************************************************************************
 * CSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CSYTRI2_TEST LAPACK_GLOBAL_SUFFIX(csytri2_test, CSYTRI2_TEST)
void CSYTRI2_TEST(const char *uplo, const lapack_int *n,
                  lapack_complex_float *a, const lapack_int *lda,
                  const lapack_int *ipiv, lapack_complex_float *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_csytri2_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                               *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("CSYTRI2", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CSYTRI2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytri2)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_csytri2_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CSYTRI2", ret);
}

/******************************************************************************
 * CSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
 ******************************************************************************/
#define CSYTRI2X_TEST LAPACK_GLOBAL_SUFFIX(csytri2x_test, CSYTRI2X_TEST)
void CSYTRI2X_TEST(const char *uplo, const lapack_int *n,
                   lapack_complex_float *a, const lapack_int *lda,
                   const lapack_int *ipiv, lapack_complex_float *work,
                   const lapack_int *nb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CSYTRI2X", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytri2x)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                       *nb);
#else
    ret = API_SUFFIX(LAPACKE_csytri2x_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                            work, *nb);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CSYTRI2X", ret);
}

/******************************************************************************
 * CSYTRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CSYTRI_3_TEST LAPACK_GLOBAL_SUFFIX(csytri_3_test, CSYTRI_3_TEST)
void CSYTRI_3_TEST(const char *uplo, const lapack_int *n,
                   lapack_complex_float *a, const lapack_int *lda,
                   const lapack_complex_float *e, const lapack_int *ipiv,
                   lapack_complex_float *work, const lapack_int *lwork,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_csytri_3_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("CSYTRI_3", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CSYTRI_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytri_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_csytri_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CSYTRI_3", ret);
}

/******************************************************************************
 * CSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CSYTRS_TEST LAPACK_GLOBAL_SUFFIX(csytrs_test, CSYTRS_TEST)
void CSYTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CSYTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                     b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_csytrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CSYTRS", ret);
}

/******************************************************************************
 * CSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )
 ******************************************************************************/
#define CSYTRS2_TEST LAPACK_GLOBAL_SUFFIX(csytrs2_test, CSYTRS2_TEST)
void CSYTRS2_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                  const lapack_complex_float *a, const lapack_int *lda,
                  const lapack_int *ipiv, lapack_complex_float *b,
                  const lapack_int *ldb, lapack_complex_float *work,
                  lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CSYTRS2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytrs2)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                      ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_csytrs2_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                           ipiv, b_r, ldb_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CSYTRS2", ret);
}

/******************************************************************************
 * CSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CSYTRS_3_TEST LAPACK_GLOBAL_SUFFIX(csytrs_3_test, CSYTRS_3_TEST)
void CSYTRS_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, const lapack_complex_float *a,
                   const lapack_int *lda, const lapack_complex_float *e,
                   const lapack_int *ipiv, lapack_complex_float *b,
                   const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CSYTRS_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytrs_3)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_csytrs_3_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CSYTRS_3", ret);
}

/******************************************************************************
 * CSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CSYTRS_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(csytrs_rook_test, CSYTRS_ROOK_TEST)
void CSYTRS_ROOK_TEST(const char *uplo, const lapack_int *n,
                      const lapack_int *nrhs, const lapack_complex_float *a,
                      const lapack_int *lda, const lapack_int *ipiv,
                      lapack_complex_float *b, const lapack_int *ldb,
                      lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CSYTRS_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csytrs_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_csytrs_rook_work)(layout, *uplo, *n, *nrhs, a_r,
                                               lda_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CSYTRS_ROOK", ret);
}

// ========================================================================== //
//                  Packed symmetric indefinite matrices (SP)                 //
// ========================================================================== //

/******************************************************************************
 * CSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define CSPCON_TEST LAPACK_GLOBAL_SUFFIX(cspcon_test, CSPCON_TEST)
void CSPCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_float *ap, const lapack_int *ipiv,
                 const float *anorm, float *rcond, lapack_complex_float *work,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_csp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CSPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cspcon)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_cspcon_work)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                          rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("CSPCON", ret);
}

/******************************************************************************
 * CSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define CSPRFS_TEST LAPACK_GLOBAL_SUFFIX(csprfs_test, CSPRFS_TEST)
void CSPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *ap,
                 const lapack_complex_float *afp, const lapack_int *ipiv,
                 const lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *ferr,
                 float *berr, lapack_complex_float *work, float *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
    const lapack_complex_float *ap_r = ap;
    const lapack_complex_float *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_csp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_csp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("CSPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                     ipiv, b_r, ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_csprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                          berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("CSPRFS", ret);
}

/******************************************************************************
 * CSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CSPSV_TEST LAPACK_GLOBAL_SUFFIX(cspsv_test, CSPSV_TEST)
void CSPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_float *ap, lapack_int *ipiv,
                lapack_complex_float *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_csp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("CSPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cspsv)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cspsv_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_csp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("CSPSV", ret);
}

/******************************************************************************
 * CSPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CSPSVX_TEST LAPACK_GLOBAL_SUFFIX(cspsvx_test, CSPSVX_TEST)
void CSPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_float *ap,
                 lapack_complex_float *afp, lapack_int *ipiv,
                 const lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *rcond,
                 float *ferr, float *berr, lapack_complex_float *work,
                 float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_float *ap_r = ap;
    lapack_complex_float *afp_r = afp;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_csp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_csp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CSPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cspsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, ipiv, b_r, ldb_r, x_r, ldx_r, rcond,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_cspsvx_work)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                          afp_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                          rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CSPSVX", ret);
}

/******************************************************************************
 * CSPTRF( UPLO, N, AP, IPIV, INFO )
 ******************************************************************************/
#define CSPTRF_TEST LAPACK_GLOBAL_SUFFIX(csptrf_test, CSPTRF_TEST)
void CSPTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_float *ap, lapack_int *ipiv, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_csp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CSPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csptrf)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_csptrf_work)(layout, *uplo, *n, ap_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("CSPTRF", ret);
}

/******************************************************************************
 * CSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
 ******************************************************************************/
#define CSPTRI_TEST LAPACK_GLOBAL_SUFFIX(csptri_test, CSPTRI_TEST)
void CSPTRI_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_float *ap, const lapack_int *ipiv,
                 lapack_complex_float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_csp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CSPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csptri)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_csptri_work)(layout, *uplo, *n, ap_r, ipiv, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_csp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("CSPTRI", ret);
}

/******************************************************************************
 * CSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CSPTRS_TEST LAPACK_GLOBAL_SUFFIX(csptrs_test, CSPTRS_TEST)
void CSPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *ap, const lapack_int *ipiv,
                 lapack_complex_float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_csp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("CSPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_csptrs)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_csptrs_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("CSPTRS", ret);
}

// ========================================================================== //
//                     Hermitian indefinite matrices (HE)                     //
// ========================================================================== //

/******************************************************************************
 * CHECON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define CHECON_TEST LAPACK_GLOBAL_SUFFIX(checon_test, CHECON_TEST)
void CHECON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_int *ipiv, const float *anorm, float *rcond,
                 lapack_complex_float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHECON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_checon)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_checon_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          *anorm, rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("CHECON", ret);
}

/******************************************************************************
 * CHECON_3( UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define CHECON_3_TEST LAPACK_GLOBAL_SUFFIX(checon_3_test, CHECON_3_TEST)
void CHECON_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_complex_float *a, const lapack_int *lda,
                   const lapack_complex_float *e, const lapack_int *ipiv,
                   const float *anorm, float *rcond, lapack_complex_float *work,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHECON_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_checon_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv,
                                       *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_checon_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, *anorm, rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("CHECON_3", ret);
}

/******************************************************************************
 * CHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define CHERFS_TEST LAPACK_GLOBAL_SUFFIX(cherfs_test, CHERFS_TEST)
void CHERFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *af, const lapack_int *ldaf,
                 const lapack_int *ipiv, const lapack_complex_float *b,
                 const lapack_int *ldb, lapack_complex_float *x,
                 const lapack_int *ldx, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_che_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CHERFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cherfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_cherfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CHERFS", ret);
}

/******************************************************************************
 * CHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHESV_TEST LAPACK_GLOBAL_SUFFIX(chesv_test, CHESV_TEST)
void CHESV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_float *a, const lapack_int *lda,
                lapack_int *ipiv, lapack_complex_float *b,
                const lapack_int *ldb, lapack_complex_float *work,
                const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chesv_work)(LAPACK_COL_MAJOR, *uplo, *n, *nrhs,
                                             a, *lda, ipiv, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("CHESV", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CHESV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chesv)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chesv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CHESV", ret);
}

/******************************************************************************
 * CHESV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHESV_AA_TEST LAPACK_GLOBAL_SUFFIX(chesv_aa_test, CHESV_AA_TEST)
void CHESV_AA_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, lapack_complex_float *a,
                   const lapack_int *lda, lapack_int *ipiv,
                   lapack_complex_float *b, const lapack_int *ldb,
                   lapack_complex_float *work, const lapack_int *lwork,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chesv_aa_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, ipiv, b, *ldb,
                                                work, *lwork);
        *info = lapacke_test_info("CHESV_AA", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CHESV_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chesv_aa)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chesv_aa_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CHESV_AA", ret);
}

/******************************************************************************
 * CHESV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHESV_RK_TEST LAPACK_GLOBAL_SUFFIX(chesv_rk_test, CHESV_RK_TEST)
void CHESV_RK_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, lapack_complex_float *a,
                   const lapack_int *lda, lapack_complex_float *e,
                   lapack_int *ipiv, lapack_complex_float *b,
                   const lapack_int *ldb, lapack_complex_float *work,
                   const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chesv_rk_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, e, ipiv, b,
                                                *ldb, work, *lwork);
        *info = lapacke_test_info("CHESV_RK", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CHESV_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chesv_rk)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chesv_rk_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r, work,
                                            *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CHESV_RK", ret);
}

/******************************************************************************
 * CHESVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND,
 * FERR, BERR, WORK, LWORK, RWORK, INFO )
 ******************************************************************************/
#define CHESVX_TEST LAPACK_GLOBAL_SUFFIX(chesvx_test, CHESVX_TEST)
void CHESVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_float *a,
                 const lapack_int *lda, lapack_complex_float *af,
                 const lapack_int *ldaf, lapack_int *ipiv,
                 const lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *rcond,
                 float *ferr, float *berr, lapack_complex_float *work,
                 const lapack_int *lwork, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chesvx_work)(
            LAPACK_COL_MAJOR, *fact, *uplo, *n, *nrhs, a, *lda, af, *ldaf, ipiv,
            b, *ldb, x, *ldx, rcond, ferr, berr, work, *lwork, rwork);
        *info = lapacke_test_info("CHESVX", ret);
        return;
    }

    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_che_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CHESVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chesvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                     ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_chesvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, b_r,
        ldb_r, x_r, ldx_r, rcond, ferr, berr, work, *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CHESVX", ret);
}

/******************************************************************************
 * CHETRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHETRF_TEST LAPACK_GLOBAL_SUFFIX(chetrf_test, CHETRF_TEST)
void CHETRF_TEST(const char *uplo, const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, lapack_int *ipiv,
                 lapack_complex_float *work, const lapack_int *lwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chetrf_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                              *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("CHETRF", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetrf)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_chetrf_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CHETRF", ret);
}

/******************************************************************************
 * CHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHETRF_AA_TEST LAPACK_GLOBAL_SUFFIX(chetrf_aa_test, CHETRF_AA_TEST)
void CHETRF_AA_TEST(const char *uplo, const lapack_int *n,
                    lapack_complex_float *a, const lapack_int *lda,
                    lapack_int *ipiv, lapack_complex_float *work,
                    const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chetrf_aa_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("CHETRF_AA", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHETRF_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetrf_aa)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_chetrf_aa_work)(layout, *uplo, *n, a_r, lda_r,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CHETRF_AA", ret);
}

/******************************************************************************
 * CHETRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHETRF_RK_TEST LAPACK_GLOBAL_SUFFIX(chetrf_rk_test, CHETRF_RK_TEST)
void CHETRF_RK_TEST(const char *uplo, const lapack_int *n,
                    lapack_complex_float *a, const lapack_int *lda,
                    lapack_complex_float *e, lapack_int *ipiv,
                    lapack_complex_float *work, const lapack_int *lwork,
                    lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chetrf_rk_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("CHETRF_RK", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHETRF_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetrf_rk)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_chetrf_rk_work)(layout, *uplo, *n, a_r, lda_r, e,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CHETRF_RK", ret);
}

/******************************************************************************
 * CHETRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHETRF_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(chetrf_rook_test, CHETRF_ROOK_TEST)
void CHETRF_ROOK_TEST(const char *uplo, const lapack_int *n,
                      lapack_complex_float *a, const lapack_int *lda,
                      lapack_int *ipiv, lapack_complex_float *work,
                      const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chetrf_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                   a, *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("CHETRF_ROOK", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHETRF_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetrf_rook)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_chetrf_rook_work)(layout, *uplo, *n, a_r, lda_r,
                                               ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CHETRF_ROOK", ret);
}

/******************************************************************************
 * CHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
 ******************************************************************************/
#define CHETRI_TEST LAPACK_GLOBAL_SUFFIX(chetri_test, CHETRI_TEST)
void CHETRI_TEST(const char *uplo, const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, const lapack_int *ipiv,
                 lapack_complex_float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetri)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_chetri_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CHETRI", ret);
}

/******************************************************************************
 * CHETRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHETRI2_TEST LAPACK_GLOBAL_SUFFIX(chetri2_test, CHETRI2_TEST)
void CHETRI2_TEST(const char *uplo, const lapack_int *n,
                  lapack_complex_float *a, const lapack_int *lda,
                  const lapack_int *ipiv, lapack_complex_float *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chetri2_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                               *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("CHETRI2", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHETRI2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetri2)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_chetri2_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CHETRI2", ret);
}

/******************************************************************************
 * CHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
 ******************************************************************************/
#define CHETRI2X_TEST LAPACK_GLOBAL_SUFFIX(chetri2x_test, CHETRI2X_TEST)
void CHETRI2X_TEST(const char *uplo, const lapack_int *n,
                   lapack_complex_float *a, const lapack_int *lda,
                   const lapack_int *ipiv, lapack_complex_float *work,
                   const lapack_int *nb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*lda, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHETRI2X", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetri2x)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                       *nb);
#else
    ret = API_SUFFIX(LAPACKE_chetri2x_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                            work, *nb);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*lda, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CHETRI2X", ret);
}

/******************************************************************************
 * CHETRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHETRI_3_TEST LAPACK_GLOBAL_SUFFIX(chetri_3_test, CHETRI_3_TEST)
void CHETRI_3_TEST(const char *uplo, const lapack_int *n,
                   lapack_complex_float *a, const lapack_int *lda,
                   const lapack_complex_float *e, const lapack_int *ipiv,
                   lapack_complex_float *work, const lapack_int *lwork,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chetri_3_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("CHETRI_3", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CHETRI_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetri_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_chetri_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_che_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CHETRI_3", ret);
}

/******************************************************************************
 * CHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CHETRS_TEST LAPACK_GLOBAL_SUFFIX(chetrs_test, CHETRS_TEST)
void CHETRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CHETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                     b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chetrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CHETRS", ret);
}

/******************************************************************************
 * CHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )
 ******************************************************************************/
#define CHETRS2_TEST LAPACK_GLOBAL_SUFFIX(chetrs2_test, CHETRS2_TEST)
void CHETRS2_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                  const lapack_complex_float *a, const lapack_int *lda,
                  const lapack_int *ipiv, lapack_complex_float *b,
                  const lapack_int *ldb, lapack_complex_float *work,
                  lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CHETRS2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetrs2)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                      ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chetrs2_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                           ipiv, b_r, ldb_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CHETRS2", ret);
}

/******************************************************************************
 * CHETRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CHETRS_3_TEST LAPACK_GLOBAL_SUFFIX(chetrs_3_test, CHETRS_3_TEST)
void CHETRS_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, const lapack_complex_float *a,
                   const lapack_int *lda, const lapack_complex_float *e,
                   const lapack_int *ipiv, lapack_complex_float *b,
                   const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CHETRS_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetrs_3)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chetrs_3_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CHETRS_3", ret);
}

/******************************************************************************
 * CHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define CHETRS_AA_TEST LAPACK_GLOBAL_SUFFIX(chetrs_aa_test, CHETRS_AA_TEST)
void CHETRS_AA_TEST(const char *uplo, const lapack_int *n,
                    const lapack_int *nrhs, const lapack_complex_float *a,
                    const lapack_int *lda, const lapack_int *ipiv,
                    lapack_complex_float *b, const lapack_int *ldb,
                    lapack_complex_float *work, const lapack_int *lwork,
                    lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_chetrs_aa_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                 *nrhs, a, *lda, ipiv, b, *ldb,
                                                 work, *lwork);
        *info = lapacke_test_info("CHETRS_AA", ret);
        return;
    }

    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CHETRS_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetrs_aa)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                        ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chetrs_aa_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CHETRS_AA", ret);
}

/******************************************************************************
 * CHETRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CHETRS_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(chetrs_rook_test, CHETRS_ROOK_TEST)
void CHETRS_ROOK_TEST(const char *uplo, const lapack_int *n,
                      const lapack_int *nrhs, const lapack_complex_float *a,
                      const lapack_int *lda, const lapack_int *ipiv,
                      lapack_complex_float *b, const lapack_int *ldb,
                      lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CHETRS_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chetrs_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chetrs_rook_work)(layout, *uplo, *n, *nrhs, a_r,
                                               lda_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CHETRS_ROOK", ret);
}

// ========================================================================== //
//                  Packed Hermitian indefinite matrices (HP)                 //
// ========================================================================== //

/******************************************************************************
 * CHPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define CHPCON_TEST LAPACK_GLOBAL_SUFFIX(chpcon_test, CHPCON_TEST)
void CHPCON_TEST(const char *uplo, const lapack_int *n,
                 const lapack_complex_float *ap, const lapack_int *ipiv,
                 const float *anorm, float *rcond, lapack_complex_float *work,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_chp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CHPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chpcon)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_chpcon_work)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                          rcond, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("CHPCON", ret);
}

/******************************************************************************
 * CHPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define CHPRFS_TEST LAPACK_GLOBAL_SUFFIX(chprfs_test, CHPRFS_TEST)
void CHPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *ap,
                 const lapack_complex_float *afp, const lapack_int *ipiv,
                 const lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *ferr,
                 float *berr, lapack_complex_float *work, float *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
    const lapack_complex_float *ap_r = ap;
    const lapack_complex_float *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_chp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_chp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("CHPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                     ipiv, b_r, ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_chprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                          berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("CHPRFS", ret);
}

/******************************************************************************
 * CHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CHPSV_TEST LAPACK_GLOBAL_SUFFIX(chpsv_test, CHPSV_TEST)
void CHPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                lapack_complex_float *ap, lapack_int *ipiv,
                lapack_complex_float *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_chp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("CHPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chpsv)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chpsv_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_chp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("CHPSV", ret);
}

/******************************************************************************
 * CHPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, RWORK, INFO )
 ******************************************************************************/
#define CHPSVX_TEST LAPACK_GLOBAL_SUFFIX(chpsvx_test, CHPSVX_TEST)
void CHPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const lapack_complex_float *ap,
                 lapack_complex_float *afp, lapack_int *ipiv,
                 const lapack_complex_float *b, const lapack_int *ldb,
                 lapack_complex_float *x, const lapack_int *ldx, float *rcond,
                 float *ferr, float *berr, lapack_complex_float *work,
                 float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_float *ap_r = ap;
    lapack_complex_float *afp_r = afp;
    lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_chp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_chp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("CHPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chpsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, ipiv, b_r, ldb_r, x_r, ldx_r, rcond,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_chpsvx_work)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                          afp_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                          rcond, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_chp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_cge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("CHPSVX", ret);
}

/******************************************************************************
 * CHPTRF( UPLO, N, AP, IPIV, INFO )
 ******************************************************************************/
#define CHPTRF_TEST LAPACK_GLOBAL_SUFFIX(chptrf_test, CHPTRF_TEST)
void CHPTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_float *ap, lapack_int *ipiv, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_chp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CHPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chptrf)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_chptrf_work)(layout, *uplo, *n, ap_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_chp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("CHPTRF", ret);
}

/******************************************************************************
 * CHPTRI( UPLO, N, AP, IPIV, WORK, INFO )
 ******************************************************************************/
#define CHPTRI_TEST LAPACK_GLOBAL_SUFFIX(chptri_test, CHPTRI_TEST)
void CHPTRI_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_float *ap, const lapack_int *ipiv,
                 lapack_complex_float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_chp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("CHPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chptri)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_chptri_work)(layout, *uplo, *n, ap_r, ipiv, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_chp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("CHPTRI", ret);
}

/******************************************************************************
 * CHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define CHPTRS_TEST LAPACK_GLOBAL_SUFFIX(chptrs_test, CHPTRS_TEST)
void CHPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *ap, const lapack_int *ipiv,
                 lapack_complex_float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_chp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("CHPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_chptrs)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_chptrs_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("CHPTRS", ret);
}

// ========================================================================== //
//                          Triangular matrices (TR)                          //
// ========================================================================== //

/******************************************************************************
 * CTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define CTRCON_TEST LAPACK_GLOBAL_SUFFIX(ctrcon_test, CTRCON_TEST)
void CTRCON_TEST(const char *norm, const char *uplo, const char *diag,
                 const lapack_int *n, const lapack_complex_float *a,
                 const lapack_int *lda, float *rcond,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ctr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CTRCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ctrcon)(layout, *norm, *uplo, *diag, *n, a_r,
                                     lda_r, rcond);
#else
    ret = API_SUFFIX(LAPACKE_ctrcon_work)(layout, *norm, *uplo, *diag, *n, a_r,
                                          lda_r, rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("CTRCON", ret);
}

/******************************************************************************
 * CTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, LDX, FERR, BERR, WORK,
 * RWORK, INFO )
 ******************************************************************************/
#define CTRRFS_TEST LAPACK_GLOBAL_SUFFIX(ctrrfs_test, CTRRFS_TEST)
void CTRRFS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *b, const lapack_int *ldb,
                 const lapack_complex_float *x, const lapack_int *ldx,
                 float *ferr, float *berr, lapack_complex_float *work,
                 float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ctr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)x_r);
        lapacke_test_report_alloc_failure("CTRRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ctrrfs)(layout, *uplo, *trans, *diag, *n, *nrhs,
                                     a_r, lda_r, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_ctrrfs_work)(layout, *uplo, *trans, *diag, *n,
                                          *nrhs, a_r, lda_r, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)x_r);
#endif
    *info = lapacke_test_info("CTRRFS", ret);
}

/******************************************************************************
 * CTRTRI( UPLO, DIAG, N, A, LDA, INFO )
 ******************************************************************************/
#define CTRTRI_TEST LAPACK_GLOBAL_SUFFIX(ctrtri_test, CTRTRI_TEST)
void CTRTRI_TEST(const char *uplo, const char *diag, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ctr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CTRTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ctrtri)(layout, *uplo, *diag, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_ctrtri_work)(layout, *uplo, *diag, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ctr_rm_to_cm(*uplo, *diag, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CTRTRI", ret);
}

/******************************************************************************
 * CTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define CTRTRS_TEST LAPACK_GLOBAL_SUFFIX(ctrtrs_test, CTRTRS_TEST)
void CTRTRS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ctr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CTRTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ctrtrs)(layout, *uplo, *trans, *diag, *n, *nrhs,
                                     a_r, lda_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ctrtrs_work)(layout, *uplo, *trans, *diag, *n,
                                          *nrhs, a_r, lda_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CTRTRS", ret);
}

// ========================================================================== //
//                        Triangular band matrices (TB)                       //
// ========================================================================== //

/******************************************************************************
 * CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, RWORK, INFO )
 ******************************************************************************/
#define CTBCON_TEST LAPACK_GLOBAL_SUFFIX(ctbcon_test, CTBCON_TEST)
void CTBCON_TEST(const char *norm, const char *uplo, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_complex_float *ab, const lapack_int *ldab,
                 float *rcond, lapack_complex_float *work, float *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_ctb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("CTBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ctbcon)(layout, *norm, *uplo, *diag, *n, *kd, ab_r,
                                     ldab_r, rcond);
#else
    ret = API_SUFFIX(LAPACKE_ctbcon_work)(layout, *norm, *uplo, *diag, *n, *kd,
                                          ab_r, ldab_r, rcond, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("CTBCON", ret);
}

/******************************************************************************
 * CTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, X, LDX, FERR, BERR,
 * WORK, RWORK, INFO )
 ******************************************************************************/
#define CTBRFS_TEST LAPACK_GLOBAL_SUFFIX(ctbrfs_test, CTBRFS_TEST)
void CTBRFS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const lapack_complex_float *ab,
                 const lapack_int *ldab, const lapack_complex_float *b,
                 const lapack_int *ldb, const lapack_complex_float *x,
                 const lapack_int *ldx, float *ferr, float *berr,
                 lapack_complex_float *work, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
    const lapack_complex_float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_ctb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)x_r);
        lapacke_test_report_alloc_failure("CTBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ctbrfs)(layout, *uplo, *trans, *diag, *n, *kd,
                                     *nrhs, ab_r, ldab_r, b_r, ldb_r, x_r,
                                     ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_ctbrfs_work)(layout, *uplo, *trans, *diag, *n, *kd,
                                          *nrhs, ab_r, ldab_r, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)x_r);
#endif
    *info = lapacke_test_info("CTBRFS", ret);
}

/******************************************************************************
 * CTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define CTBTRS_TEST LAPACK_GLOBAL_SUFFIX(ctbtrs_test, CTBTRS_TEST)
void CTBTRS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const lapack_complex_float *ab,
                 const lapack_int *ldab, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_ctb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CTBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ctbtrs)(layout, *uplo, *trans, *diag, *n, *kd,
                                     *nrhs, ab_r, ldab_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ctbtrs_work)(layout, *uplo, *trans, *diag, *n, *kd,
                                          *nrhs, ab_r, ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CTBTRS", ret);
}

// ========================================================================== //
//                   Orthogonal factorizations (QR/LQ/QL/RQ)                  //
// ========================================================================== //

/******************************************************************************
 * CGELQ2( M, N, A, LDA, TAU, WORK, INFO )
 ******************************************************************************/
#define CGELQ2_TEST LAPACK_GLOBAL_SUFFIX(cgelq2_test, CGELQ2_TEST)
void CGELQ2_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *tau, lapack_complex_float *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGELQ2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgelq2)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cgelq2_work)(layout, *m, *n, a_r, lda_r, tau,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGELQ2", ret);
}

/******************************************************************************
 * CGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGELQF_TEST LAPACK_GLOBAL_SUFFIX(cgelqf_test, CGELQF_TEST)
void CGELQF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgelqf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("CGELQF", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGELQF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgelqf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cgelqf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGELQF", ret);
}

/******************************************************************************
 * CGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGEQLF_TEST LAPACK_GLOBAL_SUFFIX(cgeqlf_test, CGEQLF_TEST)
void CGEQLF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgeqlf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("CGEQLF", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGEQLF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqlf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cgeqlf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGEQLF", ret);
}

/******************************************************************************
 * CGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGEQR_TEST LAPACK_GLOBAL_SUFFIX(cgeqr_test, CGEQR_TEST)
void CGEQR_TEST(const lapack_int *m, const lapack_int *n,
                lapack_complex_float *a, const lapack_int *lda,
                lapack_complex_float *t, const lapack_int *tsize,
                lapack_complex_float *work, const lapack_int *lwork,
                lapack_int *info)
{
    lapack_int ret = 0;
    if (*tsize == -1 || *lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgeqr_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                             t, *tsize, work, *lwork);
        *info = lapacke_test_info("CGEQR", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGEQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqr)(layout, *m, *n, a_r, lda_r, t, *tsize);
#else
    ret = API_SUFFIX(LAPACKE_cgeqr_work)(layout, *m, *n, a_r, lda_r, t, *tsize,
                                         work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGEQR", ret);
}

/******************************************************************************
 * CGEQR2( M, N, A, LDA, TAU, WORK, INFO )
 ******************************************************************************/
#define CGEQR2_TEST LAPACK_GLOBAL_SUFFIX(cgeqr2_test, CGEQR2_TEST)
void CGEQR2_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *tau, lapack_complex_float *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGEQR2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqr2)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cgeqr2_work)(layout, *m, *n, a_r, lda_r, tau,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGEQR2", ret);
}

/******************************************************************************
 * CGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGEQRF_TEST LAPACK_GLOBAL_SUFFIX(cgeqrf_test, CGEQRF_TEST)
void CGEQRF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgeqrf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("CGEQRF", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGEQRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqrf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cgeqrf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGEQRF", ret);
}

/******************************************************************************
 * CGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGEQRFP_TEST LAPACK_GLOBAL_SUFFIX(cgeqrfp_test, CGEQRFP_TEST)
void CGEQRFP_TEST(const lapack_int *m, const lapack_int *n,
                  lapack_complex_float *a, const lapack_int *lda,
                  lapack_complex_float *tau, lapack_complex_float *work,
                  const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgeqrfp_work)(LAPACK_COL_MAJOR, *m, *n, a,
                                               *lda, tau, work, *lwork);
        *info = lapacke_test_info("CGEQRFP", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGEQRFP", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqrfp)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cgeqrfp_work)(layout, *m, *n, a_r, lda_r, tau,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGEQRFP", ret);
}

/******************************************************************************
 * CGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGERQF_TEST LAPACK_GLOBAL_SUFFIX(cgerqf_test, CGERQF_TEST)
void CGERQF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgerqf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("CGERQF", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGERQF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgerqf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cgerqf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGERQF", ret);
}

/******************************************************************************
 * CUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CUNGLQ_TEST LAPACK_GLOBAL_SUFFIX(cunglq_test, CUNGLQ_TEST)
void CUNGLQ_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cunglq_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("CUNGLQ", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CUNGLQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cunglq)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cunglq_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CUNGLQ", ret);
}

/******************************************************************************
 * CUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CUNGQL_TEST LAPACK_GLOBAL_SUFFIX(cungql_test, CUNGQL_TEST)
void CUNGQL_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cungql_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("CUNGQL", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CUNGQL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cungql)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cungql_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CUNGQL", ret);
}

/******************************************************************************
 * CUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CUNGQR_TEST LAPACK_GLOBAL_SUFFIX(cungqr_test, CUNGQR_TEST)
void CUNGQR_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cungqr_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("CUNGQR", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CUNGQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cungqr)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cungqr_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CUNGQR", ret);
}

/******************************************************************************
 * CUNGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CUNGRQ_TEST LAPACK_GLOBAL_SUFFIX(cungrq_test, CUNGRQ_TEST)
void CUNGRQ_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cungrq_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("CUNGRQ", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CUNGRQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cungrq)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cungrq_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CUNGRQ", ret);
}

/******************************************************************************
 * CUNMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define CUNMLQ_TEST LAPACK_GLOBAL_SUFFIX(cunmlq_test, CUNMLQ_TEST)
void CUNMLQ_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *tau, lapack_complex_float *c,
                 const lapack_int *ldc, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cunmlq_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("CUNMLQ", ret);
        return;
    }

    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_cge_cm_to_rm(*k, nq, a, *lda, &lda_r);
    c_r = lapacke_test_cge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("CUNMLQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cunmlq)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_cunmlq_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("CUNMLQ", ret);
}

/******************************************************************************
 * CUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define CUNMQL_TEST LAPACK_GLOBAL_SUFFIX(cunmql_test, CUNMQL_TEST)
void CUNMQL_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *tau, lapack_complex_float *c,
                 const lapack_int *ldc, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cunmql_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("CUNMQL", ret);
        return;
    }

    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_cge_cm_to_rm(nq, *k, a, *lda, &lda_r);
    c_r = lapacke_test_cge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("CUNMQL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cunmql)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_cunmql_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("CUNMQL", ret);
}

/******************************************************************************
 * CUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define CUNMQR_TEST LAPACK_GLOBAL_SUFFIX(cunmqr_test, CUNMQR_TEST)
void CUNMQR_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *tau, lapack_complex_float *c,
                 const lapack_int *ldc, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cunmqr_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("CUNMQR", ret);
        return;
    }

    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_cge_cm_to_rm(nq, *k, a, *lda, &lda_r);
    c_r = lapacke_test_cge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("CUNMQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cunmqr)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_cunmqr_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("CUNMQR", ret);
}

// ========================================================================== //
//            Blocked and tall-skinny factorizations (QRT/TSQR/HR)            //
// ========================================================================== //

/******************************************************************************
 * CGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
 ******************************************************************************/
#define CGEQRT_TEST LAPACK_GLOBAL_SUFFIX(cgeqrt_test, CGEQRT_TEST)
void CGEQRT_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *nb,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *t, const lapack_int *ldt,
                 lapack_complex_float *work, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_cge_cm_to_rm(*nb, MIN(*m, *n), t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("CGEQRT", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqrt)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                     ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_cgeqrt_work)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                          ldt_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*nb, MIN(*m, *n), t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("CGEQRT", ret);
}

/******************************************************************************
 * CGEQRT2( M, N, A, LDA, T, LDT, INFO )
 ******************************************************************************/
#define CGEQRT2_TEST LAPACK_GLOBAL_SUFFIX(cgeqrt2_test, CGEQRT2_TEST)
void CGEQRT2_TEST(const lapack_int *m, const lapack_int *n,
                  lapack_complex_float *a, const lapack_int *lda,
                  lapack_complex_float *t, const lapack_int *ldt,
                  lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_cge_cm_to_rm(*n, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("CGEQRT2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqrt2)(layout, *m, *n, a_r, lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_cgeqrt2_work)(layout, *m, *n, a_r, lda_r, t_r,
                                           ldt_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("CGEQRT2", ret);
}

/******************************************************************************
 * CGEQRT3( M, N, A, LDA, T, LDT, INFO )
 ******************************************************************************/
#define CGEQRT3_TEST LAPACK_GLOBAL_SUFFIX(cgeqrt3_test, CGEQRT3_TEST)
void CGEQRT3_TEST(const lapack_int *m, const lapack_int *n,
                  lapack_complex_float *a, const lapack_int *lda,
                  lapack_complex_float *t, const lapack_int *ldt,
                  lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_cge_cm_to_rm(*n, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("CGEQRT3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqrt3)(layout, *m, *n, a_r, lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_cgeqrt3_work)(layout, *m, *n, a_r, lda_r, t_r,
                                           ldt_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*n, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("CGEQRT3", ret);
}

/******************************************************************************
 * CGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGETSQRHRT_TEST LAPACK_GLOBAL_SUFFIX(cgetsqrhrt_test, CGETSQRHRT_TEST)
void CGETSQRHRT_TEST(const lapack_int *m, const lapack_int *n,
                     const lapack_int *mb1, const lapack_int *nb1,
                     const lapack_int *nb2, lapack_complex_float *a,
                     const lapack_int *lda, lapack_complex_float *t,
                     const lapack_int *ldt, lapack_complex_float *work,
                     const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgetsqrhrt_work)(LAPACK_COL_MAJOR, *m, *n,
                                                  *mb1, *nb1, *nb2, a, *lda, t,
                                                  *ldt, work, *lwork);
        *info = lapacke_test_info("CGETSQRHRT", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_cge_cm_to_rm(*nb2, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("CGETSQRHRT", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgetsqrhrt)(layout, *m, *n, *mb1, *nb1, *nb2, a_r,
                                         lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_cgetsqrhrt_work)(
        layout, *m, *n, *mb1, *nb1, *nb2, a_r, lda_r, t_r, ldt_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*nb2, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("CGETSQRHRT", ret);
}

/******************************************************************************
 * CUNHR_COL( M, N, NB, A, LDA, T, LDT, D, INFO )
 ******************************************************************************/
#define CUNHR_COL_TEST LAPACK_GLOBAL_SUFFIX(cunhr_col_test, CUNHR_COL_TEST)
void CUNHR_COL_TEST(const lapack_int *m, const lapack_int *n,
                    const lapack_int *nb, lapack_complex_float *a,
                    const lapack_int *lda, lapack_complex_float *t,
                    const lapack_int *ldt, lapack_complex_float *d,
                    lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_cge_cm_to_rm(*ldt, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("CUNHR_COL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cunhr_col)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                        ldt_r, d);
#else
    ret = API_SUFFIX(LAPACKE_cunhr_col_work)(layout, *m, *n, *nb, a_r, lda_r,
                                             t_r, ldt_r, d);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(*ldt, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("CUNHR_COL", ret);
}

// ========================================================================== //
//                    Rank-revealing factorizations (QP/TZ)                   //
// ========================================================================== //

/******************************************************************************
 * CGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, INFO )
 ******************************************************************************/
#define CGEQP3_TEST LAPACK_GLOBAL_SUFFIX(cgeqp3_test, CGEQP3_TEST)
void CGEQP3_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_int *jpvt, lapack_complex_float *tau,
                 lapack_complex_float *work, const lapack_int *lwork,
                 float *rwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgeqp3_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              jpvt, tau, work, *lwork, rwork);
        *info = lapacke_test_info("CGEQP3", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGEQP3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqp3)(layout, *m, *n, a_r, lda_r, jpvt, tau);
#else
    ret = API_SUFFIX(LAPACKE_cgeqp3_work)(layout, *m, *n, a_r, lda_r, jpvt, tau,
                                          work, *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGEQP3", ret);
}

/******************************************************************************
 * CTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define CTZRZF_TEST LAPACK_GLOBAL_SUFFIX(ctzrzf_test, CTZRZF_TEST)
void CTZRZF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ctzrzf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("CTZRZF", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CTZRZF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ctzrzf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_ctzrzf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CTZRZF", ret);
}

// ========================================================================== //
//                             Least squares (LS)                             //
// ========================================================================== //

/******************************************************************************
 * CGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGELS_TEST LAPACK_GLOBAL_SUFFIX(cgels_test, CGELS_TEST)
void CGELS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                const lapack_int *nrhs, lapack_complex_float *a,
                const lapack_int *lda, lapack_complex_float *b,
                const lapack_int *ldb, lapack_complex_float *work,
                const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgels_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                             *nrhs, a, *lda, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("CGELS", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGELS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgels)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgels_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                         lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGELS", ret);
}

/******************************************************************************
 * CGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, RWORK,
 * IWORK, INFO )
 ******************************************************************************/
#define CGELSD_TEST LAPACK_GLOBAL_SUFFIX(cgelsd_test, CGELSD_TEST)
void CGELSD_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_float *a,
                 const lapack_int *lda, lapack_complex_float *b,
                 const lapack_int *ldb, float *s, const float *rcond,
                 lapack_int *rank, lapack_complex_float *work,
                 const lapack_int *lwork, float *rwork, lapack_int *iwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgelsd_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, s, *rcond, rank,
                                              work, *lwork, rwork, iwork);
        *info = lapacke_test_info("CGELSD", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGELSD", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgelsd)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, s, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_cgelsd_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, s, *rcond, rank, work,
                                          *lwork, rwork, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGELSD", ret);
}

/******************************************************************************
 * CGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, RWORK, INFO
 * )
 ******************************************************************************/
#define CGELSS_TEST LAPACK_GLOBAL_SUFFIX(cgelss_test, CGELSS_TEST)
void CGELSS_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_float *a,
                 const lapack_int *lda, lapack_complex_float *b,
                 const lapack_int *ldb, float *s, const float *rcond,
                 lapack_int *rank, lapack_complex_float *work,
                 const lapack_int *lwork, float *rwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgelss_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, s, *rcond, rank,
                                              work, *lwork, rwork);
        *info = lapacke_test_info("CGELSS", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGELSS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgelss)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, s, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_cgelss_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, s, *rcond, rank, work,
                                          *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGELSS", ret);
}

/******************************************************************************
 * CGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, RWORK,
 * INFO )
 ******************************************************************************/
#define CGELSY_TEST LAPACK_GLOBAL_SUFFIX(cgelsy_test, CGELSY_TEST)
void CGELSY_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, lapack_complex_float *a,
                 const lapack_int *lda, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_int *jpvt, const float *rcond,
                 lapack_int *rank, lapack_complex_float *work,
                 const lapack_int *lwork, float *rwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgelsy_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, jpvt, *rcond,
                                              rank, work, *lwork, rwork);
        *info = lapacke_test_info("CGELSY", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGELSY", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgelsy)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, jpvt, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_cgelsy_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, jpvt, *rcond, rank, work,
                                          *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGELSY", ret);
}

/******************************************************************************
 * CGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define CGETSLS_TEST LAPACK_GLOBAL_SUFFIX(cgetsls_test, CGETSLS_TEST)
void CGETSLS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                  const lapack_int *nrhs, lapack_complex_float *a,
                  const lapack_int *lda, lapack_complex_float *b,
                  const lapack_int *ldb, lapack_complex_float *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgetsls_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                               *nrhs, a, *lda, b, *ldb, work,
                                               *lwork);
        *info = lapacke_test_info("CGETSLS", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGETSLS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgetsls)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                      b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgetsls_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                           lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGETSLS", ret);
}

// ========================================================================== //
//                             SVD and bidiagonal                             //
// ========================================================================== //

/******************************************************************************
 * CBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, RWORK, INFO
 * )
 ******************************************************************************/
#define CBDSQR_TEST LAPACK_GLOBAL_SUFFIX(cbdsqr_test, CBDSQR_TEST)
void CBDSQR_TEST(const char *uplo, const lapack_int *n, const lapack_int *ncvt,
                 const lapack_int *nru, const lapack_int *ncc, float *d,
                 float *e, lapack_complex_float *vt, const lapack_int *ldvt,
                 lapack_complex_float *u, const lapack_int *ldu,
                 lapack_complex_float *c, const lapack_int *ldc, float *rwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *vt_r = vt;
    lapack_int ldvt_r = *ldvt;
    lapack_complex_float *u_r = u;
    lapack_int ldu_r = *ldu;
    lapack_complex_float *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    vt_r = lapacke_test_cge_cm_to_rm(*n, *ncvt, vt, *ldvt, &ldvt_r);
    u_r = lapacke_test_cge_cm_to_rm(*nru, *n, u, *ldu, &ldu_r);
    c_r = lapacke_test_cge_cm_to_rm(*n, *ncc, c, *ldc, &ldc_r);
    if (vt_r == NULL || u_r == NULL || c_r == NULL) {
        LAPACKE_free(vt_r);
        LAPACKE_free(u_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("CBDSQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cbdsqr)(layout, *uplo, *n, *ncvt, *nru, *ncc, d, e,
                                     vt_r, ldvt_r, u_r, ldu_r, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_cbdsqr_work)(layout, *uplo, *n, *ncvt, *nru, *ncc,
                                          d, e, vt_r, ldvt_r, u_r, ldu_r, c_r,
                                          ldc_r, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *ncvt, vt_r, ldvt_r, vt, *ldvt);
    lapacke_test_cge_rm_to_cm(*nru, *n, u_r, ldu_r, u, *ldu);
    lapacke_test_cge_rm_to_cm(*n, *ncc, c_r, ldc_r, c, *ldc);
    LAPACKE_free(vt_r);
    LAPACKE_free(u_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("CBDSQR", ret);
}

/******************************************************************************
 * CGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK,
 * INFO )
 ******************************************************************************/
#define CGESVD_TEST LAPACK_GLOBAL_SUFFIX(cgesvd_test, CGESVD_TEST)
void CGESVD_TEST(const char *jobu, const char *jobvt, const lapack_int *m,
                 const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, float *s, lapack_complex_float *u,
                 const lapack_int *ldu, lapack_complex_float *vt,
                 const lapack_int *ldvt, lapack_complex_float *work,
                 const lapack_int *lwork, float *rwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN jobu_len, FORTRAN_STRLEN jobvt_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_cgesvd_work)(LAPACK_COL_MAJOR, *jobu, *jobvt,
                                              *m, *n, a, *lda, s, u, *ldu, vt,
                                              *ldvt, work, *lwork, rwork);
        *info = lapacke_test_info("CGESVD", ret);
        return;
    }

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
    lapack_complex_float *u_r = u;
    lapack_int ldu_r = *ldu;
    lapack_complex_float *vt_r = vt;
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
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    u_r = lapacke_test_cge_cm_to_rm(urows, ucols, u, *ldu, &ldu_r);
    vt_r = lapacke_test_cge_cm_to_rm(vtrows, *n, vt, *ldvt, &ldvt_r);
    if (a_r == NULL || u_r == NULL || vt_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(u_r);
        LAPACKE_free(vt_r);
        lapacke_test_report_alloc_failure("CGESVD", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgesvd)(layout, *jobu, *jobvt, *m, *n, a_r, lda_r,
                                     s, u_r, ldu_r, vt_r, ldvt_r, rwork);
#else
    ret = API_SUFFIX(LAPACKE_cgesvd_work)(layout, *jobu, *jobvt, *m, *n, a_r,
                                          lda_r, s, u_r, ldu_r, vt_r, ldvt_r,
                                          work, *lwork, rwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(urows, ucols, u_r, ldu_r, u, *ldu);
    lapacke_test_cge_rm_to_cm(vtrows, *n, vt_r, ldvt_r, vt, *ldvt);
    LAPACKE_free(a_r);
    LAPACKE_free(u_r);
    LAPACKE_free(vt_r);
#endif
    *info = lapacke_test_info("CGESVD", ret);
}

// ========================================================================== //
//                             Auxiliary routines                             //
// ========================================================================== //

/******************************************************************************
 * CLACGV( N, X, INCX )
 ******************************************************************************/
#define CLACGV_TEST LAPACK_GLOBAL_SUFFIX(clacgv_test, CLACGV_TEST)
void CLACGV_TEST(const lapack_int *n, lapack_complex_float *x,
                 const lapack_int *incx)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_clacgv)(*n, x, *incx);
#else
    ret = API_SUFFIX(LAPACKE_clacgv_work)(*n, x, *incx);
#endif

    lapacke_test_info_unshifted("CLACGV", ret);
}

/******************************************************************************
 * REAL FUNCTION CLANGE( NORM, M, N, A, LDA, WORK )
 ******************************************************************************/
#define CLANGE_TEST LAPACK_GLOBAL_SUFFIX(clange_test, CLANGE_TEST)
lapack_float_return CLANGE_TEST(const char *norm, const lapack_int *m,
                                const lapack_int *n,
                                const lapack_complex_float *a,
                                const lapack_int *lda, float *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                                ,
                                FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_float_return res = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("CLANGE", &info);
        return (lapack_float_return)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_clange)(layout, *norm, *m, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_clange_work)(layout, *norm, *m, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * REAL FUNCTION CLANHE( NORM, UPLO, N, A, LDA, WORK )
 ******************************************************************************/
#define CLANHE_TEST LAPACK_GLOBAL_SUFFIX(clanhe_test, CLANHE_TEST)
lapack_float_return CLANHE_TEST(const char *norm, const char *uplo,
                                const lapack_int *n,
                                const lapack_complex_float *a,
                                const lapack_int *lda, float *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                                ,
                                FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_float_return res = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_che_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("CLANHE", &info);
        return (lapack_float_return)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_clanhe)(layout, *norm, *uplo, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_clanhe_work)(layout, *norm, *uplo, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * REAL FUNCTION CLANSY( NORM, UPLO, N, A, LDA, WORK )
 ******************************************************************************/
#define CLANSY_TEST LAPACK_GLOBAL_SUFFIX(clansy_test, CLANSY_TEST)
lapack_float_return CLANSY_TEST(const char *norm, const char *uplo,
                                const lapack_int *n,
                                const lapack_complex_float *a,
                                const lapack_int *lda, float *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                                ,
                                FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_float_return res = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_csy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("CLANSY", &info);
        return (lapack_float_return)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_clansy)(layout, *norm, *uplo, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_clansy_work)(layout, *norm, *uplo, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * REAL FUNCTION CLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
 ******************************************************************************/
#define CLANTR_TEST LAPACK_GLOBAL_SUFFIX(clantr_test, CLANTR_TEST)
lapack_float_return CLANTR_TEST(const char *norm, const char *uplo,
                                const char *diag, const lapack_int *m,
                                const lapack_int *n,
                                const lapack_complex_float *a,
                                const lapack_int *lda, float *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                                ,
                                FORTRAN_STRLEN norm_len,
                                FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_float_return res = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("CLANTR", &info);
        return (lapack_float_return)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_clantr)(layout, *norm, *uplo, *diag, *m, *n, a_r,
                                     lda_r);
#else
    res = API_SUFFIX(LAPACKE_clantr_work)(layout, *norm, *uplo, *diag, *m, *n,
                                          a_r, lda_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * CLARNV( IDIST, ISEED, N, X )
 ******************************************************************************/
#define CLARNV_TEST LAPACK_GLOBAL_SUFFIX(clarnv_test, CLARNV_TEST)
void CLARNV_TEST(const lapack_int *idist, lapack_int *iseed,
                 const lapack_int *n, lapack_complex_float *x)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_clarnv)(*idist, iseed, *n, x);
#else
    ret = API_SUFFIX(LAPACKE_clarnv_work)(*idist, iseed, *n, x);
#endif

    lapacke_test_info_unshifted("CLARNV", ret);
}

/******************************************************************************
 * CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
 ******************************************************************************/
#define CLASCL_TEST LAPACK_GLOBAL_SUFFIX(clascl_test, CLASCL_TEST)
void CLASCL_TEST(const char *type, const lapack_int *kl, const lapack_int *ku,
                 const float *cfrom, const float *cto, const lapack_int *m,
                 const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN type_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int arows = API_SUFFIX(LAPACKE_lsame)(*type, 'b')   ? *kl + 1
                             : API_SUFFIX(LAPACKE_lsame)(*type, 'q') ? *ku + 1
                             : API_SUFFIX(LAPACKE_lsame)(*type, 'z')
                                 ? 2 * *kl + *ku + 1
                                 : *m;
    a_r = lapacke_test_cge_cm_to_rm(arows, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CLASCL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_clascl)(layout, *type, *kl, *ku, *cfrom, *cto, *m,
                                     *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_clascl_work)(layout, *type, *kl, *ku, *cfrom, *cto,
                                          *m, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(arows, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CLASCL", ret);
}

/******************************************************************************
 * CLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
 ******************************************************************************/
#define CLASET_TEST LAPACK_GLOBAL_SUFFIX(claset_test, CLASET_TEST)
void CLASET_TEST(const char *uplo, const lapack_int *m, const lapack_int *n,
                 const lapack_complex_float *alpha,
                 const lapack_complex_float *beta, lapack_complex_float *a,
                 const lapack_int *lda
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("CLASET", &info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_claset)(layout, *uplo, *m, *n, *alpha, *beta, a_r,
                                     lda_r);
#else
    ret = API_SUFFIX(LAPACKE_claset_work)(layout, *uplo, *m, *n, *alpha, *beta,
                                          a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    lapacke_test_info("CLASET", ret);
}

/******************************************************************************
 * CLASWP( N, A, LDA, K1, K2, IPIV, INCX )
 ******************************************************************************/
#define CLASWP_TEST LAPACK_GLOBAL_SUFFIX(claswp_test, CLASWP_TEST)
void CLASWP_TEST(const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, const lapack_int *k1,
                 const lapack_int *k2, const lapack_int *ipiv,
                 const lapack_int *incx)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(
        lapacke_test_laswp_rows(*k1, *k2, ipiv, *incx), *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("CLASWP", &info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_claswp)(layout, *n, a_r, lda_r, *k1, *k2, ipiv,
                                     *incx);
#else
    ret = API_SUFFIX(LAPACKE_claswp_work)(layout, *n, a_r, lda_r, *k1, *k2,
                                          ipiv, *incx);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(lapacke_test_laswp_rows(*k1, *k2, ipiv, *incx),
                              *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    lapacke_test_info("CLASWP", ret);
}

/******************************************************************************
 * CROT( N, CX, INCX, CY, INCY, C, S )
 ******************************************************************************/
#define CROT_TEST LAPACK_GLOBAL_SUFFIX(crot_test, CROT_TEST)
void CROT_TEST(const lapack_int *n, lapack_complex_float *cx,
               const lapack_int *incx, lapack_complex_float *cy,
               const lapack_int *incy, const float *c,
               const lapack_complex_float *s)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_crot)(*n, cx, *incx, cy, *incy, *c, *s);
#else
    ret = API_SUFFIX(LAPACKE_crot_work)(*n, cx, *incx, cy, *incy, *c, *s);
#endif

    lapacke_test_info_unshifted("CROT", ret);
}
