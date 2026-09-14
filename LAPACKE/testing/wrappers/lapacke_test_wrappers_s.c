/******************************************************************************
 * LAPACKE test wrappers (single precision real)
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
 * SGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define SGECON_TEST LAPACK_GLOBAL_SUFFIX(sgecon_test, SGECON_TEST)
void SGECON_TEST(const char *norm, const lapack_int *n, const float *a,
                 const lapack_int *lda, const float *anorm, float *rcond,
                 float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGECON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgecon)(layout, *norm, *n, a_r, lda_r, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_sgecon_work)(layout, *norm, *n, a_r, lda_r, *anorm,
                                          rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("SGECON", ret);
}

/******************************************************************************
 * SGECXX( FACT, USESD, M, N, SESEL_ROWS, SEL_DESEL_COLS, KMAXFREE, ABSTOL,
 * RELTOL, A, LDA, K, MAXC2NRMK, RELMAXC2NRMK, FNRMK, IPIV, JPIV, TAU, C, LDC,
 * QRC, LDQRC, X, LDX, WORK, LWORK, IWORK, LIWORK, INFO )
 ******************************************************************************/
#define SGECXX_TEST LAPACK_GLOBAL_SUFFIX(sgecxx_test, SGECXX_TEST)
void SGECXX_TEST(const char *fact, const char *usesd, const lapack_int *m,
                 const lapack_int *n, const lapack_int *sesel_rows,
                 const lapack_int *sel_desel_cols, const lapack_int *kmaxfree,
                 const float *abstol, const float *reltol, float *a,
                 const lapack_int *lda, lapack_int *k, float *maxc2nrmk,
                 float *relmaxc2nrmk, float *fnrmk, lapack_int *ipiv,
                 lapack_int *jpiv, float *tau, float *c, const lapack_int *ldc,
                 float *qrc, const lapack_int *ldqrc, float *x,
                 const lapack_int *ldx, float *work, const lapack_int *lwork,
                 lapack_int *iwork, const lapack_int *liwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN usesd_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgecxx_work)(
            LAPACK_COL_MAJOR, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
            (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a, *lda,
            k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c, *ldc, qrc,
            *ldqrc, x, *ldx, work, *lwork, iwork, *liwork);
        *info = lapacke_test_info("SGECXX", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *c_r = c;
    lapack_int ldc_r = *ldc;
    float *qrc_r = qrc;
    lapack_int ldqrc_r = *ldqrc;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    c_r = lapacke_test_sge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    qrc_r = lapacke_test_sge_cm_to_rm(*m, MIN(*m, *n), qrc, *ldqrc, &ldqrc_r);
    x_r = lapacke_test_sge_cm_to_rm(*m, *n, x, *ldx, &ldx_r);
    if (a_r == NULL || c_r == NULL || qrc_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(c_r);
        LAPACKE_free(qrc_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SGECXX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgecxx)(
        layout, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
        (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a_r, lda_r,
        k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c_r, ldc_r, qrc_r,
        ldqrc_r, x_r, ldx_r);
#else
    ret = API_SUFFIX(LAPACKE_sgecxx_work)(
        layout, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
        (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a_r, lda_r,
        k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c_r, ldc_r, qrc_r,
        ldqrc_r, x_r, ldx_r, work, *lwork, iwork, *liwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    lapacke_test_sge_rm_to_cm(*m, MIN(*m, *n), qrc_r, ldqrc_r, qrc, *ldqrc);
    lapacke_test_sge_rm_to_cm(*m, *n, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(c_r);
    LAPACKE_free(qrc_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SGECXX", ret);
}

/******************************************************************************
 * SGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )
 ******************************************************************************/
#define SGEEQU_TEST LAPACK_GLOBAL_SUFFIX(sgeequ_test, SGEEQU_TEST)
void SGEEQU_TEST(const lapack_int *m, const lapack_int *n, const float *a,
                 const lapack_int *lda, float *r, float *c, float *rowcnd,
                 float *colcnd, float *amax, lapack_int *info)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGEEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeequ)(layout, *m, *n, a_r, lda_r, r, c, rowcnd,
                                     colcnd, amax);
#else
    ret = API_SUFFIX(LAPACKE_sgeequ_work)(layout, *m, *n, a_r, lda_r, r, c,
                                          rowcnd, colcnd, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("SGEEQU", ret);
}

/******************************************************************************
 * SGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, IWORK, INFO )
 ******************************************************************************/
#define SGERFS_TEST LAPACK_GLOBAL_SUFFIX(sgerfs_test, SGERFS_TEST)
void SGERFS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const float *a, const lapack_int *lda, const float *af,
                 const lapack_int *ldaf, const lapack_int *ipiv, const float *b,
                 const lapack_int *ldb, float *x, const lapack_int *ldx,
                 float *ferr, float *berr, float *work, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    const float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    af_r = lapacke_test_sge_cm_to_rm(*n, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SGERFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgerfs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     af_r, ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_sgerfs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SGERFS", ret);
}

/******************************************************************************
 * SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SGESV_TEST LAPACK_GLOBAL_SUFFIX(sgesv_test, SGESV_TEST)
void SGESV_TEST(const lapack_int *n, const lapack_int *nrhs, float *a,
                const lapack_int *lda, lapack_int *ipiv, float *b,
                const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGESV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgesv)(layout, *n, *nrhs, a_r, lda_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgesv_work)(layout, *n, *nrhs, a_r, lda_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGESV", ret);
}

/******************************************************************************
 * SGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define SGESVX_TEST LAPACK_GLOBAL_SUFFIX(sgesvx_test, SGESVX_TEST)
void SGESVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *nrhs, float *a, const lapack_int *lda,
                 float *af, const lapack_int *ldaf, lapack_int *ipiv,
                 char *equed, float *r, float *c, float *b,
                 const lapack_int *ldb, float *x, const lapack_int *ldx,
                 float *rcond, float *ferr, float *berr, float *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
    float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    af_r = lapacke_test_sge_cm_to_rm(*n, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(af_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SGESVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgesvx)(
        layout, *fact, *trans, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, equed,
        r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work);
#else
    ret = API_SUFFIX(LAPACKE_sgesvx_work)(
        layout, *fact, *trans, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, equed,
        r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*n, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(af_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SGESVX", ret);
}

/******************************************************************************
 * SGETF2( M, N, A, LDA, IPIV, INFO )
 ******************************************************************************/
#define SGETF2_TEST LAPACK_GLOBAL_SUFFIX(sgetf2_test, SGETF2_TEST)
void SGETF2_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGETF2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgetf2)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_sgetf2_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGETF2", ret);
}

/******************************************************************************
 * SGETRF( M, N, A, LDA, IPIV, INFO )
 ******************************************************************************/
#define SGETRF_TEST LAPACK_GLOBAL_SUFFIX(sgetrf_test, SGETRF_TEST)
void SGETRF_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgetrf)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_sgetrf_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGETRF", ret);
}

/******************************************************************************
 * SGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGETRI_TEST LAPACK_GLOBAL_SUFFIX(sgetri_test, SGETRI_TEST)
void SGETRI_TEST(const lapack_int *n, float *a, const lapack_int *lda,
                 const lapack_int *ipiv, float *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgetri_work)(LAPACK_COL_MAJOR, *n, a, *lda,
                                              ipiv, work, *lwork);
        *info = lapacke_test_info("SGETRI", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgetri)(layout, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_sgetri_work)(layout, *n, a_r, lda_r, ipiv, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGETRI", ret);
}

/******************************************************************************
 * SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SGETRS_TEST LAPACK_GLOBAL_SUFFIX(sgetrs_test, SGETRS_TEST)
void SGETRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const float *a, const lapack_int *lda, const lapack_int *ipiv,
                 float *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgetrs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgetrs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGETRS", ret);
}

// ========================================================================== //
//                         General band matrices (GB)                         //
// ========================================================================== //

/******************************************************************************
 * SGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define SGBCON_TEST LAPACK_GLOBAL_SUFFIX(sgbcon_test, SGBCON_TEST)
void SGBCON_TEST(const char *norm, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const float *ab, const lapack_int *ldab,
                 const lapack_int *ipiv, const float *anorm, float *rcond,
                 float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_sgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("SGBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgbcon)(layout, *norm, *n, *kl, *ku, ab_r, ldab_r,
                                     ipiv, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_sgbcon_work)(layout, *norm, *n, *kl, *ku, ab_r,
                                          ldab_r, ipiv, *anorm, rcond, work,
                                          iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("SGBCON", ret);
}

/******************************************************************************
 * SGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFO )
 ******************************************************************************/
#define SGBEQU_TEST LAPACK_GLOBAL_SUFFIX(sgbequ_test, SGBEQU_TEST)
void SGBEQU_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const float *ab, const lapack_int *ldab,
                 float *r, float *c, float *rowcnd, float *colcnd, float *amax,
                 lapack_int *info)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_sgb_cm_to_rm(*m, *n, *kl, *ku, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("SGBEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgbequ)(layout, *m, *n, *kl, *ku, ab_r, ldab_r, r,
                                     c, rowcnd, colcnd, amax);
#else
    ret = API_SUFFIX(LAPACKE_sgbequ_work)(layout, *m, *n, *kl, *ku, ab_r,
                                          ldab_r, r, c, rowcnd, colcnd, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("SGBEQU", ret);
}

/******************************************************************************
 * SGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX,
 * FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define SGBRFS_TEST LAPACK_GLOBAL_SUFFIX(sgbrfs_test, SGBRFS_TEST)
void SGBRFS_TEST(const char *trans, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_int *nrhs, const float *ab,
                 const lapack_int *ldab, const float *afb,
                 const lapack_int *ldafb, const lapack_int *ipiv,
                 const float *b, const lapack_int *ldb, float *x,
                 const lapack_int *ldx, float *ferr, float *berr, float *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const float *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_sgb_cm_to_rm(*n, *n, *kl, *ku, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_sgb_cm_to_rm(*n, *n, *kl, *kl + *ku, afb, *ldafb,
                                      &ldafb_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)afb_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SGBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgbrfs)(layout, *trans, *n, *kl, *ku, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, ipiv, b_r, ldb_r,
                                     x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_sgbrfs_work)(
        layout, *trans, *n, *kl, *ku, *nrhs, ab_r, ldab_r, afb_r, ldafb_r, ipiv,
        b_r, ldb_r, x_r, ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)afb_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SGBRFS", ret);
}

/******************************************************************************
 * SGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SGBSV_TEST LAPACK_GLOBAL_SUFFIX(sgbsv_test, SGBSV_TEST)
void SGBSV_TEST(const lapack_int *n, const lapack_int *kl, const lapack_int *ku,
                const lapack_int *nrhs, float *ab, const lapack_int *ldab,
                lapack_int *ipiv, float *b, const lapack_int *ldb,
                lapack_int *info)
{
    lapack_int ret = 0;
    float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_sgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGBSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgbsv)(layout, *n, *kl, *ku, *nrhs, ab_r, ldab_r,
                                    ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgbsv_work)(layout, *n, *kl, *ku, *nrhs, ab_r,
                                         ldab_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sgb_rm_to_cm(*n, *n, *kl, *kl + *ku, ab_r, ldab_r, ab, *ldab);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGBSV", ret);
}

/******************************************************************************
 * SGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, EQUED, R,
 * C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define SGBSVX_TEST LAPACK_GLOBAL_SUFFIX(sgbsvx_test, SGBSVX_TEST)
void SGBSVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *kl, const lapack_int *ku,
                 const lapack_int *nrhs, float *ab, const lapack_int *ldab,
                 float *afb, const lapack_int *ldafb, lapack_int *ipiv,
                 char *equed, float *r, float *c, float *b,
                 const lapack_int *ldb, float *x, const lapack_int *ldx,
                 float *rcond, float *ferr, float *berr, float *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    float *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_sgb_cm_to_rm(*n, *n, *kl, *ku, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_sgb_cm_to_rm(*n, *n, *kl, *kl + *ku, afb, *ldafb,
                                      &ldafb_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(afb_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SGBSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgbsvx)(layout, *fact, *trans, *n, *kl, *ku, *nrhs,
                                     ab_r, ldab_r, afb_r, ldafb_r, ipiv, equed,
                                     r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr,
                                     berr, work);
#else
    ret = API_SUFFIX(LAPACKE_sgbsvx_work)(
        layout, *fact, *trans, *n, *kl, *ku, *nrhs, ab_r, ldab_r, afb_r,
        ldafb_r, ipiv, equed, r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr,
        work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sgb_rm_to_cm(*n, *n, *kl, *ku, ab_r, ldab_r, ab, *ldab);
    lapacke_test_sgb_rm_to_cm(*n, *n, *kl, *kl + *ku, afb_r, ldafb_r, afb,
                              *ldafb);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(ab_r);
    LAPACKE_free(afb_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SGBSVX", ret);
}

/******************************************************************************
 * SGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
 ******************************************************************************/
#define SGBTRF_TEST LAPACK_GLOBAL_SUFFIX(sgbtrf_test, SGBTRF_TEST)
void SGBTRF_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, float *ab, const lapack_int *ldab,
                 lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_sgb_cm_to_rm(*m, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("SGBTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgbtrf)(layout, *m, *n, *kl, *ku, ab_r, ldab_r,
                                     ipiv);
#else
    ret = API_SUFFIX(LAPACKE_sgbtrf_work)(layout, *m, *n, *kl, *ku, ab_r,
                                          ldab_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sgb_rm_to_cm(*m, *n, *kl, *kl + *ku, ab_r, ldab_r, ab, *ldab);
    LAPACKE_free(ab_r);
#endif
    *info = lapacke_test_info("SGBTRF", ret);
}

/******************************************************************************
 * SGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SGBTRS_TEST LAPACK_GLOBAL_SUFFIX(sgbtrs_test, SGBTRS_TEST)
void SGBTRS_TEST(const char *trans, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_int *nrhs, const float *ab,
                 const lapack_int *ldab, const lapack_int *ipiv, float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_sgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgbtrs)(layout, *trans, *n, *kl, *ku, *nrhs, ab_r,
                                     ldab_r, ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgbtrs_work)(layout, *trans, *n, *kl, *ku, *nrhs,
                                          ab_r, ldab_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGBTRS", ret);
}

// ========================================================================== //
//                      General tridiagonal matrices (GT)                     //
// ========================================================================== //

/******************************************************************************
 * SGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define SGTCON_TEST LAPACK_GLOBAL_SUFFIX(sgtcon_test, SGTCON_TEST)
void SGTCON_TEST(const char *norm, const lapack_int *n, const float *dl,
                 const float *d, const float *du, const float *du2,
                 const lapack_int *ipiv, const float *anorm, float *rcond,
                 float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgtcon)(*norm, *n, dl, d, du, du2, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_sgtcon_work)(*norm, *n, dl, d, du, du2, ipiv,
                                          *anorm, rcond, work, iwork);
#endif

    *info = lapacke_test_info_unshifted("SGTCON", ret);
}

/******************************************************************************
 * SGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX,
 * FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define SGTRFS_TEST LAPACK_GLOBAL_SUFFIX(sgtrfs_test, SGTRFS_TEST)
void SGTRFS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const float *dl, const float *d, const float *du,
                 const float *dlf, const float *df, const float *duf,
                 const float *du2, const lapack_int *ipiv, const float *b,
                 const lapack_int *ldb, float *x, const lapack_int *ldx,
                 float *ferr, float *berr, float *work, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SGTRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgtrfs)(layout, *trans, *n, *nrhs, dl, d, du, dlf,
                                     df, duf, du2, ipiv, b_r, ldb_r, x_r, ldx_r,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_sgtrfs_work)(layout, *trans, *n, *nrhs, dl, d, du,
                                          dlf, df, duf, du2, ipiv, b_r, ldb_r,
                                          x_r, ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SGTRFS", ret);
}

/******************************************************************************
 * SGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
 ******************************************************************************/
#define SGTSV_TEST LAPACK_GLOBAL_SUFFIX(sgtsv_test, SGTSV_TEST)
void SGTSV_TEST(const lapack_int *n, const lapack_int *nrhs, float *dl,
                float *d, float *du, float *b, const lapack_int *ldb,
                lapack_int *info)
{
    lapack_int ret = 0;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("SGTSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgtsv)(layout, *n, *nrhs, dl, d, du, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgtsv_work)(layout, *n, *nrhs, dl, d, du, b_r,
                                         ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGTSV", ret);
}

/******************************************************************************
 * SGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define SGTSVX_TEST LAPACK_GLOBAL_SUFFIX(sgtsvx_test, SGTSVX_TEST)
void SGTSVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *nrhs, const float *dl, const float *d,
                 const float *du, float *dlf, float *df, float *duf, float *du2,
                 lapack_int *ipiv, const float *b, const lapack_int *ldb,
                 float *x, const lapack_int *ldx, float *rcond, float *ferr,
                 float *berr, float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SGTSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgtsvx)(layout, *fact, *trans, *n, *nrhs, dl, d,
                                     du, dlf, df, duf, du2, ipiv, b_r, ldb_r,
                                     x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_sgtsvx_work)(
        layout, *fact, *trans, *n, *nrhs, dl, d, du, dlf, df, duf, du2, ipiv,
        b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SGTSVX", ret);
}

/******************************************************************************
 * SGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
 ******************************************************************************/
#define SGTTRF_TEST LAPACK_GLOBAL_SUFFIX(sgttrf_test, SGTTRF_TEST)
void SGTTRF_TEST(const lapack_int *n, float *dl, float *d, float *du,
                 float *du2, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgttrf)(*n, dl, d, du, du2, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_sgttrf_work)(*n, dl, d, du, du2, ipiv);
#endif

    *info = lapacke_test_info_unshifted("SGTTRF", ret);
}

/******************************************************************************
 * SGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SGTTRS_TEST LAPACK_GLOBAL_SUFFIX(sgttrs_test, SGTTRS_TEST)
void SGTTRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const float *dl, const float *d, const float *du,
                 const float *du2, const lapack_int *ipiv, float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("SGTTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgttrs)(layout, *trans, *n, *nrhs, dl, d, du, du2,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgttrs_work)(layout, *trans, *n, *nrhs, dl, d, du,
                                          du2, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGTTRS", ret);
}

// ========================================================================== //
//             Symmetric/Hermitian positive definite matrices (PO)            //
// ========================================================================== //

/******************************************************************************
 * SPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define SPOCON_TEST LAPACK_GLOBAL_SUFFIX(spocon_test, SPOCON_TEST)
void SPOCON_TEST(const char *uplo, const lapack_int *n, const float *a,
                 const lapack_int *lda, const float *anorm, float *rcond,
                 float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SPOCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spocon)(layout, *uplo, *n, a_r, lda_r, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_spocon_work)(layout, *uplo, *n, a_r, lda_r, *anorm,
                                          rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("SPOCON", ret);
}

/******************************************************************************
 * SPOEQU( N, A, LDA, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define SPOEQU_TEST LAPACK_GLOBAL_SUFFIX(spoequ_test, SPOEQU_TEST)
void SPOEQU_TEST(const lapack_int *n, const float *a, const lapack_int *lda,
                 float *s, float *scond, float *amax, lapack_int *info)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SPOEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spoequ)(layout, *n, a_r, lda_r, s, scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_spoequ_work)(layout, *n, a_r, lda_r, s, scond,
                                          amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("SPOEQU", ret);
}

/******************************************************************************
 * SPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX, FERR, BERR, WORK,
 * IWORK, INFO )
 ******************************************************************************/
#define SPORFS_TEST LAPACK_GLOBAL_SUFFIX(sporfs_test, SPORFS_TEST)
void SPORFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *a, const lapack_int *lda, const float *af,
                 const lapack_int *ldaf, const float *b, const lapack_int *ldb,
                 float *x, const lapack_int *ldx, float *ferr, float *berr,
                 float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    const float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_spo_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SPORFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sporfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_sporfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, b_r, ldb_r, x_r, ldx_r,
                                          ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SPORFS", ret);
}

/******************************************************************************
 * SPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define SPOSV_TEST LAPACK_GLOBAL_SUFFIX(sposv_test, SPOSV_TEST)
void SPOSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                float *a, const lapack_int *lda, float *b,
                const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SPOSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sposv)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sposv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SPOSV", ret);
}

/******************************************************************************
 * SPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S, B, LDB, X, LDX,
 * RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define SPOSVX_TEST LAPACK_GLOBAL_SUFFIX(sposvx_test, SPOSVX_TEST)
void SPOSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, float *a, const lapack_int *lda,
                 float *af, const lapack_int *ldaf, char *equed, float *s,
                 float *b, const lapack_int *ldb, float *x,
                 const lapack_int *ldx, float *rcond, float *ferr, float *berr,
                 float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
    float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_spo_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(af_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SPOSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sposvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, equed, s, b_r, ldb_r,
                                     x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_sposvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, equed, s,
        b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_spo_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(af_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SPOSVX", ret);
}

/******************************************************************************
 * SPOTRF( UPLO, N, A, LDA, INFO )
 ******************************************************************************/
#define SPOTRF_TEST LAPACK_GLOBAL_SUFFIX(spotrf_test, SPOTRF_TEST)
void SPOTRF_TEST(const char *uplo, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SPOTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spotrf)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_spotrf_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SPOTRF", ret);
}

/******************************************************************************
 * SPOTRI( UPLO, N, A, LDA, INFO )
 ******************************************************************************/
#define SPOTRI_TEST LAPACK_GLOBAL_SUFFIX(spotri_test, SPOTRI_TEST)
void SPOTRI_TEST(const char *uplo, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SPOTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spotri)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_spotri_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SPOTRI", ret);
}

/******************************************************************************
 * SPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define SPOTRS_TEST LAPACK_GLOBAL_SUFFIX(spotrs_test, SPOTRS_TEST)
void SPOTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *a, const lapack_int *lda, float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SPOTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spotrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_spotrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SPOTRS", ret);
}

/******************************************************************************
 * SPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
 ******************************************************************************/
#define SPSTRF_TEST LAPACK_GLOBAL_SUFFIX(spstrf_test, SPSTRF_TEST)
void SPSTRF_TEST(const char *uplo, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *piv, lapack_int *rank,
                 const float *tol, float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SPSTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spstrf)(layout, *uplo, *n, a_r, lda_r, piv, rank,
                                     *tol);
#else
    ret = API_SUFFIX(LAPACKE_spstrf_work)(layout, *uplo, *n, a_r, lda_r, piv,
                                          rank, *tol, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SPSTRF", ret);
}

// ========================================================================== //
//                   Packed positive definite matrices (PP)                   //
// ========================================================================== //

/******************************************************************************
 * SPPCON( UPLO, N, AP, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define SPPCON_TEST LAPACK_GLOBAL_SUFFIX(sppcon_test, SPPCON_TEST)
void SPPCON_TEST(const char *uplo, const lapack_int *n, const float *ap,
                 const float *anorm, float *rcond, float *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_spp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("SPPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sppcon)(layout, *uplo, *n, ap_r, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_sppcon_work)(layout, *uplo, *n, ap_r, *anorm,
                                          rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("SPPCON", ret);
}

/******************************************************************************
 * SPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define SPPEQU_TEST LAPACK_GLOBAL_SUFFIX(sppequ_test, SPPEQU_TEST)
void SPPEQU_TEST(const char *uplo, const lapack_int *n, const float *ap,
                 float *s, float *scond, float *amax, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_spp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("SPPEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sppequ)(layout, *uplo, *n, ap_r, s, scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_sppequ_work)(layout, *uplo, *n, ap_r, s, scond,
                                          amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("SPPEQU", ret);
}

/******************************************************************************
 * SPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO
 * )
 ******************************************************************************/
#define SPPRFS_TEST LAPACK_GLOBAL_SUFFIX(spprfs_test, SPPRFS_TEST)
void SPPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *ap, const float *afp, const float *b,
                 const lapack_int *ldb, float *x, const lapack_int *ldx,
                 float *ferr, float *berr, float *work, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
    const float *ap_r = ap;
    const float *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_spp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_spp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("SPPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r, b_r,
                                     ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_spprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          b_r, ldb_r, x_r, ldx_r, ferr, berr,
                                          work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("SPPRFS", ret);
}

/******************************************************************************
 * SPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
 ******************************************************************************/
#define SPPSV_TEST LAPACK_GLOBAL_SUFFIX(sppsv_test, SPPSV_TEST)
void SPPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                float *ap, float *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_spp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("SPPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sppsv)(layout, *uplo, *n, *nrhs, ap_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sppsv_work)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                         ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_spp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("SPPSV", ret);
}

/******************************************************************************
 * SPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define SPPSVX_TEST LAPACK_GLOBAL_SUFFIX(sppsvx_test, SPPSVX_TEST)
void SPPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, float *ap, float *afp, char *equed,
                 float *s, float *b, const lapack_int *ldb, float *x,
                 const lapack_int *ldx, float *rcond, float *ferr, float *berr,
                 float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *ap_r = ap;
    float *afp_r = afp;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_spp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_spp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SPPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sppsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, equed, s, b_r, ldb_r, x_r, ldx_r,
                                     rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_sppsvx_work)(
        layout, *fact, *uplo, *n, *nrhs, ap_r, afp_r, equed, s, b_r, ldb_r, x_r,
        ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_spp_rm_to_cm(*uplo, *n, ap_r, ap);
    lapacke_test_spp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SPPSVX", ret);
}

/******************************************************************************
 * SPPTRF( UPLO, N, AP, INFO )
 ******************************************************************************/
#define SPPTRF_TEST LAPACK_GLOBAL_SUFFIX(spptrf_test, SPPTRF_TEST)
void SPPTRF_TEST(const char *uplo, const lapack_int *n, float *ap,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_spp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("SPPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spptrf)(layout, *uplo, *n, ap_r);
#else
    ret = API_SUFFIX(LAPACKE_spptrf_work)(layout, *uplo, *n, ap_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("SPPTRF", ret);
}

/******************************************************************************
 * SPPTRI( UPLO, N, AP, INFO )
 ******************************************************************************/
#define SPPTRI_TEST LAPACK_GLOBAL_SUFFIX(spptri_test, SPPTRI_TEST)
void SPPTRI_TEST(const char *uplo, const lapack_int *n, float *ap,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_spp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("SPPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spptri)(layout, *uplo, *n, ap_r);
#else
    ret = API_SUFFIX(LAPACKE_spptri_work)(layout, *uplo, *n, ap_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("SPPTRI", ret);
}

/******************************************************************************
 * SPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
 ******************************************************************************/
#define SPPTRS_TEST LAPACK_GLOBAL_SUFFIX(spptrs_test, SPPTRS_TEST)
void SPPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *ap, float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
    const float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_spp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("SPPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spptrs)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_spptrs_work)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                          ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("SPPTRS", ret);
}

// ========================================================================== //
//                    Positive definite band matrices (PB)                    //
// ========================================================================== //

/******************************************************************************
 * SPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define SPBCON_TEST LAPACK_GLOBAL_SUFFIX(spbcon_test, SPBCON_TEST)
void SPBCON_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const float *ab, const lapack_int *ldab, const float *anorm,
                 float *rcond, float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_spb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("SPBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spbcon)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_spbcon_work)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                          *anorm, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("SPBCON", ret);
}

/******************************************************************************
 * SPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define SPBEQU_TEST LAPACK_GLOBAL_SUFFIX(spbequ_test, SPBEQU_TEST)
void SPBEQU_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const float *ab, const lapack_int *ldab, float *s,
                 float *scond, float *amax, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_spb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("SPBEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spbequ)(layout, *uplo, *n, *kd, ab_r, ldab_r, s,
                                     scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_spbequ_work)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                          s, scond, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("SPBEQU", ret);
}

/******************************************************************************
 * SPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR,
 * WORK, IWORK, INFO )
 ******************************************************************************/
#define SPBRFS_TEST LAPACK_GLOBAL_SUFFIX(spbrfs_test, SPBRFS_TEST)
void SPBRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const float *ab,
                 const lapack_int *ldab, const float *afb,
                 const lapack_int *ldafb, const float *b, const lapack_int *ldb,
                 float *x, const lapack_int *ldx, float *ferr, float *berr,
                 float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const float *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_spb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_spb_cm_to_rm(*uplo, *n, *kd, afb, *ldafb, &ldafb_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)afb_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SPBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spbrfs)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, b_r, ldb_r, x_r,
                                     ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_spbrfs_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                          ldab_r, afb_r, ldafb_r, b_r, ldb_r,
                                          x_r, ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)afb_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SPBRFS", ret);
}

/******************************************************************************
 * SPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define SPBSV_TEST LAPACK_GLOBAL_SUFFIX(spbsv_test, SPBSV_TEST)
void SPBSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                const lapack_int *nrhs, float *ab, const lapack_int *ldab,
                float *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_spb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SPBSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spbsv)(layout, *uplo, *n, *kd, *nrhs, ab_r, ldab_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_spbsv_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                         ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SPBSV", ret);
}

/******************************************************************************
 * SPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define SPBSVX_TEST LAPACK_GLOBAL_SUFFIX(spbsvx_test, SPBSVX_TEST)
void SPBSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *kd, const lapack_int *nrhs, float *ab,
                 const lapack_int *ldab, float *afb, const lapack_int *ldafb,
                 char *equed, float *s, float *b, const lapack_int *ldb,
                 float *x, const lapack_int *ldx, float *rcond, float *ferr,
                 float *berr, float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    float *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_spb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_spb_cm_to_rm(*uplo, *n, *kd, afb, *ldafb, &ldafb_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(afb_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SPBSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spbsvx)(layout, *fact, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, equed, s, b_r,
                                     ldb_r, x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_spbsvx_work)(
        layout, *fact, *uplo, *n, *kd, *nrhs, ab_r, ldab_r, afb_r, ldafb_r,
        equed, s, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    lapacke_test_spb_rm_to_cm(*uplo, *n, *kd, afb_r, ldafb_r, afb, *ldafb);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(ab_r);
    LAPACKE_free(afb_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SPBSVX", ret);
}

/******************************************************************************
 * SPBTRF( UPLO, N, KD, AB, LDAB, INFO )
 ******************************************************************************/
#define SPBTRF_TEST LAPACK_GLOBAL_SUFFIX(spbtrf_test, SPBTRF_TEST)
void SPBTRF_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 float *ab, const lapack_int *ldab, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_spb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("SPBTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spbtrf)(layout, *uplo, *n, *kd, ab_r, ldab_r);
#else
    ret = API_SUFFIX(LAPACKE_spbtrf_work)(layout, *uplo, *n, *kd, ab_r, ldab_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    LAPACKE_free(ab_r);
#endif
    *info = lapacke_test_info("SPBTRF", ret);
}

/******************************************************************************
 * SPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define SPBTRS_TEST LAPACK_GLOBAL_SUFFIX(spbtrs_test, SPBTRS_TEST)
void SPBTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const float *ab,
                 const lapack_int *ldab, float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_spb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SPBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spbtrs)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_spbtrs_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                          ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SPBTRS", ret);
}

// ========================================================================== //
//                 Positive definite tridiagonal matrices (PT)                //
// ========================================================================== //

/******************************************************************************
 * SPTCON( N, D, E, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define SPTCON_TEST LAPACK_GLOBAL_SUFFIX(sptcon_test, SPTCON_TEST)
void SPTCON_TEST(const lapack_int *n, const float *d, const float *e,
                 const float *anorm, float *rcond, float *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sptcon)(*n, d, e, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_sptcon_work)(*n, d, e, *anorm, rcond, work);
#endif

    *info = lapacke_test_info_unshifted("SPTCON", ret);
}

/******************************************************************************
 * SPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, INFO )
 ******************************************************************************/
#define SPTRFS_TEST LAPACK_GLOBAL_SUFFIX(sptrfs_test, SPTRFS_TEST)
void SPTRFS_TEST(const lapack_int *n, const lapack_int *nrhs, const float *d,
                 const float *e, const float *df, const float *ef,
                 const float *b, const lapack_int *ldb, float *x,
                 const lapack_int *ldx, float *ferr, float *berr, float *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SPTRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sptrfs)(layout, *n, *nrhs, d, e, df, ef, b_r,
                                     ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_sptrfs_work)(layout, *n, *nrhs, d, e, df, ef, b_r,
                                          ldb_r, x_r, ldx_r, ferr, berr, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SPTRFS", ret);
}

/******************************************************************************
 * SPTSV( N, NRHS, D, E, B, LDB, INFO )
 ******************************************************************************/
#define SPTSV_TEST LAPACK_GLOBAL_SUFFIX(sptsv_test, SPTSV_TEST)
void SPTSV_TEST(const lapack_int *n, const lapack_int *nrhs, float *d, float *e,
                float *b, const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("SPTSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sptsv)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sptsv_work)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SPTSV", ret);
}

/******************************************************************************
 * SPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK,
 * INFO )
 ******************************************************************************/
#define SPTSVX_TEST LAPACK_GLOBAL_SUFFIX(sptsvx_test, SPTSVX_TEST)
void SPTSVX_TEST(const char *fact, const lapack_int *n, const lapack_int *nrhs,
                 const float *d, const float *e, float *df, float *ef,
                 const float *b, const lapack_int *ldb, float *x,
                 const lapack_int *ldx, float *rcond, float *ferr, float *berr,
                 float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len
#endif
)
{
    lapack_int ret = 0;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SPTSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sptsvx)(layout, *fact, *n, *nrhs, d, e, df, ef,
                                     b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_sptsvx_work)(layout, *fact, *n, *nrhs, d, e, df,
                                          ef, b_r, ldb_r, x_r, ldx_r, rcond,
                                          ferr, berr, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SPTSVX", ret);
}

/******************************************************************************
 * SPTTRF( N, D, E, INFO )
 ******************************************************************************/
#define SPTTRF_TEST LAPACK_GLOBAL_SUFFIX(spttrf_test, SPTTRF_TEST)
void SPTTRF_TEST(const lapack_int *n, float *d, float *e, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spttrf)(*n, d, e);
#else
    ret = API_SUFFIX(LAPACKE_spttrf_work)(*n, d, e);
#endif

    *info = lapacke_test_info_unshifted("SPTTRF", ret);
}

/******************************************************************************
 * SPTTRS( N, NRHS, D, E, B, LDB, INFO )
 ******************************************************************************/
#define SPTTRS_TEST LAPACK_GLOBAL_SUFFIX(spttrs_test, SPTTRS_TEST)
void SPTTRS_TEST(const lapack_int *n, const lapack_int *nrhs, const float *d,
                 const float *e, float *b, const lapack_int *ldb,
                 lapack_int *info)
{
    lapack_int ret = 0;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("SPTTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spttrs)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_spttrs_work)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SPTTRS", ret);
}

// ========================================================================== //
//                     Symmetric indefinite matrices (SY)                     //
// ========================================================================== //

/******************************************************************************
 * SSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define SSYCON_TEST LAPACK_GLOBAL_SUFFIX(ssycon_test, SSYCON_TEST)
void SSYCON_TEST(const char *uplo, const lapack_int *n, const float *a,
                 const lapack_int *lda, const lapack_int *ipiv,
                 const float *anorm, float *rcond, float *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssycon)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_ssycon_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          *anorm, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("SSYCON", ret);
}

/******************************************************************************
 * SSYCON_3( UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define SSYCON_3_TEST LAPACK_GLOBAL_SUFFIX(ssycon_3_test, SSYCON_3_TEST)
void SSYCON_3_TEST(const char *uplo, const lapack_int *n, const float *a,
                   const lapack_int *lda, const float *e,
                   const lapack_int *ipiv, const float *anorm, float *rcond,
                   float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYCON_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssycon_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv,
                                       *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_ssycon_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, *anorm, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("SSYCON_3", ret);
}

/******************************************************************************
 * SSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, IWORK, INFO )
 ******************************************************************************/
#define SSYRFS_TEST LAPACK_GLOBAL_SUFFIX(ssyrfs_test, SSYRFS_TEST)
void SSYRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *a, const lapack_int *lda, const float *af,
                 const lapack_int *ldaf, const lapack_int *ipiv, const float *b,
                 const lapack_int *ldb, float *x, const lapack_int *ldx,
                 float *ferr, float *berr, float *work, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    const float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SSYRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssyrfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_ssyrfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SSYRFS", ret);
}

/******************************************************************************
 * SSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYSV_TEST LAPACK_GLOBAL_SUFFIX(ssysv_test, SSYSV_TEST)
void SSYSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                float *a, const lapack_int *lda, lapack_int *ipiv, float *b,
                const lapack_int *ldb, float *work, const lapack_int *lwork,
                lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssysv_work)(LAPACK_COL_MAJOR, *uplo, *n, *nrhs,
                                             a, *lda, ipiv, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("SSYSV", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SSYSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssysv)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssysv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SSYSV", ret);
}

/******************************************************************************
 * SSYSV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYSV_AA_TEST LAPACK_GLOBAL_SUFFIX(ssysv_aa_test, SSYSV_AA_TEST)
void SSYSV_AA_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, float *a, const lapack_int *lda,
                   lapack_int *ipiv, float *b, const lapack_int *ldb,
                   float *work, const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssysv_aa_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, ipiv, b, *ldb,
                                                work, *lwork);
        *info = lapacke_test_info("SSYSV_AA", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SSYSV_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssysv_aa)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssysv_aa_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SSYSV_AA", ret);
}

/******************************************************************************
 * SSYSV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYSV_RK_TEST LAPACK_GLOBAL_SUFFIX(ssysv_rk_test, SSYSV_RK_TEST)
void SSYSV_RK_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, float *a, const lapack_int *lda,
                   float *e, lapack_int *ipiv, float *b, const lapack_int *ldb,
                   float *work, const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssysv_rk_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, e, ipiv, b,
                                                *ldb, work, *lwork);
        *info = lapacke_test_info("SSYSV_RK", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SSYSV_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssysv_rk)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssysv_rk_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r, work,
                                            *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SSYSV_RK", ret);
}

/******************************************************************************
 * SSYSV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYSV_ROOK_TEST LAPACK_GLOBAL_SUFFIX(ssysv_rook_test, SSYSV_ROOK_TEST)
void SSYSV_ROOK_TEST(const char *uplo, const lapack_int *n,
                     const lapack_int *nrhs, float *a, const lapack_int *lda,
                     lapack_int *ipiv, float *b, const lapack_int *ldb,
                     float *work, const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                     ,
                     FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssysv_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                  *nrhs, a, *lda, ipiv, b, *ldb,
                                                  work, *lwork);
        *info = lapacke_test_info("SSYSV_ROOK", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SSYSV_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssysv_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssysv_rook_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SSYSV_ROOK", ret);
}

/******************************************************************************
 * SSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND,
 * FERR, BERR, WORK, LWORK, IWORK, INFO )
 ******************************************************************************/
#define SSYSVX_TEST LAPACK_GLOBAL_SUFFIX(ssysvx_test, SSYSVX_TEST)
void SSYSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const float *a, const lapack_int *lda,
                 float *af, const lapack_int *ldaf, lapack_int *ipiv,
                 const float *b, const lapack_int *ldb, float *x,
                 const lapack_int *ldx, float *rcond, float *ferr, float *berr,
                 float *work, const lapack_int *lwork, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssysvx_work)(
            LAPACK_COL_MAJOR, *fact, *uplo, *n, *nrhs, a, *lda, af, *ldaf, ipiv,
            b, *ldb, x, *ldx, rcond, ferr, berr, work, *lwork, iwork);
        *info = lapacke_test_info("SSYSVX", ret);
        return;
    }

    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SSYSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssysvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                     ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_ssysvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, b_r,
        ldb_r, x_r, ldx_r, rcond, ferr, berr, work, *lwork, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SSYSVX", ret);
}

/******************************************************************************
 * SSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYTRF_TEST LAPACK_GLOBAL_SUFFIX(ssytrf_test, SSYTRF_TEST)
void SSYTRF_TEST(const char *uplo, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *ipiv, float *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssytrf_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                              *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("SSYTRF", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytrf)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_ssytrf_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SSYTRF", ret);
}

/******************************************************************************
 * SSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYTRF_AA_TEST LAPACK_GLOBAL_SUFFIX(ssytrf_aa_test, SSYTRF_AA_TEST)
void SSYTRF_AA_TEST(const char *uplo, const lapack_int *n, float *a,
                    const lapack_int *lda, lapack_int *ipiv, float *work,
                    const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssytrf_aa_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("SSYTRF_AA", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYTRF_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytrf_aa)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_ssytrf_aa_work)(layout, *uplo, *n, a_r, lda_r,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SSYTRF_AA", ret);
}

/******************************************************************************
 * SSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYTRF_RK_TEST LAPACK_GLOBAL_SUFFIX(ssytrf_rk_test, SSYTRF_RK_TEST)
void SSYTRF_RK_TEST(const char *uplo, const lapack_int *n, float *a,
                    const lapack_int *lda, float *e, lapack_int *ipiv,
                    float *work, const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssytrf_rk_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("SSYTRF_RK", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYTRF_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytrf_rk)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_ssytrf_rk_work)(layout, *uplo, *n, a_r, lda_r, e,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SSYTRF_RK", ret);
}

/******************************************************************************
 * SSYTRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYTRF_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(ssytrf_rook_test, SSYTRF_ROOK_TEST)
void SSYTRF_ROOK_TEST(const char *uplo, const lapack_int *n, float *a,
                      const lapack_int *lda, lapack_int *ipiv, float *work,
                      const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssytrf_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                   a, *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("SSYTRF_ROOK", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYTRF_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytrf_rook)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_ssytrf_rook_work)(layout, *uplo, *n, a_r, lda_r,
                                               ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SSYTRF_ROOK", ret);
}

/******************************************************************************
 * SSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
 ******************************************************************************/
#define SSYTRI_TEST LAPACK_GLOBAL_SUFFIX(ssytri_test, SSYTRI_TEST)
void SSYTRI_TEST(const char *uplo, const lapack_int *n, float *a,
                 const lapack_int *lda, const lapack_int *ipiv, float *work,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytri)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_ssytri_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SSYTRI", ret);
}

/******************************************************************************
 * SSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYTRI2_TEST LAPACK_GLOBAL_SUFFIX(ssytri2_test, SSYTRI2_TEST)
void SSYTRI2_TEST(const char *uplo, const lapack_int *n, float *a,
                  const lapack_int *lda, const lapack_int *ipiv, float *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssytri2_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                               *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("SSYTRI2", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYTRI2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytri2)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_ssytri2_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SSYTRI2", ret);
}

/******************************************************************************
 * SSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
 ******************************************************************************/
#define SSYTRI2X_TEST LAPACK_GLOBAL_SUFFIX(ssytri2x_test, SSYTRI2X_TEST)
void SSYTRI2X_TEST(const char *uplo, const lapack_int *n, float *a,
                   const lapack_int *lda, const lapack_int *ipiv, float *work,
                   const lapack_int *nb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYTRI2X", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytri2x)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                       *nb);
#else
    ret = API_SUFFIX(LAPACKE_ssytri2x_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                            work, *nb);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SSYTRI2X", ret);
}

/******************************************************************************
 * SSYTRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYTRI_3_TEST LAPACK_GLOBAL_SUFFIX(ssytri_3_test, SSYTRI_3_TEST)
void SSYTRI_3_TEST(const char *uplo, const lapack_int *n, float *a,
                   const lapack_int *lda, const float *e,
                   const lapack_int *ipiv, float *work, const lapack_int *lwork,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssytri_3_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("SSYTRI_3", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SSYTRI_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytri_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_ssytri_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SSYTRI_3", ret);
}

/******************************************************************************
 * SSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SSYTRS_TEST LAPACK_GLOBAL_SUFFIX(ssytrs_test, SSYTRS_TEST)
void SSYTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *a, const lapack_int *lda, const lapack_int *ipiv,
                 float *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SSYTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                     b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssytrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SSYTRS", ret);
}

/******************************************************************************
 * SSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )
 ******************************************************************************/
#define SSYTRS2_TEST LAPACK_GLOBAL_SUFFIX(ssytrs2_test, SSYTRS2_TEST)
void SSYTRS2_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                  const float *a, const lapack_int *lda, const lapack_int *ipiv,
                  float *b, const lapack_int *ldb, float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SSYTRS2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytrs2)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                      ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssytrs2_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                           ipiv, b_r, ldb_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SSYTRS2", ret);
}

/******************************************************************************
 * SSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SSYTRS_3_TEST LAPACK_GLOBAL_SUFFIX(ssytrs_3_test, SSYTRS_3_TEST)
void SSYTRS_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, const float *a,
                   const lapack_int *lda, const float *e,
                   const lapack_int *ipiv, float *b, const lapack_int *ldb,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SSYTRS_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytrs_3)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssytrs_3_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SSYTRS_3", ret);
}

/******************************************************************************
 * SSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define SSYTRS_AA_TEST LAPACK_GLOBAL_SUFFIX(ssytrs_aa_test, SSYTRS_AA_TEST)
void SSYTRS_AA_TEST(const char *uplo, const lapack_int *n,
                    const lapack_int *nrhs, const float *a,
                    const lapack_int *lda, const lapack_int *ipiv, float *b,
                    const lapack_int *ldb, float *work, const lapack_int *lwork,
                    lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_ssytrs_aa_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                 *nrhs, a, *lda, ipiv, b, *ldb,
                                                 work, *lwork);
        *info = lapacke_test_info("SSYTRS_AA", ret);
        return;
    }

    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SSYTRS_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytrs_aa)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                        ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssytrs_aa_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SSYTRS_AA", ret);
}

/******************************************************************************
 * SSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SSYTRS_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(ssytrs_rook_test, SSYTRS_ROOK_TEST)
void SSYTRS_ROOK_TEST(const char *uplo, const lapack_int *n,
                      const lapack_int *nrhs, const float *a,
                      const lapack_int *lda, const lapack_int *ipiv, float *b,
                      const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SSYTRS_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssytrs_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssytrs_rook_work)(layout, *uplo, *n, *nrhs, a_r,
                                               lda_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SSYTRS_ROOK", ret);
}

// ========================================================================== //
//                  Packed symmetric indefinite matrices (SP)                 //
// ========================================================================== //

/******************************************************************************
 * SSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define SSPCON_TEST LAPACK_GLOBAL_SUFFIX(sspcon_test, SSPCON_TEST)
void SSPCON_TEST(const char *uplo, const lapack_int *n, const float *ap,
                 const lapack_int *ipiv, const float *anorm, float *rcond,
                 float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_ssp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("SSPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sspcon)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_sspcon_work)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                          rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("SSPCON", ret);
}

/******************************************************************************
 * SSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK,
 * IWORK, INFO )
 ******************************************************************************/
#define SSPRFS_TEST LAPACK_GLOBAL_SUFFIX(ssprfs_test, SSPRFS_TEST)
void SSPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *ap, const float *afp, const lapack_int *ipiv,
                 const float *b, const lapack_int *ldb, float *x,
                 const lapack_int *ldx, float *ferr, float *berr, float *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
    const float *ap_r = ap;
    const float *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_ssp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_ssp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("SSPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                     ipiv, b_r, ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_ssprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                          berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("SSPRFS", ret);
}

/******************************************************************************
 * SSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SSPSV_TEST LAPACK_GLOBAL_SUFFIX(sspsv_test, SSPSV_TEST)
void SSPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                float *ap, lapack_int *ipiv, float *b, const lapack_int *ldb,
                lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
    float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_ssp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("SSPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sspsv)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sspsv_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_ssp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("SSPSV", ret);
}

/******************************************************************************
 * SSPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define SSPSVX_TEST LAPACK_GLOBAL_SUFFIX(sspsvx_test, SSPSVX_TEST)
void SSPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const float *ap, float *afp,
                 lapack_int *ipiv, const float *b, const lapack_int *ldb,
                 float *x, const lapack_int *ldx, float *rcond, float *ferr,
                 float *berr, float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    const float *ap_r = ap;
    float *afp_r = afp;
    float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_ssp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_ssp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("SSPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sspsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, ipiv, b_r, ldb_r, x_r, ldx_r, rcond,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_sspsvx_work)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                          afp_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                          rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_sge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("SSPSVX", ret);
}

/******************************************************************************
 * SSPTRF( UPLO, N, AP, IPIV, INFO )
 ******************************************************************************/
#define SSPTRF_TEST LAPACK_GLOBAL_SUFFIX(ssptrf_test, SSPTRF_TEST)
void SSPTRF_TEST(const char *uplo, const lapack_int *n, float *ap,
                 lapack_int *ipiv, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_ssp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("SSPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssptrf)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_ssptrf_work)(layout, *uplo, *n, ap_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("SSPTRF", ret);
}

/******************************************************************************
 * SSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
 ******************************************************************************/
#define SSPTRI_TEST LAPACK_GLOBAL_SUFFIX(ssptri_test, SSPTRI_TEST)
void SSPTRI_TEST(const char *uplo, const lapack_int *n, float *ap,
                 const lapack_int *ipiv, float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_ssp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("SSPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssptri)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_ssptri_work)(layout, *uplo, *n, ap_r, ipiv, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_ssp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("SSPTRI", ret);
}

/******************************************************************************
 * SSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define SSPTRS_TEST LAPACK_GLOBAL_SUFFIX(ssptrs_test, SSPTRS_TEST)
void SSPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *ap, const lapack_int *ipiv, float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
    const float *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_ssp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("SSPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_ssptrs)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_ssptrs_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("SSPTRS", ret);
}

// ========================================================================== //
//                          Triangular matrices (TR)                          //
// ========================================================================== //

/******************************************************************************
 * STRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define STRCON_TEST LAPACK_GLOBAL_SUFFIX(strcon_test, STRCON_TEST)
void STRCON_TEST(const char *norm, const char *uplo, const char *diag,
                 const lapack_int *n, const float *a, const lapack_int *lda,
                 float *rcond, float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_str_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("STRCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_strcon)(layout, *norm, *uplo, *diag, *n, a_r,
                                     lda_r, rcond);
#else
    ret = API_SUFFIX(LAPACKE_strcon_work)(layout, *norm, *uplo, *diag, *n, a_r,
                                          lda_r, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("STRCON", ret);
}

/******************************************************************************
 * STRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, LDX, FERR, BERR, WORK,
 * IWORK, INFO )
 ******************************************************************************/
#define STRRFS_TEST LAPACK_GLOBAL_SUFFIX(strrfs_test, STRRFS_TEST)
void STRRFS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *nrhs, const float *a,
                 const lapack_int *lda, const float *b, const lapack_int *ldb,
                 const float *x, const lapack_int *ldx, float *ferr,
                 float *berr, float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    const float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_str_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)x_r);
        lapacke_test_report_alloc_failure("STRRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_strrfs)(layout, *uplo, *trans, *diag, *n, *nrhs,
                                     a_r, lda_r, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_strrfs_work)(layout, *uplo, *trans, *diag, *n,
                                          *nrhs, a_r, lda_r, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)x_r);
#endif
    *info = lapacke_test_info("STRRFS", ret);
}

/******************************************************************************
 * STRTRI( UPLO, DIAG, N, A, LDA, INFO )
 ******************************************************************************/
#define STRTRI_TEST LAPACK_GLOBAL_SUFFIX(strtri_test, STRTRI_TEST)
void STRTRI_TEST(const char *uplo, const char *diag, const lapack_int *n,
                 float *a, const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_str_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("STRTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_strtri)(layout, *uplo, *diag, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_strtri_work)(layout, *uplo, *diag, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_str_rm_to_cm(*uplo, *diag, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("STRTRI", ret);
}

/******************************************************************************
 * STRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define STRTRS_TEST LAPACK_GLOBAL_SUFFIX(strtrs_test, STRTRS_TEST)
void STRTRS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *nrhs, const float *a,
                 const lapack_int *lda, float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_str_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("STRTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_strtrs)(layout, *uplo, *trans, *diag, *n, *nrhs,
                                     a_r, lda_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_strtrs_work)(layout, *uplo, *trans, *diag, *n,
                                          *nrhs, a_r, lda_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("STRTRS", ret);
}

// ========================================================================== //
//                        Triangular band matrices (TB)                       //
// ========================================================================== //

/******************************************************************************
 * STBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define STBCON_TEST LAPACK_GLOBAL_SUFFIX(stbcon_test, STBCON_TEST)
void STBCON_TEST(const char *norm, const char *uplo, const char *diag,
                 const lapack_int *n, const lapack_int *kd, const float *ab,
                 const lapack_int *ldab, float *rcond, float *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_stb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("STBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_stbcon)(layout, *norm, *uplo, *diag, *n, *kd, ab_r,
                                     ldab_r, rcond);
#else
    ret = API_SUFFIX(LAPACKE_stbcon_work)(layout, *norm, *uplo, *diag, *n, *kd,
                                          ab_r, ldab_r, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("STBCON", ret);
}

/******************************************************************************
 * STBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, X, LDX, FERR, BERR,
 * WORK, IWORK, INFO )
 ******************************************************************************/
#define STBRFS_TEST LAPACK_GLOBAL_SUFFIX(stbrfs_test, STBRFS_TEST)
void STBRFS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const float *ab,
                 const lapack_int *ldab, const float *b, const lapack_int *ldb,
                 const float *x, const lapack_int *ldx, float *ferr,
                 float *berr, float *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const float *b_r = b;
    lapack_int ldb_r = *ldb;
    const float *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_stb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)x_r);
        lapacke_test_report_alloc_failure("STBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_stbrfs)(layout, *uplo, *trans, *diag, *n, *kd,
                                     *nrhs, ab_r, ldab_r, b_r, ldb_r, x_r,
                                     ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_stbrfs_work)(layout, *uplo, *trans, *diag, *n, *kd,
                                          *nrhs, ab_r, ldab_r, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)x_r);
#endif
    *info = lapacke_test_info("STBRFS", ret);
}

/******************************************************************************
 * STBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define STBTRS_TEST LAPACK_GLOBAL_SUFFIX(stbtrs_test, STBTRS_TEST)
void STBTRS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const float *ab,
                 const lapack_int *ldab, float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const float *ab_r = ab;
    lapack_int ldab_r = *ldab;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_stb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("STBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_stbtrs)(layout, *uplo, *trans, *diag, *n, *kd,
                                     *nrhs, ab_r, ldab_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_stbtrs_work)(layout, *uplo, *trans, *diag, *n, *kd,
                                          *nrhs, ab_r, ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("STBTRS", ret);
}

// ========================================================================== //
//                   Orthogonal factorizations (QR/LQ/QL/RQ)                  //
// ========================================================================== //

/******************************************************************************
 * SGELQ2( M, N, A, LDA, TAU, WORK, INFO )
 ******************************************************************************/
#define SGELQ2_TEST LAPACK_GLOBAL_SUFFIX(sgelq2_test, SGELQ2_TEST)
void SGELQ2_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, float *tau, float *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGELQ2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgelq2)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sgelq2_work)(layout, *m, *n, a_r, lda_r, tau,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGELQ2", ret);
}

/******************************************************************************
 * SGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGELQF_TEST LAPACK_GLOBAL_SUFFIX(sgelqf_test, SGELQF_TEST)
void SGELQF_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgelqf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("SGELQF", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGELQF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgelqf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sgelqf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGELQF", ret);
}

/******************************************************************************
 * SGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGEQLF_TEST LAPACK_GLOBAL_SUFFIX(sgeqlf_test, SGEQLF_TEST)
void SGEQLF_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgeqlf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("SGEQLF", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGEQLF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqlf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sgeqlf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGEQLF", ret);
}

/******************************************************************************
 * SGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGEQR_TEST LAPACK_GLOBAL_SUFFIX(sgeqr_test, SGEQR_TEST)
void SGEQR_TEST(const lapack_int *m, const lapack_int *n, float *a,
                const lapack_int *lda, float *t, const lapack_int *tsize,
                float *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*tsize == -1 || *lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgeqr_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                             t, *tsize, work, *lwork);
        *info = lapacke_test_info("SGEQR", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGEQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqr)(layout, *m, *n, a_r, lda_r, t, *tsize);
#else
    ret = API_SUFFIX(LAPACKE_sgeqr_work)(layout, *m, *n, a_r, lda_r, t, *tsize,
                                         work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGEQR", ret);
}

/******************************************************************************
 * SGEQR2( M, N, A, LDA, TAU, WORK, INFO )
 ******************************************************************************/
#define SGEQR2_TEST LAPACK_GLOBAL_SUFFIX(sgeqr2_test, SGEQR2_TEST)
void SGEQR2_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, float *tau, float *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGEQR2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqr2)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sgeqr2_work)(layout, *m, *n, a_r, lda_r, tau,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGEQR2", ret);
}

/******************************************************************************
 * SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGEQRF_TEST LAPACK_GLOBAL_SUFFIX(sgeqrf_test, SGEQRF_TEST)
void SGEQRF_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgeqrf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("SGEQRF", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGEQRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqrf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sgeqrf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGEQRF", ret);
}

/******************************************************************************
 * SGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGEQRFP_TEST LAPACK_GLOBAL_SUFFIX(sgeqrfp_test, SGEQRFP_TEST)
void SGEQRFP_TEST(const lapack_int *m, const lapack_int *n, float *a,
                  const lapack_int *lda, float *tau, float *work,
                  const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgeqrfp_work)(LAPACK_COL_MAJOR, *m, *n, a,
                                               *lda, tau, work, *lwork);
        *info = lapacke_test_info("SGEQRFP", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGEQRFP", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqrfp)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sgeqrfp_work)(layout, *m, *n, a_r, lda_r, tau,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGEQRFP", ret);
}

/******************************************************************************
 * SGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGERQF_TEST LAPACK_GLOBAL_SUFFIX(sgerqf_test, SGERQF_TEST)
void SGERQF_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgerqf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("SGERQF", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGERQF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgerqf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sgerqf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGERQF", ret);
}

/******************************************************************************
 * SORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SORGLQ_TEST LAPACK_GLOBAL_SUFFIX(sorglq_test, SORGLQ_TEST)
void SORGLQ_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 float *a, const lapack_int *lda, const float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sorglq_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("SORGLQ", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SORGLQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sorglq)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sorglq_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SORGLQ", ret);
}

/******************************************************************************
 * SORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SORGQL_TEST LAPACK_GLOBAL_SUFFIX(sorgql_test, SORGQL_TEST)
void SORGQL_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 float *a, const lapack_int *lda, const float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sorgql_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("SORGQL", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SORGQL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sorgql)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sorgql_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SORGQL", ret);
}

/******************************************************************************
 * SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SORGQR_TEST LAPACK_GLOBAL_SUFFIX(sorgqr_test, SORGQR_TEST)
void SORGQR_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 float *a, const lapack_int *lda, const float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sorgqr_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("SORGQR", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SORGQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sorgqr)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sorgqr_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SORGQR", ret);
}

/******************************************************************************
 * SORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SORGRQ_TEST LAPACK_GLOBAL_SUFFIX(sorgrq_test, SORGRQ_TEST)
void SORGRQ_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 float *a, const lapack_int *lda, const float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sorgrq_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("SORGRQ", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SORGRQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sorgrq)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sorgrq_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SORGRQ", ret);
}

/******************************************************************************
 * SORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define SORMLQ_TEST LAPACK_GLOBAL_SUFFIX(sormlq_test, SORMLQ_TEST)
void SORMLQ_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k, const float *a,
                 const lapack_int *lda, const float *tau, float *c,
                 const lapack_int *ldc, float *work, const lapack_int *lwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sormlq_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("SORMLQ", ret);
        return;
    }

    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_sge_cm_to_rm(*k, nq, a, *lda, &lda_r);
    c_r = lapacke_test_sge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("SORMLQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sormlq)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_sormlq_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("SORMLQ", ret);
}

/******************************************************************************
 * SORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define SORMQL_TEST LAPACK_GLOBAL_SUFFIX(sormql_test, SORMQL_TEST)
void SORMQL_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k, const float *a,
                 const lapack_int *lda, const float *tau, float *c,
                 const lapack_int *ldc, float *work, const lapack_int *lwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sormql_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("SORMQL", ret);
        return;
    }

    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_sge_cm_to_rm(nq, *k, a, *lda, &lda_r);
    c_r = lapacke_test_sge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("SORMQL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sormql)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_sormql_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("SORMQL", ret);
}

/******************************************************************************
 * SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define SORMQR_TEST LAPACK_GLOBAL_SUFFIX(sormqr_test, SORMQR_TEST)
void SORMQR_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k, const float *a,
                 const lapack_int *lda, const float *tau, float *c,
                 const lapack_int *ldc, float *work, const lapack_int *lwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sormqr_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("SORMQR", ret);
        return;
    }

    const float *a_r = a;
    lapack_int lda_r = *lda;
    float *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_sge_cm_to_rm(nq, *k, a, *lda, &lda_r);
    c_r = lapacke_test_sge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("SORMQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sormqr)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_sormqr_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("SORMQR", ret);
}

// ========================================================================== //
//            Blocked and tall-skinny factorizations (QRT/TSQR/HR)            //
// ========================================================================== //

/******************************************************************************
 * SGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
 ******************************************************************************/
#define SGEQRT_TEST LAPACK_GLOBAL_SUFFIX(sgeqrt_test, SGEQRT_TEST)
void SGEQRT_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *nb,
                 float *a, const lapack_int *lda, float *t,
                 const lapack_int *ldt, float *work, lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
    float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_sge_cm_to_rm(*nb, MIN(*m, *n), t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("SGEQRT", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqrt)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                     ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_sgeqrt_work)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                          ldt_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*nb, MIN(*m, *n), t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("SGEQRT", ret);
}

/******************************************************************************
 * SGEQRT2( M, N, A, LDA, T, LDT, INFO )
 ******************************************************************************/
#define SGEQRT2_TEST LAPACK_GLOBAL_SUFFIX(sgeqrt2_test, SGEQRT2_TEST)
void SGEQRT2_TEST(const lapack_int *m, const lapack_int *n, float *a,
                  const lapack_int *lda, float *t, const lapack_int *ldt,
                  lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
    float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_sge_cm_to_rm(*n, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("SGEQRT2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqrt2)(layout, *m, *n, a_r, lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_sgeqrt2_work)(layout, *m, *n, a_r, lda_r, t_r,
                                           ldt_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*n, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("SGEQRT2", ret);
}

/******************************************************************************
 * SGEQRT3( M, N, A, LDA, T, LDT, INFO )
 ******************************************************************************/
#define SGEQRT3_TEST LAPACK_GLOBAL_SUFFIX(sgeqrt3_test, SGEQRT3_TEST)
void SGEQRT3_TEST(const lapack_int *m, const lapack_int *n, float *a,
                  const lapack_int *lda, float *t, const lapack_int *ldt,
                  lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
    float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_sge_cm_to_rm(*n, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("SGEQRT3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqrt3)(layout, *m, *n, a_r, lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_sgeqrt3_work)(layout, *m, *n, a_r, lda_r, t_r,
                                           ldt_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*n, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("SGEQRT3", ret);
}

/******************************************************************************
 * SGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGETSQRHRT_TEST LAPACK_GLOBAL_SUFFIX(sgetsqrhrt_test, SGETSQRHRT_TEST)
void SGETSQRHRT_TEST(const lapack_int *m, const lapack_int *n,
                     const lapack_int *mb1, const lapack_int *nb1,
                     const lapack_int *nb2, float *a, const lapack_int *lda,
                     float *t, const lapack_int *ldt, float *work,
                     const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgetsqrhrt_work)(LAPACK_COL_MAJOR, *m, *n,
                                                  *mb1, *nb1, *nb2, a, *lda, t,
                                                  *ldt, work, *lwork);
        *info = lapacke_test_info("SGETSQRHRT", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_sge_cm_to_rm(*nb2, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("SGETSQRHRT", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgetsqrhrt)(layout, *m, *n, *mb1, *nb1, *nb2, a_r,
                                         lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_sgetsqrhrt_work)(
        layout, *m, *n, *mb1, *nb1, *nb2, a_r, lda_r, t_r, ldt_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*nb2, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("SGETSQRHRT", ret);
}

/******************************************************************************
 * SORHR_COL( M, N, NB, A, LDA, T, LDT, D, INFO )
 ******************************************************************************/
#define SORHR_COL_TEST LAPACK_GLOBAL_SUFFIX(sorhr_col_test, SORHR_COL_TEST)
void SORHR_COL_TEST(const lapack_int *m, const lapack_int *n,
                    const lapack_int *nb, float *a, const lapack_int *lda,
                    float *t, const lapack_int *ldt, float *d, lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
    float *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_sge_cm_to_rm(*ldt, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("SORHR_COL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sorhr_col)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                        ldt_r, d);
#else
    ret = API_SUFFIX(LAPACKE_sorhr_col_work)(layout, *m, *n, *nb, a_r, lda_r,
                                             t_r, ldt_r, d);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(*ldt, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("SORHR_COL", ret);
}

// ========================================================================== //
//                    Rank-revealing factorizations (QP/TZ)                   //
// ========================================================================== //

/******************************************************************************
 * SGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGEQP3_TEST LAPACK_GLOBAL_SUFFIX(sgeqp3_test, SGEQP3_TEST)
void SGEQP3_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *jpvt, float *tau,
                 float *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgeqp3_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              jpvt, tau, work, *lwork);
        *info = lapacke_test_info("SGEQP3", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGEQP3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqp3)(layout, *m, *n, a_r, lda_r, jpvt, tau);
#else
    ret = API_SUFFIX(LAPACKE_sgeqp3_work)(layout, *m, *n, a_r, lda_r, jpvt, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGEQP3", ret);
}

/******************************************************************************
 * STZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define STZRZF_TEST LAPACK_GLOBAL_SUFFIX(stzrzf_test, STZRZF_TEST)
void STZRZF_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_stzrzf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("STZRZF", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("STZRZF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_stzrzf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_stzrzf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("STZRZF", ret);
}

// ========================================================================== //
//                             Least squares (LS)                             //
// ========================================================================== //

/******************************************************************************
 * SGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGELS_TEST LAPACK_GLOBAL_SUFFIX(sgels_test, SGELS_TEST)
void SGELS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                const lapack_int *nrhs, float *a, const lapack_int *lda,
                float *b, const lapack_int *ldb, float *work,
                const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgels_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                             *nrhs, a, *lda, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("SGELS", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGELS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgels)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgels_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                         lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGELS", ret);
}

/******************************************************************************
 * SGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO
 * )
 ******************************************************************************/
#define SGELSD_TEST LAPACK_GLOBAL_SUFFIX(sgelsd_test, SGELSD_TEST)
void SGELSD_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, float *a, const lapack_int *lda,
                 float *b, const lapack_int *ldb, float *s, const float *rcond,
                 lapack_int *rank, float *work, const lapack_int *lwork,
                 lapack_int *iwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgelsd_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, s, *rcond, rank,
                                              work, *lwork, iwork);
        *info = lapacke_test_info("SGELSD", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGELSD", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgelsd)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, s, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_sgelsd_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, s, *rcond, rank, work,
                                          *lwork, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGELSD", ret);
}

/******************************************************************************
 * SGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGELSS_TEST LAPACK_GLOBAL_SUFFIX(sgelss_test, SGELSS_TEST)
void SGELSS_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, float *a, const lapack_int *lda,
                 float *b, const lapack_int *ldb, float *s, const float *rcond,
                 lapack_int *rank, float *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgelss_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, s, *rcond, rank,
                                              work, *lwork);
        *info = lapacke_test_info("SGELSS", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGELSS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgelss)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, s, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_sgelss_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, s, *rcond, rank, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGELSS", ret);
}

/******************************************************************************
 * SGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGELSY_TEST LAPACK_GLOBAL_SUFFIX(sgelsy_test, SGELSY_TEST)
void SGELSY_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, float *a, const lapack_int *lda,
                 float *b, const lapack_int *ldb, lapack_int *jpvt,
                 const float *rcond, lapack_int *rank, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgelsy_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, jpvt, *rcond,
                                              rank, work, *lwork);
        *info = lapacke_test_info("SGELSY", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGELSY", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgelsy)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, jpvt, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_sgelsy_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, jpvt, *rcond, rank, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGELSY", ret);
}

/******************************************************************************
 * SGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGETSLS_TEST LAPACK_GLOBAL_SUFFIX(sgetsls_test, SGETSLS_TEST)
void SGETSLS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                  const lapack_int *nrhs, float *a, const lapack_int *lda,
                  float *b, const lapack_int *ldb, float *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgetsls_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                               *nrhs, a, *lda, b, *ldb, work,
                                               *lwork);
        *info = lapacke_test_info("SGETSLS", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGETSLS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgetsls)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                      b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgetsls_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                           lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGETSLS", ret);
}

// ========================================================================== //
//                             SVD and bidiagonal                             //
// ========================================================================== //

/******************************************************************************
 * SBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )
 ******************************************************************************/
#define SBDSQR_TEST LAPACK_GLOBAL_SUFFIX(sbdsqr_test, SBDSQR_TEST)
void SBDSQR_TEST(const char *uplo, const lapack_int *n, const lapack_int *ncvt,
                 const lapack_int *nru, const lapack_int *ncc, float *d,
                 float *e, float *vt, const lapack_int *ldvt, float *u,
                 const lapack_int *ldu, float *c, const lapack_int *ldc,
                 float *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *vt_r = vt;
    lapack_int ldvt_r = *ldvt;
    float *u_r = u;
    lapack_int ldu_r = *ldu;
    float *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    vt_r = lapacke_test_sge_cm_to_rm(*n, *ncvt, vt, *ldvt, &ldvt_r);
    u_r = lapacke_test_sge_cm_to_rm(*nru, *n, u, *ldu, &ldu_r);
    c_r = lapacke_test_sge_cm_to_rm(*n, *ncc, c, *ldc, &ldc_r);
    if (vt_r == NULL || u_r == NULL || c_r == NULL) {
        LAPACKE_free(vt_r);
        LAPACKE_free(u_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("SBDSQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sbdsqr)(layout, *uplo, *n, *ncvt, *nru, *ncc, d, e,
                                     vt_r, ldvt_r, u_r, ldu_r, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_sbdsqr_work)(layout, *uplo, *n, *ncvt, *nru, *ncc,
                                          d, e, vt_r, ldvt_r, u_r, ldu_r, c_r,
                                          ldc_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *ncvt, vt_r, ldvt_r, vt, *ldvt);
    lapacke_test_sge_rm_to_cm(*nru, *n, u_r, ldu_r, u, *ldu);
    lapacke_test_sge_rm_to_cm(*n, *ncc, c_r, ldc_r, c, *ldc);
    LAPACKE_free(vt_r);
    LAPACKE_free(u_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("SBDSQR", ret);
}

/******************************************************************************
 * SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
 ******************************************************************************/
#define SGESVD_TEST LAPACK_GLOBAL_SUFFIX(sgesvd_test, SGESVD_TEST)
void SGESVD_TEST(const char *jobu, const char *jobvt, const lapack_int *m,
                 const lapack_int *n, float *a, const lapack_int *lda, float *s,
                 float *u, const lapack_int *ldu, float *vt,
                 const lapack_int *ldvt, float *work, const lapack_int *lwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN jobu_len, FORTRAN_STRLEN jobvt_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_sgesvd_work)(LAPACK_COL_MAJOR, *jobu, *jobvt,
                                              *m, *n, a, *lda, s, u, *ldu, vt,
                                              *ldvt, work, *lwork);
        *info = lapacke_test_info("SGESVD", ret);
        return;
    }

    float *a_r = a;
    lapack_int lda_r = *lda;
    float *u_r = u;
    lapack_int ldu_r = *ldu;
    float *vt_r = vt;
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
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    u_r = lapacke_test_sge_cm_to_rm(urows, ucols, u, *ldu, &ldu_r);
    vt_r = lapacke_test_sge_cm_to_rm(vtrows, *n, vt, *ldvt, &ldvt_r);
    if (a_r == NULL || u_r == NULL || vt_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(u_r);
        LAPACKE_free(vt_r);
        lapacke_test_report_alloc_failure("SGESVD", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgesvd)(layout, *jobu, *jobvt, *m, *n, a_r, lda_r,
                                     s, u_r, ldu_r, vt_r, ldvt_r, work);
#else
    ret = API_SUFFIX(LAPACKE_sgesvd_work)(layout, *jobu, *jobvt, *m, *n, a_r,
                                          lda_r, s, u_r, ldu_r, vt_r, ldvt_r,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(urows, ucols, u_r, ldu_r, u, *ldu);
    lapacke_test_sge_rm_to_cm(vtrows, *n, vt_r, ldvt_r, vt, *ldvt);
    LAPACKE_free(a_r);
    LAPACKE_free(u_r);
    LAPACKE_free(vt_r);
#endif
    *info = lapacke_test_info("SGESVD", ret);
}

// ========================================================================== //
//                             Auxiliary routines                             //
// ========================================================================== //

/******************************************************************************
 * REAL FUNCTION SLAMCH( CMACH )
 ******************************************************************************/
#define SLAMCH_TEST LAPACK_GLOBAL_SUFFIX(slamch_test, SLAMCH_TEST)
lapack_float_return SLAMCH_TEST(const char *cmach
#ifdef LAPACK_FORTRAN_STRLEN_END
                                ,
                                FORTRAN_STRLEN cmach_len
#endif
)
{
    lapack_float_return res = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_slamch)(*cmach);
#else
    res = API_SUFFIX(LAPACKE_slamch_work)(*cmach);
#endif

    return res;
}

/******************************************************************************
 * REAL FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )
 ******************************************************************************/
#define SLANGE_TEST LAPACK_GLOBAL_SUFFIX(slange_test, SLANGE_TEST)
lapack_float_return SLANGE_TEST(const char *norm, const lapack_int *m,
                                const lapack_int *n, const float *a,
                                const lapack_int *lda, float *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                                ,
                                FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_float_return res = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("SLANGE", &info);
        return (lapack_float_return)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_slange)(layout, *norm, *m, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_slange_work)(layout, *norm, *m, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * REAL FUNCTION SLANSY( NORM, UPLO, N, A, LDA, WORK )
 ******************************************************************************/
#define SLANSY_TEST LAPACK_GLOBAL_SUFFIX(slansy_test, SLANSY_TEST)
lapack_float_return SLANSY_TEST(const char *norm, const char *uplo,
                                const lapack_int *n, const float *a,
                                const lapack_int *lda, float *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                                ,
                                FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_float_return res = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_ssy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("SLANSY", &info);
        return (lapack_float_return)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_slansy)(layout, *norm, *uplo, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_slansy_work)(layout, *norm, *uplo, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * REAL FUNCTION SLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
 ******************************************************************************/
#define SLANTR_TEST LAPACK_GLOBAL_SUFFIX(slantr_test, SLANTR_TEST)
lapack_float_return SLANTR_TEST(const char *norm, const char *uplo,
                                const char *diag, const lapack_int *m,
                                const lapack_int *n, const float *a,
                                const lapack_int *lda, float *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                                ,
                                FORTRAN_STRLEN norm_len,
                                FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_float_return res = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("SLANTR", &info);
        return (lapack_float_return)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_slantr)(layout, *norm, *uplo, *diag, *m, *n, a_r,
                                     lda_r);
#else
    res = API_SUFFIX(LAPACKE_slantr_work)(layout, *norm, *uplo, *diag, *m, *n,
                                          a_r, lda_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * SLARNV( IDIST, ISEED, N, X )
 ******************************************************************************/
#define SLARNV_TEST LAPACK_GLOBAL_SUFFIX(slarnv_test, SLARNV_TEST)
void SLARNV_TEST(const lapack_int *idist, lapack_int *iseed,
                 const lapack_int *n, float *x)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_slarnv)(*idist, iseed, *n, x);
#else
    ret = API_SUFFIX(LAPACKE_slarnv_work)(*idist, iseed, *n, x);
#endif

    lapacke_test_info_unshifted("SLARNV", ret);
}

/******************************************************************************
 * SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
 ******************************************************************************/
#define SLASCL_TEST LAPACK_GLOBAL_SUFFIX(slascl_test, SLASCL_TEST)
void SLASCL_TEST(const char *type, const lapack_int *kl, const lapack_int *ku,
                 const float *cfrom, const float *cto, const lapack_int *m,
                 const lapack_int *n, float *a, const lapack_int *lda,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN type_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int arows = API_SUFFIX(LAPACKE_lsame)(*type, 'b')   ? *kl + 1
                             : API_SUFFIX(LAPACKE_lsame)(*type, 'q') ? *ku + 1
                             : API_SUFFIX(LAPACKE_lsame)(*type, 'z')
                                 ? 2 * *kl + *ku + 1
                                 : *m;
    a_r = lapacke_test_sge_cm_to_rm(arows, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SLASCL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_slascl)(layout, *type, *kl, *ku, *cfrom, *cto, *m,
                                     *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_slascl_work)(layout, *type, *kl, *ku, *cfrom, *cto,
                                          *m, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(arows, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SLASCL", ret);
}

/******************************************************************************
 * SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
 ******************************************************************************/
#define SLASET_TEST LAPACK_GLOBAL_SUFFIX(slaset_test, SLASET_TEST)
void SLASET_TEST(const char *uplo, const lapack_int *m, const lapack_int *n,
                 const float *alpha, const float *beta, float *a,
                 const lapack_int *lda
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("SLASET", &info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_slaset)(layout, *uplo, *m, *n, *alpha, *beta, a_r,
                                     lda_r);
#else
    ret = API_SUFFIX(LAPACKE_slaset_work)(layout, *uplo, *m, *n, *alpha, *beta,
                                          a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    lapacke_test_info("SLASET", ret);
}

/******************************************************************************
 * SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
 ******************************************************************************/
#define SLASWP_TEST LAPACK_GLOBAL_SUFFIX(slaswp_test, SLASWP_TEST)
void SLASWP_TEST(const lapack_int *n, float *a, const lapack_int *lda,
                 const lapack_int *k1, const lapack_int *k2,
                 const lapack_int *ipiv, const lapack_int *incx)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(
        lapacke_test_laswp_rows(*k1, *k2, ipiv, *incx), *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("SLASWP", &info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_slaswp)(layout, *n, a_r, lda_r, *k1, *k2, ipiv,
                                     *incx);
#else
    ret = API_SUFFIX(LAPACKE_slaswp_work)(layout, *n, a_r, lda_r, *k1, *k2,
                                          ipiv, *incx);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(lapacke_test_laswp_rows(*k1, *k2, ipiv, *incx),
                              *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    lapacke_test_info("SLASWP", ret);
}
