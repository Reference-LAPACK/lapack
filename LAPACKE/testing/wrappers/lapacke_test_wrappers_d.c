/******************************************************************************
 * LAPACKE test wrappers (double precision real)
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
 * DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DGECON_TEST LAPACK_GLOBAL_SUFFIX(dgecon_test, DGECON_TEST)
void DGECON_TEST(const char *norm, const lapack_int *n, const double *a,
                 const lapack_int *lda, const double *anorm, double *rcond,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGECON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgecon)(layout, *norm, *n, a_r, lda_r, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_dgecon_work)(layout, *norm, *n, a_r, lda_r, *anorm,
                                          rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("DGECON", ret);
}

/******************************************************************************
 * DGECXX( FACT, USESD, M, N, SESEL_ROWS, SEL_DESEL_COLS, KMAXFREE, ABSTOL,
 * RELTOL, A, LDA, K, MAXC2NRMK, RELMAXC2NRMK, FNRMK, IPIV, JPIV, TAU, C, LDC,
 * QRC, LDQRC, X, LDX, WORK, LWORK, IWORK, LIWORK, INFO )
 ******************************************************************************/
#define DGECXX_TEST LAPACK_GLOBAL_SUFFIX(dgecxx_test, DGECXX_TEST)
void DGECXX_TEST(const char *fact, const char *usesd, const lapack_int *m,
                 const lapack_int *n, const lapack_int *sesel_rows,
                 const lapack_int *sel_desel_cols, const lapack_int *kmaxfree,
                 const double *abstol, const double *reltol, double *a,
                 const lapack_int *lda, lapack_int *k, double *maxc2nrmk,
                 double *relmaxc2nrmk, double *fnrmk, lapack_int *ipiv,
                 lapack_int *jpiv, double *tau, double *c,
                 const lapack_int *ldc, double *qrc, const lapack_int *ldqrc,
                 double *x, const lapack_int *ldx, double *work,
                 const lapack_int *lwork, lapack_int *iwork,
                 const lapack_int *liwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN usesd_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgecxx_work)(
            LAPACK_COL_MAJOR, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
            (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a, *lda,
            k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c, *ldc, qrc,
            *ldqrc, x, *ldx, work, *lwork, iwork, *liwork);
        *info = lapacke_test_info("DGECXX", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *c_r = c;
    lapack_int ldc_r = *ldc;
    double *qrc_r = qrc;
    lapack_int ldqrc_r = *ldqrc;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    c_r = lapacke_test_dge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    qrc_r = lapacke_test_dge_cm_to_rm(*m, MIN(*m, *n), qrc, *ldqrc, &ldqrc_r);
    x_r = lapacke_test_dge_cm_to_rm(*m, *n, x, *ldx, &ldx_r);
    if (a_r == NULL || c_r == NULL || qrc_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(c_r);
        LAPACKE_free(qrc_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DGECXX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgecxx)(
        layout, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
        (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a_r, lda_r,
        k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c_r, ldc_r, qrc_r,
        ldqrc_r, x_r, ldx_r);
#else
    ret = API_SUFFIX(LAPACKE_dgecxx_work)(
        layout, *fact, *usesd, *m, *n, (lapack_int *)sesel_rows,
        (lapack_int *)sel_desel_cols, *kmaxfree, *abstol, *reltol, a_r, lda_r,
        k, maxc2nrmk, relmaxc2nrmk, fnrmk, ipiv, jpiv, tau, c_r, ldc_r, qrc_r,
        ldqrc_r, x_r, ldx_r, work, *lwork, iwork, *liwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    lapacke_test_dge_rm_to_cm(*m, MIN(*m, *n), qrc_r, ldqrc_r, qrc, *ldqrc);
    lapacke_test_dge_rm_to_cm(*m, *n, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(c_r);
    LAPACKE_free(qrc_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DGECXX", ret);
}

/******************************************************************************
 * DGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )
 ******************************************************************************/
#define DGEEQU_TEST LAPACK_GLOBAL_SUFFIX(dgeequ_test, DGEEQU_TEST)
void DGEEQU_TEST(const lapack_int *m, const lapack_int *n, const double *a,
                 const lapack_int *lda, double *r, double *c, double *rowcnd,
                 double *colcnd, double *amax, lapack_int *info)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGEEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeequ)(layout, *m, *n, a_r, lda_r, r, c, rowcnd,
                                     colcnd, amax);
#else
    ret = API_SUFFIX(LAPACKE_dgeequ_work)(layout, *m, *n, a_r, lda_r, r, c,
                                          rowcnd, colcnd, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("DGEEQU", ret);
}

/******************************************************************************
 * DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, IWORK, INFO )
 ******************************************************************************/
#define DGERFS_TEST LAPACK_GLOBAL_SUFFIX(dgerfs_test, DGERFS_TEST)
void DGERFS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, const double *af,
                 const lapack_int *ldaf, const lapack_int *ipiv,
                 const double *b, const lapack_int *ldb, double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    const double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    af_r = lapacke_test_dge_cm_to_rm(*n, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DGERFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgerfs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     af_r, ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dgerfs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DGERFS", ret);
}

/******************************************************************************
 * DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DGESV_TEST LAPACK_GLOBAL_SUFFIX(dgesv_test, DGESV_TEST)
void DGESV_TEST(const lapack_int *n, const lapack_int *nrhs, double *a,
                const lapack_int *lda, lapack_int *ipiv, double *b,
                const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGESV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgesv)(layout, *n, *nrhs, a_r, lda_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgesv_work)(layout, *n, *nrhs, a_r, lda_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGESV", ret);
}

/******************************************************************************
 * DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define DGESVX_TEST LAPACK_GLOBAL_SUFFIX(dgesvx_test, DGESVX_TEST)
void DGESVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *nrhs, double *a, const lapack_int *lda,
                 double *af, const lapack_int *ldaf, lapack_int *ipiv,
                 char *equed, double *r, double *c, double *b,
                 const lapack_int *ldb, double *x, const lapack_int *ldx,
                 double *rcond, double *ferr, double *berr, double *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
    double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    af_r = lapacke_test_dge_cm_to_rm(*n, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(af_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DGESVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgesvx)(
        layout, *fact, *trans, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, equed,
        r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work);
#else
    ret = API_SUFFIX(LAPACKE_dgesvx_work)(
        layout, *fact, *trans, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, equed,
        r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*n, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(af_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DGESVX", ret);
}

/******************************************************************************
 * DGETF2( M, N, A, LDA, IPIV, INFO )
 ******************************************************************************/
#define DGETF2_TEST LAPACK_GLOBAL_SUFFIX(dgetf2_test, DGETF2_TEST)
void DGETF2_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGETF2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgetf2)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dgetf2_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGETF2", ret);
}

/******************************************************************************
 * DGETRF( M, N, A, LDA, IPIV, INFO )
 ******************************************************************************/
#define DGETRF_TEST LAPACK_GLOBAL_SUFFIX(dgetrf_test, DGETRF_TEST)
void DGETRF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgetrf)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dgetrf_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGETRF", ret);
}

/******************************************************************************
 * DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGETRI_TEST LAPACK_GLOBAL_SUFFIX(dgetri_test, DGETRI_TEST)
void DGETRI_TEST(const lapack_int *n, double *a, const lapack_int *lda,
                 const lapack_int *ipiv, double *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgetri_work)(LAPACK_COL_MAJOR, *n, a, *lda,
                                              ipiv, work, *lwork);
        *info = lapacke_test_info("DGETRI", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgetri)(layout, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dgetri_work)(layout, *n, a_r, lda_r, ipiv, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGETRI", ret);
}

/******************************************************************************
 * DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DGETRS_TEST LAPACK_GLOBAL_SUFFIX(dgetrs_test, DGETRS_TEST)
void DGETRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, const lapack_int *ipiv,
                 double *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgetrs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgetrs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGETRS", ret);
}

// ========================================================================== //
//                         General band matrices (GB)                         //
// ========================================================================== //

/******************************************************************************
 * DGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DGBCON_TEST LAPACK_GLOBAL_SUFFIX(dgbcon_test, DGBCON_TEST)
void DGBCON_TEST(const char *norm, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const double *ab, const lapack_int *ldab,
                 const lapack_int *ipiv, const double *anorm, double *rcond,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("DGBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgbcon)(layout, *norm, *n, *kl, *ku, ab_r, ldab_r,
                                     ipiv, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_dgbcon_work)(layout, *norm, *n, *kl, *ku, ab_r,
                                          ldab_r, ipiv, *anorm, rcond, work,
                                          iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("DGBCON", ret);
}

/******************************************************************************
 * DGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFO )
 ******************************************************************************/
#define DGBEQU_TEST LAPACK_GLOBAL_SUFFIX(dgbequ_test, DGBEQU_TEST)
void DGBEQU_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const double *ab, const lapack_int *ldab,
                 double *r, double *c, double *rowcnd, double *colcnd,
                 double *amax, lapack_int *info)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dgb_cm_to_rm(*m, *n, *kl, *ku, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("DGBEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgbequ)(layout, *m, *n, *kl, *ku, ab_r, ldab_r, r,
                                     c, rowcnd, colcnd, amax);
#else
    ret = API_SUFFIX(LAPACKE_dgbequ_work)(layout, *m, *n, *kl, *ku, ab_r,
                                          ldab_r, r, c, rowcnd, colcnd, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("DGBEQU", ret);
}

/******************************************************************************
 * DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX,
 * FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define DGBRFS_TEST LAPACK_GLOBAL_SUFFIX(dgbrfs_test, DGBRFS_TEST)
void DGBRFS_TEST(const char *trans, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_int *nrhs, const double *ab,
                 const lapack_int *ldab, const double *afb,
                 const lapack_int *ldafb, const lapack_int *ipiv,
                 const double *b, const lapack_int *ldb, double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const double *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dgb_cm_to_rm(*n, *n, *kl, *ku, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_dgb_cm_to_rm(*n, *n, *kl, *kl + *ku, afb, *ldafb,
                                      &ldafb_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)afb_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DGBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgbrfs)(layout, *trans, *n, *kl, *ku, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, ipiv, b_r, ldb_r,
                                     x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dgbrfs_work)(
        layout, *trans, *n, *kl, *ku, *nrhs, ab_r, ldab_r, afb_r, ldafb_r, ipiv,
        b_r, ldb_r, x_r, ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)afb_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DGBRFS", ret);
}

/******************************************************************************
 * DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DGBSV_TEST LAPACK_GLOBAL_SUFFIX(dgbsv_test, DGBSV_TEST)
void DGBSV_TEST(const lapack_int *n, const lapack_int *kl, const lapack_int *ku,
                const lapack_int *nrhs, double *ab, const lapack_int *ldab,
                lapack_int *ipiv, double *b, const lapack_int *ldb,
                lapack_int *info)
{
    lapack_int ret = 0;
    double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGBSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgbsv)(layout, *n, *kl, *ku, *nrhs, ab_r, ldab_r,
                                    ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgbsv_work)(layout, *n, *kl, *ku, *nrhs, ab_r,
                                         ldab_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dgb_rm_to_cm(*n, *n, *kl, *kl + *ku, ab_r, ldab_r, ab, *ldab);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGBSV", ret);
}

/******************************************************************************
 * DGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, EQUED, R,
 * C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define DGBSVX_TEST LAPACK_GLOBAL_SUFFIX(dgbsvx_test, DGBSVX_TEST)
void DGBSVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *kl, const lapack_int *ku,
                 const lapack_int *nrhs, double *ab, const lapack_int *ldab,
                 double *afb, const lapack_int *ldafb, lapack_int *ipiv,
                 char *equed, double *r, double *c, double *b,
                 const lapack_int *ldb, double *x, const lapack_int *ldx,
                 double *rcond, double *ferr, double *berr, double *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    double *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dgb_cm_to_rm(*n, *n, *kl, *ku, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_dgb_cm_to_rm(*n, *n, *kl, *kl + *ku, afb, *ldafb,
                                      &ldafb_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(afb_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DGBSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgbsvx)(layout, *fact, *trans, *n, *kl, *ku, *nrhs,
                                     ab_r, ldab_r, afb_r, ldafb_r, ipiv, equed,
                                     r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr,
                                     berr, work);
#else
    ret = API_SUFFIX(LAPACKE_dgbsvx_work)(
        layout, *fact, *trans, *n, *kl, *ku, *nrhs, ab_r, ldab_r, afb_r,
        ldafb_r, ipiv, equed, r, c, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr,
        work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dgb_rm_to_cm(*n, *n, *kl, *ku, ab_r, ldab_r, ab, *ldab);
    lapacke_test_dgb_rm_to_cm(*n, *n, *kl, *kl + *ku, afb_r, ldafb_r, afb,
                              *ldafb);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(ab_r);
    LAPACKE_free(afb_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DGBSVX", ret);
}

/******************************************************************************
 * DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
 ******************************************************************************/
#define DGBTRF_TEST LAPACK_GLOBAL_SUFFIX(dgbtrf_test, DGBTRF_TEST)
void DGBTRF_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, double *ab, const lapack_int *ldab,
                 lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dgb_cm_to_rm(*m, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("DGBTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgbtrf)(layout, *m, *n, *kl, *ku, ab_r, ldab_r,
                                     ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dgbtrf_work)(layout, *m, *n, *kl, *ku, ab_r,
                                          ldab_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dgb_rm_to_cm(*m, *n, *kl, *kl + *ku, ab_r, ldab_r, ab, *ldab);
    LAPACKE_free(ab_r);
#endif
    *info = lapacke_test_info("DGBTRF", ret);
}

/******************************************************************************
 * DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DGBTRS_TEST LAPACK_GLOBAL_SUFFIX(dgbtrs_test, DGBTRS_TEST)
void DGBTRS_TEST(const char *trans, const lapack_int *n, const lapack_int *kl,
                 const lapack_int *ku, const lapack_int *nrhs, const double *ab,
                 const lapack_int *ldab, const lapack_int *ipiv, double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dgb_cm_to_rm(*n, *n, *kl, *kl + *ku, ab, *ldab,
                                     &ldab_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgbtrs)(layout, *trans, *n, *kl, *ku, *nrhs, ab_r,
                                     ldab_r, ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgbtrs_work)(layout, *trans, *n, *kl, *ku, *nrhs,
                                          ab_r, ldab_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGBTRS", ret);
}

// ========================================================================== //
//                      General tridiagonal matrices (GT)                     //
// ========================================================================== //

/******************************************************************************
 * DGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DGTCON_TEST LAPACK_GLOBAL_SUFFIX(dgtcon_test, DGTCON_TEST)
void DGTCON_TEST(const char *norm, const lapack_int *n, const double *dl,
                 const double *d, const double *du, const double *du2,
                 const lapack_int *ipiv, const double *anorm, double *rcond,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len
#endif
)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgtcon)(*norm, *n, dl, d, du, du2, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_dgtcon_work)(*norm, *n, dl, d, du, du2, ipiv,
                                          *anorm, rcond, work, iwork);
#endif

    *info = lapacke_test_info_unshifted("DGTCON", ret);
}

/******************************************************************************
 * DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX,
 * FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define DGTRFS_TEST LAPACK_GLOBAL_SUFFIX(dgtrfs_test, DGTRFS_TEST)
void DGTRFS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const double *dl, const double *d, const double *du,
                 const double *dlf, const double *df, const double *duf,
                 const double *du2, const lapack_int *ipiv, const double *b,
                 const lapack_int *ldb, double *x, const lapack_int *ldx,
                 double *ferr, double *berr, double *work, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DGTRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgtrfs)(layout, *trans, *n, *nrhs, dl, d, du, dlf,
                                     df, duf, du2, ipiv, b_r, ldb_r, x_r, ldx_r,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dgtrfs_work)(layout, *trans, *n, *nrhs, dl, d, du,
                                          dlf, df, duf, du2, ipiv, b_r, ldb_r,
                                          x_r, ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DGTRFS", ret);
}

/******************************************************************************
 * DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
 ******************************************************************************/
#define DGTSV_TEST LAPACK_GLOBAL_SUFFIX(dgtsv_test, DGTSV_TEST)
void DGTSV_TEST(const lapack_int *n, const lapack_int *nrhs, double *dl,
                double *d, double *du, double *b, const lapack_int *ldb,
                lapack_int *info)
{
    lapack_int ret = 0;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("DGTSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgtsv)(layout, *n, *nrhs, dl, d, du, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgtsv_work)(layout, *n, *nrhs, dl, d, du, b_r,
                                         ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGTSV", ret);
}

/******************************************************************************
 * DGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define DGTSVX_TEST LAPACK_GLOBAL_SUFFIX(dgtsvx_test, DGTSVX_TEST)
void DGTSVX_TEST(const char *fact, const char *trans, const lapack_int *n,
                 const lapack_int *nrhs, const double *dl, const double *d,
                 const double *du, double *dlf, double *df, double *duf,
                 double *du2, lapack_int *ipiv, const double *b,
                 const lapack_int *ldb, double *x, const lapack_int *ldx,
                 double *rcond, double *ferr, double *berr, double *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DGTSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgtsvx)(layout, *fact, *trans, *n, *nrhs, dl, d,
                                     du, dlf, df, duf, du2, ipiv, b_r, ldb_r,
                                     x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dgtsvx_work)(
        layout, *fact, *trans, *n, *nrhs, dl, d, du, dlf, df, duf, du2, ipiv,
        b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DGTSVX", ret);
}

/******************************************************************************
 * DGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
 ******************************************************************************/
#define DGTTRF_TEST LAPACK_GLOBAL_SUFFIX(dgttrf_test, DGTTRF_TEST)
void DGTTRF_TEST(const lapack_int *n, double *dl, double *d, double *du,
                 double *du2, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgttrf)(*n, dl, d, du, du2, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dgttrf_work)(*n, dl, d, du, du2, ipiv);
#endif

    *info = lapacke_test_info_unshifted("DGTTRF", ret);
}

/******************************************************************************
 * DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DGTTRS_TEST LAPACK_GLOBAL_SUFFIX(dgttrs_test, DGTTRS_TEST)
void DGTTRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const double *dl, const double *d, const double *du,
                 const double *du2, const lapack_int *ipiv, double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("DGTTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgttrs)(layout, *trans, *n, *nrhs, dl, d, du, du2,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgttrs_work)(layout, *trans, *n, *nrhs, dl, d, du,
                                          du2, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGTTRS", ret);
}

// ========================================================================== //
//             Symmetric/Hermitian positive definite matrices (PO)            //
// ========================================================================== //

/******************************************************************************
 * DPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DPOCON_TEST LAPACK_GLOBAL_SUFFIX(dpocon_test, DPOCON_TEST)
void DPOCON_TEST(const char *uplo, const lapack_int *n, const double *a,
                 const lapack_int *lda, const double *anorm, double *rcond,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DPOCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpocon)(layout, *uplo, *n, a_r, lda_r, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_dpocon_work)(layout, *uplo, *n, a_r, lda_r, *anorm,
                                          rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("DPOCON", ret);
}

/******************************************************************************
 * DPOEQU( N, A, LDA, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define DPOEQU_TEST LAPACK_GLOBAL_SUFFIX(dpoequ_test, DPOEQU_TEST)
void DPOEQU_TEST(const lapack_int *n, const double *a, const lapack_int *lda,
                 double *s, double *scond, double *amax, lapack_int *info)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DPOEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpoequ)(layout, *n, a_r, lda_r, s, scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_dpoequ_work)(layout, *n, a_r, lda_r, s, scond,
                                          amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("DPOEQU", ret);
}

/******************************************************************************
 * DPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX, FERR, BERR, WORK,
 * IWORK, INFO )
 ******************************************************************************/
#define DPORFS_TEST LAPACK_GLOBAL_SUFFIX(dporfs_test, DPORFS_TEST)
void DPORFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, const double *af,
                 const lapack_int *ldaf, const double *b, const lapack_int *ldb,
                 double *x, const lapack_int *ldx, double *ferr, double *berr,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    const double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DPORFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dporfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_dporfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, b_r, ldb_r, x_r, ldx_r,
                                          ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DPORFS", ret);
}

/******************************************************************************
 * DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define DPOSV_TEST LAPACK_GLOBAL_SUFFIX(dposv_test, DPOSV_TEST)
void DPOSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                double *a, const lapack_int *lda, double *b,
                const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DPOSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dposv)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dposv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DPOSV", ret);
}

/******************************************************************************
 * DPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S, B, LDB, X, LDX,
 * RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define DPOSVX_TEST LAPACK_GLOBAL_SUFFIX(dposvx_test, DPOSVX_TEST)
void DPOSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, double *a, const lapack_int *lda,
                 double *af, const lapack_int *ldaf, char *equed, double *s,
                 double *b, const lapack_int *ldb, double *x,
                 const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
    double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(af_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DPOSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dposvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, equed, s, b_r, ldb_r,
                                     x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dposvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, equed, s,
        b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_dpo_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(a_r);
    LAPACKE_free(af_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DPOSVX", ret);
}

/******************************************************************************
 * DPOTRF( UPLO, N, A, LDA, INFO )
 ******************************************************************************/
#define DPOTRF_TEST LAPACK_GLOBAL_SUFFIX(dpotrf_test, DPOTRF_TEST)
void DPOTRF_TEST(const char *uplo, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DPOTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpotrf)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_dpotrf_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DPOTRF", ret);
}

/******************************************************************************
 * DPOTRI( UPLO, N, A, LDA, INFO )
 ******************************************************************************/
#define DPOTRI_TEST LAPACK_GLOBAL_SUFFIX(dpotri_test, DPOTRI_TEST)
void DPOTRI_TEST(const char *uplo, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DPOTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpotri)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_dpotri_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DPOTRI", ret);
}

/******************************************************************************
 * DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define DPOTRS_TEST LAPACK_GLOBAL_SUFFIX(dpotrs_test, DPOTRS_TEST)
void DPOTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DPOTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpotrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dpotrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DPOTRS", ret);
}

/******************************************************************************
 * DPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
 ******************************************************************************/
#define DPSTRF_TEST LAPACK_GLOBAL_SUFFIX(dpstrf_test, DPSTRF_TEST)
void DPSTRF_TEST(const char *uplo, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *piv, lapack_int *rank,
                 const double *tol, double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DPSTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpstrf)(layout, *uplo, *n, a_r, lda_r, piv, rank,
                                     *tol);
#else
    ret = API_SUFFIX(LAPACKE_dpstrf_work)(layout, *uplo, *n, a_r, lda_r, piv,
                                          rank, *tol, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DPSTRF", ret);
}

// ========================================================================== //
//                   Packed positive definite matrices (PP)                   //
// ========================================================================== //

/******************************************************************************
 * DPPCON( UPLO, N, AP, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DPPCON_TEST LAPACK_GLOBAL_SUFFIX(dppcon_test, DPPCON_TEST)
void DPPCON_TEST(const char *uplo, const lapack_int *n, const double *ap,
                 const double *anorm, double *rcond, double *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("DPPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dppcon)(layout, *uplo, *n, ap_r, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_dppcon_work)(layout, *uplo, *n, ap_r, *anorm,
                                          rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("DPPCON", ret);
}

/******************************************************************************
 * DPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define DPPEQU_TEST LAPACK_GLOBAL_SUFFIX(dppequ_test, DPPEQU_TEST)
void DPPEQU_TEST(const char *uplo, const lapack_int *n, const double *ap,
                 double *s, double *scond, double *amax, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("DPPEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dppequ)(layout, *uplo, *n, ap_r, s, scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_dppequ_work)(layout, *uplo, *n, ap_r, s, scond,
                                          amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("DPPEQU", ret);
}

/******************************************************************************
 * DPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO
 * )
 ******************************************************************************/
#define DPPRFS_TEST LAPACK_GLOBAL_SUFFIX(dpprfs_test, DPPRFS_TEST)
void DPPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *ap, const double *afp, const double *b,
                 const lapack_int *ldb, double *x, const lapack_int *ldx,
                 double *ferr, double *berr, double *work, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
    const double *ap_r = ap;
    const double *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("DPPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r, b_r,
                                     ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dpprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          b_r, ldb_r, x_r, ldx_r, ferr, berr,
                                          work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("DPPRFS", ret);
}

/******************************************************************************
 * DPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
 ******************************************************************************/
#define DPPSV_TEST LAPACK_GLOBAL_SUFFIX(dppsv_test, DPPSV_TEST)
void DPPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                double *ap, double *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("DPPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dppsv)(layout, *uplo, *n, *nrhs, ap_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dppsv_work)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                         ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_dpp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("DPPSV", ret);
}

/******************************************************************************
 * DPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define DPPSVX_TEST LAPACK_GLOBAL_SUFFIX(dppsvx_test, DPPSVX_TEST)
void DPPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, double *ap, double *afp, char *equed,
                 double *s, double *b, const lapack_int *ldb, double *x,
                 const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *ap_r = ap;
    double *afp_r = afp;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DPPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dppsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, equed, s, b_r, ldb_r, x_r, ldx_r,
                                     rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dppsvx_work)(
        layout, *fact, *uplo, *n, *nrhs, ap_r, afp_r, equed, s, b_r, ldb_r, x_r,
        ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_dpp_rm_to_cm(*uplo, *n, ap_r, ap);
    lapacke_test_dpp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DPPSVX", ret);
}

/******************************************************************************
 * DPPTRF( UPLO, N, AP, INFO )
 ******************************************************************************/
#define DPPTRF_TEST LAPACK_GLOBAL_SUFFIX(dpptrf_test, DPPTRF_TEST)
void DPPTRF_TEST(const char *uplo, const lapack_int *n, double *ap,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("DPPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpptrf)(layout, *uplo, *n, ap_r);
#else
    ret = API_SUFFIX(LAPACKE_dpptrf_work)(layout, *uplo, *n, ap_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("DPPTRF", ret);
}

/******************************************************************************
 * DPPTRI( UPLO, N, AP, INFO )
 ******************************************************************************/
#define DPPTRI_TEST LAPACK_GLOBAL_SUFFIX(dpptri_test, DPPTRI_TEST)
void DPPTRI_TEST(const char *uplo, const lapack_int *n, double *ap,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("DPPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpptri)(layout, *uplo, *n, ap_r);
#else
    ret = API_SUFFIX(LAPACKE_dpptri_work)(layout, *uplo, *n, ap_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("DPPTRI", ret);
}

/******************************************************************************
 * DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
 ******************************************************************************/
#define DPPTRS_TEST LAPACK_GLOBAL_SUFFIX(dpptrs_test, DPPTRS_TEST)
void DPPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *ap, double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
    const double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_dpp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("DPPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpptrs)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dpptrs_work)(layout, *uplo, *n, *nrhs, ap_r, b_r,
                                          ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("DPPTRS", ret);
}

// ========================================================================== //
//                    Positive definite band matrices (PB)                    //
// ========================================================================== //

/******************************************************************************
 * DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DPBCON_TEST LAPACK_GLOBAL_SUFFIX(dpbcon_test, DPBCON_TEST)
void DPBCON_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const double *ab, const lapack_int *ldab, const double *anorm,
                 double *rcond, double *work, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("DPBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpbcon)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_dpbcon_work)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                          *anorm, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("DPBCON", ret);
}

/******************************************************************************
 * DPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )
 ******************************************************************************/
#define DPBEQU_TEST LAPACK_GLOBAL_SUFFIX(dpbequ_test, DPBEQU_TEST)
void DPBEQU_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const double *ab, const lapack_int *ldab, double *s,
                 double *scond, double *amax, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("DPBEQU", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpbequ)(layout, *uplo, *n, *kd, ab_r, ldab_r, s,
                                     scond, amax);
#else
    ret = API_SUFFIX(LAPACKE_dpbequ_work)(layout, *uplo, *n, *kd, ab_r, ldab_r,
                                          s, scond, amax);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("DPBEQU", ret);
}

/******************************************************************************
 * DPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR,
 * WORK, IWORK, INFO )
 ******************************************************************************/
#define DPBRFS_TEST LAPACK_GLOBAL_SUFFIX(dpbrfs_test, DPBRFS_TEST)
void DPBRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const double *ab,
                 const lapack_int *ldab, const double *afb,
                 const lapack_int *ldafb, const double *b,
                 const lapack_int *ldb, double *x, const lapack_int *ldx,
                 double *ferr, double *berr, double *work, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const double *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_dpb_cm_to_rm(*uplo, *n, *kd, afb, *ldafb, &ldafb_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)afb_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DPBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpbrfs)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, b_r, ldb_r, x_r,
                                     ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dpbrfs_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                          ldab_r, afb_r, ldafb_r, b_r, ldb_r,
                                          x_r, ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)afb_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DPBRFS", ret);
}

/******************************************************************************
 * DPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define DPBSV_TEST LAPACK_GLOBAL_SUFFIX(dpbsv_test, DPBSV_TEST)
void DPBSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                const lapack_int *nrhs, double *ab, const lapack_int *ldab,
                double *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DPBSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpbsv)(layout, *uplo, *n, *kd, *nrhs, ab_r, ldab_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dpbsv_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                         ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DPBSV", ret);
}

/******************************************************************************
 * DPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X,
 * LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define DPBSVX_TEST LAPACK_GLOBAL_SUFFIX(dpbsvx_test, DPBSVX_TEST)
void DPBSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *kd, const lapack_int *nrhs, double *ab,
                 const lapack_int *ldab, double *afb, const lapack_int *ldafb,
                 char *equed, double *s, double *b, const lapack_int *ldb,
                 double *x, const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN equed_len
#endif
)
{
    lapack_int ret = 0;
    double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    double *afb_r = afb;
    lapack_int ldafb_r = *ldafb;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    afb_r = lapacke_test_dpb_cm_to_rm(*uplo, *n, *kd, afb, *ldafb, &ldafb_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || afb_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free(ab_r);
        LAPACKE_free(afb_r);
        LAPACKE_free(b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DPBSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpbsvx)(layout, *fact, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, afb_r, ldafb_r, equed, s, b_r,
                                     ldb_r, x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dpbsvx_work)(
        layout, *fact, *uplo, *n, *kd, *nrhs, ab_r, ldab_r, afb_r, ldafb_r,
        equed, s, b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    lapacke_test_dpb_rm_to_cm(*uplo, *n, *kd, afb_r, ldafb_r, afb, *ldafb);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free(ab_r);
    LAPACKE_free(afb_r);
    LAPACKE_free(b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DPBSVX", ret);
}

/******************************************************************************
 * DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
 ******************************************************************************/
#define DPBTRF_TEST LAPACK_GLOBAL_SUFFIX(dpbtrf_test, DPBTRF_TEST)
void DPBTRF_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 double *ab, const lapack_int *ldab, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("DPBTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpbtrf)(layout, *uplo, *n, *kd, ab_r, ldab_r);
#else
    ret = API_SUFFIX(LAPACKE_dpbtrf_work)(layout, *uplo, *n, *kd, ab_r, ldab_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpb_rm_to_cm(*uplo, *n, *kd, ab_r, ldab_r, ab, *ldab);
    LAPACKE_free(ab_r);
#endif
    *info = lapacke_test_info("DPBTRF", ret);
}

/******************************************************************************
 * DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define DPBTRS_TEST LAPACK_GLOBAL_SUFFIX(dpbtrs_test, DPBTRS_TEST)
void DPBTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const double *ab,
                 const lapack_int *ldab, double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dpb_cm_to_rm(*uplo, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DPBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpbtrs)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                     ldab_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dpbtrs_work)(layout, *uplo, *n, *kd, *nrhs, ab_r,
                                          ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DPBTRS", ret);
}

// ========================================================================== //
//                 Positive definite tridiagonal matrices (PT)                //
// ========================================================================== //

/******************************************************************************
 * DPTCON( N, D, E, ANORM, RCOND, WORK, INFO )
 ******************************************************************************/
#define DPTCON_TEST LAPACK_GLOBAL_SUFFIX(dptcon_test, DPTCON_TEST)
void DPTCON_TEST(const lapack_int *n, const double *d, const double *e,
                 const double *anorm, double *rcond, double *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dptcon)(*n, d, e, *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_dptcon_work)(*n, d, e, *anorm, rcond, work);
#endif

    *info = lapacke_test_info_unshifted("DPTCON", ret);
}

/******************************************************************************
 * DPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, INFO )
 ******************************************************************************/
#define DPTRFS_TEST LAPACK_GLOBAL_SUFFIX(dptrfs_test, DPTRFS_TEST)
void DPTRFS_TEST(const lapack_int *n, const lapack_int *nrhs, const double *d,
                 const double *e, const double *df, const double *ef,
                 const double *b, const lapack_int *ldb, double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 double *work, lapack_int *info)
{
    lapack_int ret = 0;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DPTRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dptrfs)(layout, *n, *nrhs, d, e, df, ef, b_r,
                                     ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dptrfs_work)(layout, *n, *nrhs, d, e, df, ef, b_r,
                                          ldb_r, x_r, ldx_r, ferr, berr, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DPTRFS", ret);
}

/******************************************************************************
 * DPTSV( N, NRHS, D, E, B, LDB, INFO )
 ******************************************************************************/
#define DPTSV_TEST LAPACK_GLOBAL_SUFFIX(dptsv_test, DPTSV_TEST)
void DPTSV_TEST(const lapack_int *n, const lapack_int *nrhs, double *d,
                double *e, double *b, const lapack_int *ldb, lapack_int *info)
{
    lapack_int ret = 0;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("DPTSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dptsv)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dptsv_work)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DPTSV", ret);
}

/******************************************************************************
 * DPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK,
 * INFO )
 ******************************************************************************/
#define DPTSVX_TEST LAPACK_GLOBAL_SUFFIX(dptsvx_test, DPTSVX_TEST)
void DPTSVX_TEST(const char *fact, const lapack_int *n, const lapack_int *nrhs,
                 const double *d, const double *e, double *df, double *ef,
                 const double *b, const lapack_int *ldb, double *x,
                 const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len
#endif
)
{
    lapack_int ret = 0;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DPTSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dptsvx)(layout, *fact, *n, *nrhs, d, e, df, ef,
                                     b_r, ldb_r, x_r, ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dptsvx_work)(layout, *fact, *n, *nrhs, d, e, df,
                                          ef, b_r, ldb_r, x_r, ldx_r, rcond,
                                          ferr, berr, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DPTSVX", ret);
}

/******************************************************************************
 * DPTTRF( N, D, E, INFO )
 ******************************************************************************/
#define DPTTRF_TEST LAPACK_GLOBAL_SUFFIX(dpttrf_test, DPTTRF_TEST)
void DPTTRF_TEST(const lapack_int *n, double *d, double *e, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpttrf)(*n, d, e);
#else
    ret = API_SUFFIX(LAPACKE_dpttrf_work)(*n, d, e);
#endif

    *info = lapacke_test_info_unshifted("DPTTRF", ret);
}

/******************************************************************************
 * DPTTRS( N, NRHS, D, E, B, LDB, INFO )
 ******************************************************************************/
#define DPTTRS_TEST LAPACK_GLOBAL_SUFFIX(dpttrs_test, DPTTRS_TEST)
void DPTTRS_TEST(const lapack_int *n, const lapack_int *nrhs, const double *d,
                 const double *e, double *b, const lapack_int *ldb,
                 lapack_int *info)
{
    lapack_int ret = 0;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (b_r == NULL) {
        lapacke_test_report_alloc_failure("DPTTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpttrs)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dpttrs_work)(layout, *n, *nrhs, d, e, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DPTTRS", ret);
}

// ========================================================================== //
//                     Symmetric indefinite matrices (SY)                     //
// ========================================================================== //

/******************************************************************************
 * DSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DSYCON_TEST LAPACK_GLOBAL_SUFFIX(dsycon_test, DSYCON_TEST)
void DSYCON_TEST(const char *uplo, const lapack_int *n, const double *a,
                 const lapack_int *lda, const lapack_int *ipiv,
                 const double *anorm, double *rcond, double *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsycon)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                     *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_dsycon_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          *anorm, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("DSYCON", ret);
}

/******************************************************************************
 * DSYCON_3( UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DSYCON_3_TEST LAPACK_GLOBAL_SUFFIX(dsycon_3_test, DSYCON_3_TEST)
void DSYCON_3_TEST(const char *uplo, const lapack_int *n, const double *a,
                   const lapack_int *lda, const double *e,
                   const lapack_int *ipiv, const double *anorm, double *rcond,
                   double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYCON_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsycon_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv,
                                       *anorm, rcond);
#else
    ret = API_SUFFIX(LAPACKE_dsycon_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, *anorm, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("DSYCON_3", ret);
}

/******************************************************************************
 * DSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
 * WORK, IWORK, INFO )
 ******************************************************************************/
#define DSYRFS_TEST LAPACK_GLOBAL_SUFFIX(dsyrfs_test, DSYRFS_TEST)
void DSYRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, const double *af,
                 const lapack_int *ldaf, const lapack_int *ipiv,
                 const double *b, const lapack_int *ldb, double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    const double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DSYRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsyrfs)(layout, *uplo, *n, *nrhs, a_r, lda_r, af_r,
                                     ldaf_r, ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_dsyrfs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DSYRFS", ret);
}

/******************************************************************************
 * DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYSV_TEST LAPACK_GLOBAL_SUFFIX(dsysv_test, DSYSV_TEST)
void DSYSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                double *a, const lapack_int *lda, lapack_int *ipiv, double *b,
                const lapack_int *ldb, double *work, const lapack_int *lwork,
                lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsysv_work)(LAPACK_COL_MAJOR, *uplo, *n, *nrhs,
                                             a, *lda, ipiv, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("DSYSV", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DSYSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsysv)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsysv_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DSYSV", ret);
}

/******************************************************************************
 * DSYSV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYSV_AA_TEST LAPACK_GLOBAL_SUFFIX(dsysv_aa_test, DSYSV_AA_TEST)
void DSYSV_AA_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, double *a, const lapack_int *lda,
                   lapack_int *ipiv, double *b, const lapack_int *ldb,
                   double *work, const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsysv_aa_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, ipiv, b, *ldb,
                                                work, *lwork);
        *info = lapacke_test_info("DSYSV_AA", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DSYSV_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsysv_aa)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsysv_aa_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DSYSV_AA", ret);
}

/******************************************************************************
 * DSYSV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYSV_RK_TEST LAPACK_GLOBAL_SUFFIX(dsysv_rk_test, DSYSV_RK_TEST)
void DSYSV_RK_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, double *a, const lapack_int *lda,
                   double *e, lapack_int *ipiv, double *b,
                   const lapack_int *ldb, double *work, const lapack_int *lwork,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsysv_rk_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                *nrhs, a, *lda, e, ipiv, b,
                                                *ldb, work, *lwork);
        *info = lapacke_test_info("DSYSV_RK", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DSYSV_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsysv_rk)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsysv_rk_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r, work,
                                            *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DSYSV_RK", ret);
}

/******************************************************************************
 * DSYSV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYSV_ROOK_TEST LAPACK_GLOBAL_SUFFIX(dsysv_rook_test, DSYSV_ROOK_TEST)
void DSYSV_ROOK_TEST(const char *uplo, const lapack_int *n,
                     const lapack_int *nrhs, double *a, const lapack_int *lda,
                     lapack_int *ipiv, double *b, const lapack_int *ldb,
                     double *work, const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                     ,
                     FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsysv_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                  *nrhs, a, *lda, ipiv, b, *ldb,
                                                  work, *lwork);
        *info = lapacke_test_info("DSYSV_ROOK", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DSYSV_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsysv_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                         ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsysv_rook_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DSYSV_ROOK", ret);
}

/******************************************************************************
 * DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND,
 * FERR, BERR, WORK, LWORK, IWORK, INFO )
 ******************************************************************************/
#define DSYSVX_TEST LAPACK_GLOBAL_SUFFIX(dsysvx_test, DSYSVX_TEST)
void DSYSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const double *a, const lapack_int *lda,
                 double *af, const lapack_int *ldaf, lapack_int *ipiv,
                 const double *b, const lapack_int *ldb, double *x,
                 const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, double *work, const lapack_int *lwork,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsysvx_work)(
            LAPACK_COL_MAJOR, *fact, *uplo, *n, *nrhs, a, *lda, af, *ldaf, ipiv,
            b, *ldb, x, *ldx, rcond, ferr, berr, work, *lwork, iwork);
        *info = lapacke_test_info("DSYSVX", ret);
        return;
    }

    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *af_r = af;
    lapack_int ldaf_r = *ldaf;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    af_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, af, *ldaf, &ldaf_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || af_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(af_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DSYSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsysvx)(layout, *fact, *uplo, *n, *nrhs, a_r,
                                     lda_r, af_r, ldaf_r, ipiv, b_r, ldb_r, x_r,
                                     ldx_r, rcond, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dsysvx_work)(
        layout, *fact, *uplo, *n, *nrhs, a_r, lda_r, af_r, ldaf_r, ipiv, b_r,
        ldb_r, x_r, ldx_r, rcond, ferr, berr, work, *lwork, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, af_r, ldaf_r, af, *ldaf);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(af_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DSYSVX", ret);
}

/******************************************************************************
 * DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYTRF_TEST LAPACK_GLOBAL_SUFFIX(dsytrf_test, DSYTRF_TEST)
void DSYTRF_TEST(const char *uplo, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *ipiv, double *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsytrf_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                              *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("DSYTRF", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytrf)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dsytrf_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DSYTRF", ret);
}

/******************************************************************************
 * DSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYTRF_AA_TEST LAPACK_GLOBAL_SUFFIX(dsytrf_aa_test, DSYTRF_AA_TEST)
void DSYTRF_AA_TEST(const char *uplo, const lapack_int *n, double *a,
                    const lapack_int *lda, lapack_int *ipiv, double *work,
                    const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsytrf_aa_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("DSYTRF_AA", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYTRF_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytrf_aa)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dsytrf_aa_work)(layout, *uplo, *n, a_r, lda_r,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DSYTRF_AA", ret);
}

/******************************************************************************
 * DSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYTRF_RK_TEST LAPACK_GLOBAL_SUFFIX(dsytrf_rk_test, DSYTRF_RK_TEST)
void DSYTRF_RK_TEST(const char *uplo, const lapack_int *n, double *a,
                    const lapack_int *lda, double *e, lapack_int *ipiv,
                    double *work, const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsytrf_rk_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                 *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("DSYTRF_RK", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYTRF_RK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytrf_rk)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dsytrf_rk_work)(layout, *uplo, *n, a_r, lda_r, e,
                                             ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DSYTRF_RK", ret);
}

/******************************************************************************
 * DSYTRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYTRF_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(dsytrf_rook_test, DSYTRF_ROOK_TEST)
void DSYTRF_ROOK_TEST(const char *uplo, const lapack_int *n, double *a,
                      const lapack_int *lda, lapack_int *ipiv, double *work,
                      const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsytrf_rook_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                   a, *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("DSYTRF_ROOK", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYTRF_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytrf_rook)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dsytrf_rook_work)(layout, *uplo, *n, a_r, lda_r,
                                               ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DSYTRF_ROOK", ret);
}

/******************************************************************************
 * DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
 ******************************************************************************/
#define DSYTRI_TEST LAPACK_GLOBAL_SUFFIX(dsytri_test, DSYTRI_TEST)
void DSYTRI_TEST(const char *uplo, const lapack_int *n, double *a,
                 const lapack_int *lda, const lapack_int *ipiv, double *work,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytri)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dsytri_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DSYTRI", ret);
}

/******************************************************************************
 * DSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYTRI2_TEST LAPACK_GLOBAL_SUFFIX(dsytri2_test, DSYTRI2_TEST)
void DSYTRI2_TEST(const char *uplo, const lapack_int *n, double *a,
                  const lapack_int *lda, const lapack_int *ipiv, double *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsytri2_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                               *lda, ipiv, work, *lwork);
        *info = lapacke_test_info("DSYTRI2", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYTRI2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytri2)(layout, *uplo, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dsytri2_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DSYTRI2", ret);
}

/******************************************************************************
 * DSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )
 ******************************************************************************/
#define DSYTRI2X_TEST LAPACK_GLOBAL_SUFFIX(dsytri2x_test, DSYTRI2X_TEST)
void DSYTRI2X_TEST(const char *uplo, const lapack_int *n, double *a,
                   const lapack_int *lda, const lapack_int *ipiv, double *work,
                   const lapack_int *nb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYTRI2X", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytri2x)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                       *nb);
#else
    ret = API_SUFFIX(LAPACKE_dsytri2x_work)(layout, *uplo, *n, a_r, lda_r, ipiv,
                                            work, *nb);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DSYTRI2X", ret);
}

/******************************************************************************
 * DSYTRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYTRI_3_TEST LAPACK_GLOBAL_SUFFIX(dsytri_3_test, DSYTRI_3_TEST)
void DSYTRI_3_TEST(const char *uplo, const lapack_int *n, double *a,
                   const lapack_int *lda, const double *e,
                   const lapack_int *ipiv, double *work,
                   const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsytri_3_work)(LAPACK_COL_MAJOR, *uplo, *n, a,
                                                *lda, e, ipiv, work, *lwork);
        *info = lapacke_test_info("DSYTRI_3", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DSYTRI_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytri_3)(layout, *uplo, *n, a_r, lda_r, e, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dsytri_3_work)(layout, *uplo, *n, a_r, lda_r, e,
                                            ipiv, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsy_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DSYTRI_3", ret);
}

/******************************************************************************
 * DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DSYTRS_TEST LAPACK_GLOBAL_SUFFIX(dsytrs_test, DSYTRS_TEST)
void DSYTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, const lapack_int *ipiv,
                 double *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DSYTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv,
                                     b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsytrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DSYTRS", ret);
}

/******************************************************************************
 * DSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )
 ******************************************************************************/
#define DSYTRS2_TEST LAPACK_GLOBAL_SUFFIX(dsytrs2_test, DSYTRS2_TEST)
void DSYTRS2_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                  const double *a, const lapack_int *lda,
                  const lapack_int *ipiv, double *b, const lapack_int *ldb,
                  double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DSYTRS2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytrs2)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                      ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsytrs2_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                           ipiv, b_r, ldb_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DSYTRS2", ret);
}

/******************************************************************************
 * DSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DSYTRS_3_TEST LAPACK_GLOBAL_SUFFIX(dsytrs_3_test, DSYTRS_3_TEST)
void DSYTRS_3_TEST(const char *uplo, const lapack_int *n,
                   const lapack_int *nrhs, const double *a,
                   const lapack_int *lda, const double *e,
                   const lapack_int *ipiv, double *b, const lapack_int *ldb,
                   lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DSYTRS_3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytrs_3)(layout, *uplo, *n, *nrhs, a_r, lda_r, e,
                                       ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsytrs_3_work)(layout, *uplo, *n, *nrhs, a_r,
                                            lda_r, e, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DSYTRS_3", ret);
}

/******************************************************************************
 * DSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define DSYTRS_AA_TEST LAPACK_GLOBAL_SUFFIX(dsytrs_aa_test, DSYTRS_AA_TEST)
void DSYTRS_AA_TEST(const char *uplo, const lapack_int *n,
                    const lapack_int *nrhs, const double *a,
                    const lapack_int *lda, const lapack_int *ipiv, double *b,
                    const lapack_int *ldb, double *work,
                    const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dsytrs_aa_work)(LAPACK_COL_MAJOR, *uplo, *n,
                                                 *nrhs, a, *lda, ipiv, b, *ldb,
                                                 work, *lwork);
        *info = lapacke_test_info("DSYTRS_AA", ret);
        return;
    }

    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DSYTRS_AA", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytrs_aa)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                        ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsytrs_aa_work)(
        layout, *uplo, *n, *nrhs, a_r, lda_r, ipiv, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DSYTRS_AA", ret);
}

/******************************************************************************
 * DSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DSYTRS_ROOK_TEST                                                       \
    LAPACK_GLOBAL_SUFFIX(dsytrs_rook_test, DSYTRS_ROOK_TEST)
void DSYTRS_ROOK_TEST(const char *uplo, const lapack_int *n,
                      const lapack_int *nrhs, const double *a,
                      const lapack_int *lda, const lapack_int *ipiv, double *b,
                      const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                      ,
                      FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DSYTRS_ROOK", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsytrs_rook)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsytrs_rook_work)(layout, *uplo, *n, *nrhs, a_r,
                                               lda_r, ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DSYTRS_ROOK", ret);
}

// ========================================================================== //
//                  Packed symmetric indefinite matrices (SP)                 //
// ========================================================================== //

/******************************************************************************
 * DSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DSPCON_TEST LAPACK_GLOBAL_SUFFIX(dspcon_test, DSPCON_TEST)
void DSPCON_TEST(const char *uplo, const lapack_int *n, const double *ap,
                 const lapack_int *ipiv, const double *anorm, double *rcond,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_dsp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("DSPCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dspcon)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                     rcond);
#else
    ret = API_SUFFIX(LAPACKE_dspcon_work)(layout, *uplo, *n, ap_r, ipiv, *anorm,
                                          rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("DSPCON", ret);
}

/******************************************************************************
 * DSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK,
 * IWORK, INFO )
 ******************************************************************************/
#define DSPRFS_TEST LAPACK_GLOBAL_SUFFIX(dsprfs_test, DSPRFS_TEST)
void DSPRFS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *ap, const double *afp, const lapack_int *ipiv,
                 const double *b, const lapack_int *ldb, double *x,
                 const lapack_int *ldx, double *ferr, double *berr,
                 double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
    const double *ap_r = ap;
    const double *afp_r = afp;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    ap_r = lapacke_test_dsp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_dsp_cm_to_rm(*uplo, *n, afp);
    if (b_r == NULL || x_r == NULL || ap_r == NULL || afp_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free(x_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free((void *)afp_r);
        lapacke_test_report_alloc_failure("DSPRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsprfs)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                     ipiv, b_r, ldb_r, x_r, ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dsprfs_work)(layout, *uplo, *n, *nrhs, ap_r, afp_r,
                                          ipiv, b_r, ldb_r, x_r, ldx_r, ferr,
                                          berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free(x_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free((void *)afp_r);
#endif
    *info = lapacke_test_info("DSPRFS", ret);
}

/******************************************************************************
 * DSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DSPSV_TEST LAPACK_GLOBAL_SUFFIX(dspsv_test, DSPSV_TEST)
void DSPSV_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                double *ap, lapack_int *ipiv, double *b, const lapack_int *ldb,
                lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
    double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_dsp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free(ap_r);
        lapacke_test_report_alloc_failure("DSPSV", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dspsv)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                    ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dspsv_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                         b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    lapacke_test_dsp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(b_r);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("DSPSV", ret);
}

/******************************************************************************
 * DSPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR,
 * BERR, WORK, IWORK, INFO )
 ******************************************************************************/
#define DSPSVX_TEST LAPACK_GLOBAL_SUFFIX(dspsvx_test, DSPSVX_TEST)
void DSPSVX_TEST(const char *fact, const char *uplo, const lapack_int *n,
                 const lapack_int *nrhs, const double *ap, double *afp,
                 lapack_int *ipiv, const double *b, const lapack_int *ldb,
                 double *x, const lapack_int *ldx, double *rcond, double *ferr,
                 double *berr, double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN fact_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    const double *ap_r = ap;
    double *afp_r = afp;
    double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_dsp_cm_to_rm(*uplo, *n, ap);
    afp_r = lapacke_test_dsp_cm_to_rm(*uplo, *n, afp);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (b_r == NULL || ap_r == NULL || afp_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)ap_r);
        LAPACKE_free(afp_r);
        LAPACKE_free(x_r);
        lapacke_test_report_alloc_failure("DSPSVX", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dspsvx)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                     afp_r, ipiv, b_r, ldb_r, x_r, ldx_r, rcond,
                                     ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dspsvx_work)(layout, *fact, *uplo, *n, *nrhs, ap_r,
                                          afp_r, ipiv, b_r, ldb_r, x_r, ldx_r,
                                          rcond, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsp_rm_to_cm(*uplo, *n, afp_r, afp);
    lapacke_test_dge_rm_to_cm(*n, *nrhs, x_r, ldx_r, x, *ldx);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)ap_r);
    LAPACKE_free(afp_r);
    LAPACKE_free(x_r);
#endif
    *info = lapacke_test_info("DSPSVX", ret);
}

/******************************************************************************
 * DSPTRF( UPLO, N, AP, IPIV, INFO )
 ******************************************************************************/
#define DSPTRF_TEST LAPACK_GLOBAL_SUFFIX(dsptrf_test, DSPTRF_TEST)
void DSPTRF_TEST(const char *uplo, const lapack_int *n, double *ap,
                 lapack_int *ipiv, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_dsp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("DSPTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsptrf)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dsptrf_work)(layout, *uplo, *n, ap_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("DSPTRF", ret);
}

/******************************************************************************
 * DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
 ******************************************************************************/
#define DSPTRI_TEST LAPACK_GLOBAL_SUFFIX(dsptri_test, DSPTRI_TEST)
void DSPTRI_TEST(const char *uplo, const lapack_int *n, double *ap,
                 const lapack_int *ipiv, double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    ap_r = lapacke_test_dsp_cm_to_rm(*uplo, *n, ap);
    if (ap_r == NULL) {
        lapacke_test_report_alloc_failure("DSPTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsptri)(layout, *uplo, *n, ap_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dsptri_work)(layout, *uplo, *n, ap_r, ipiv, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dsp_rm_to_cm(*uplo, *n, ap_r, ap);
    LAPACKE_free(ap_r);
#endif
    *info = lapacke_test_info("DSPTRI", ret);
}

/******************************************************************************
 * DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 ******************************************************************************/
#define DSPTRS_TEST LAPACK_GLOBAL_SUFFIX(dsptrs_test, DSPTRS_TEST)
void DSPTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *ap, const lapack_int *ipiv, double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
    const double *ap_r = ap;
#if LAPACKE_TEST_ROW_MAJOR
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    ap_r = lapacke_test_dsp_cm_to_rm(*uplo, *n, ap);
    if (b_r == NULL || ap_r == NULL) {
        LAPACKE_free(b_r);
        LAPACKE_free((void *)ap_r);
        lapacke_test_report_alloc_failure("DSPTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dsptrs)(layout, *uplo, *n, *nrhs, ap_r, ipiv, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dsptrs_work)(layout, *uplo, *n, *nrhs, ap_r, ipiv,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(b_r);
    LAPACKE_free((void *)ap_r);
#endif
    *info = lapacke_test_info("DSPTRS", ret);
}

// ========================================================================== //
//                          Triangular matrices (TR)                          //
// ========================================================================== //

/******************************************************************************
 * DTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DTRCON_TEST LAPACK_GLOBAL_SUFFIX(dtrcon_test, DTRCON_TEST)
void DTRCON_TEST(const char *norm, const char *uplo, const char *diag,
                 const lapack_int *n, const double *a, const lapack_int *lda,
                 double *rcond, double *work, lapack_int *iwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dtr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DTRCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dtrcon)(layout, *norm, *uplo, *diag, *n, a_r,
                                     lda_r, rcond);
#else
    ret = API_SUFFIX(LAPACKE_dtrcon_work)(layout, *norm, *uplo, *diag, *n, a_r,
                                          lda_r, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    *info = lapacke_test_info("DTRCON", ret);
}

/******************************************************************************
 * DTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, LDX, FERR, BERR, WORK,
 * IWORK, INFO )
 ******************************************************************************/
#define DTRRFS_TEST LAPACK_GLOBAL_SUFFIX(dtrrfs_test, DTRRFS_TEST)
void DTRRFS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *nrhs, const double *a,
                 const lapack_int *lda, const double *b, const lapack_int *ldb,
                 const double *x, const lapack_int *ldx, double *ferr,
                 double *berr, double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    const double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dtr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (a_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)x_r);
        lapacke_test_report_alloc_failure("DTRRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dtrrfs)(layout, *uplo, *trans, *diag, *n, *nrhs,
                                     a_r, lda_r, b_r, ldb_r, x_r, ldx_r, ferr,
                                     berr);
#else
    ret = API_SUFFIX(LAPACKE_dtrrfs_work)(layout, *uplo, *trans, *diag, *n,
                                          *nrhs, a_r, lda_r, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)x_r);
#endif
    *info = lapacke_test_info("DTRRFS", ret);
}

/******************************************************************************
 * DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
 ******************************************************************************/
#define DTRTRI_TEST LAPACK_GLOBAL_SUFFIX(dtrtri_test, DTRTRI_TEST)
void DTRTRI_TEST(const char *uplo, const char *diag, const lapack_int *n,
                 double *a, const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dtr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DTRTRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dtrtri)(layout, *uplo, *diag, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_dtrtri_work)(layout, *uplo, *diag, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dtr_rm_to_cm(*uplo, *diag, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DTRTRI", ret);
}

/******************************************************************************
 * DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
 ******************************************************************************/
#define DTRTRS_TEST LAPACK_GLOBAL_SUFFIX(dtrtrs_test, DTRTRS_TEST)
void DTRTRS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *nrhs, const double *a,
                 const lapack_int *lda, double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dtr_cm_to_rm(*uplo, *diag, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DTRTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dtrtrs)(layout, *uplo, *trans, *diag, *n, *nrhs,
                                     a_r, lda_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dtrtrs_work)(layout, *uplo, *trans, *diag, *n,
                                          *nrhs, a_r, lda_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DTRTRS", ret);
}

// ========================================================================== //
//                        Triangular band matrices (TB)                       //
// ========================================================================== //

/******************************************************************************
 * DTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, IWORK, INFO )
 ******************************************************************************/
#define DTBCON_TEST LAPACK_GLOBAL_SUFFIX(dtbcon_test, DTBCON_TEST)
void DTBCON_TEST(const char *norm, const char *uplo, const char *diag,
                 const lapack_int *n, const lapack_int *kd, const double *ab,
                 const lapack_int *ldab, double *rcond, double *work,
                 lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dtb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    if (ab_r == NULL) {
        lapacke_test_report_alloc_failure("DTBCON", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dtbcon)(layout, *norm, *uplo, *diag, *n, *kd, ab_r,
                                     ldab_r, rcond);
#else
    ret = API_SUFFIX(LAPACKE_dtbcon_work)(layout, *norm, *uplo, *diag, *n, *kd,
                                          ab_r, ldab_r, rcond, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
#endif
    *info = lapacke_test_info("DTBCON", ret);
}

/******************************************************************************
 * DTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, X, LDX, FERR, BERR,
 * WORK, IWORK, INFO )
 ******************************************************************************/
#define DTBRFS_TEST LAPACK_GLOBAL_SUFFIX(dtbrfs_test, DTBRFS_TEST)
void DTBRFS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const double *ab,
                 const lapack_int *ldab, const double *b, const lapack_int *ldb,
                 const double *x, const lapack_int *ldx, double *ferr,
                 double *berr, double *work, lapack_int *iwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    const double *b_r = b;
    lapack_int ldb_r = *ldb;
    const double *x_r = x;
    lapack_int ldx_r = *ldx;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dtb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    x_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, x, *ldx, &ldx_r);
    if (ab_r == NULL || b_r == NULL || x_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free((void *)b_r);
        LAPACKE_free((void *)x_r);
        lapacke_test_report_alloc_failure("DTBRFS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dtbrfs)(layout, *uplo, *trans, *diag, *n, *kd,
                                     *nrhs, ab_r, ldab_r, b_r, ldb_r, x_r,
                                     ldx_r, ferr, berr);
#else
    ret = API_SUFFIX(LAPACKE_dtbrfs_work)(layout, *uplo, *trans, *diag, *n, *kd,
                                          *nrhs, ab_r, ldab_r, b_r, ldb_r, x_r,
                                          ldx_r, ferr, berr, work, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)ab_r);
    LAPACKE_free((void *)b_r);
    LAPACKE_free((void *)x_r);
#endif
    *info = lapacke_test_info("DTBRFS", ret);
}

/******************************************************************************
 * DTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 ******************************************************************************/
#define DTBTRS_TEST LAPACK_GLOBAL_SUFFIX(dtbtrs_test, DTBTRS_TEST)
void DTBTRS_TEST(const char *uplo, const char *trans, const char *diag,
                 const lapack_int *n, const lapack_int *kd,
                 const lapack_int *nrhs, const double *ab,
                 const lapack_int *ldab, double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len, FORTRAN_STRLEN trans_len,
                 FORTRAN_STRLEN diag_len
#endif
)
{
    lapack_int ret = 0;
    const double *ab_r = ab;
    lapack_int ldab_r = *ldab;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    ab_r = lapacke_test_dtb_cm_to_rm(*uplo, *diag, *n, *kd, ab, *ldab, &ldab_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (ab_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)ab_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DTBTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dtbtrs)(layout, *uplo, *trans, *diag, *n, *kd,
                                     *nrhs, ab_r, ldab_r, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dtbtrs_work)(layout, *uplo, *trans, *diag, *n, *kd,
                                          *nrhs, ab_r, ldab_r, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)ab_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DTBTRS", ret);
}

// ========================================================================== //
//                   Orthogonal factorizations (QR/LQ/QL/RQ)                  //
// ========================================================================== //

/******************************************************************************
 * DGELQ2( M, N, A, LDA, TAU, WORK, INFO )
 ******************************************************************************/
#define DGELQ2_TEST LAPACK_GLOBAL_SUFFIX(dgelq2_test, DGELQ2_TEST)
void DGELQ2_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, double *tau, double *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGELQ2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgelq2)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dgelq2_work)(layout, *m, *n, a_r, lda_r, tau,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGELQ2", ret);
}

/******************************************************************************
 * DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGELQF_TEST LAPACK_GLOBAL_SUFFIX(dgelqf_test, DGELQF_TEST)
void DGELQF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, double *tau, double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgelqf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("DGELQF", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGELQF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgelqf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dgelqf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGELQF", ret);
}

/******************************************************************************
 * DGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGEQLF_TEST LAPACK_GLOBAL_SUFFIX(dgeqlf_test, DGEQLF_TEST)
void DGEQLF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, double *tau, double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgeqlf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("DGEQLF", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGEQLF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqlf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dgeqlf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGEQLF", ret);
}

/******************************************************************************
 * DGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGEQR_TEST LAPACK_GLOBAL_SUFFIX(dgeqr_test, DGEQR_TEST)
void DGEQR_TEST(const lapack_int *m, const lapack_int *n, double *a,
                const lapack_int *lda, double *t, const lapack_int *tsize,
                double *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*tsize == -1 || *lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgeqr_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                             t, *tsize, work, *lwork);
        *info = lapacke_test_info("DGEQR", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGEQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqr)(layout, *m, *n, a_r, lda_r, t, *tsize);
#else
    ret = API_SUFFIX(LAPACKE_dgeqr_work)(layout, *m, *n, a_r, lda_r, t, *tsize,
                                         work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGEQR", ret);
}

/******************************************************************************
 * DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
 ******************************************************************************/
#define DGEQR2_TEST LAPACK_GLOBAL_SUFFIX(dgeqr2_test, DGEQR2_TEST)
void DGEQR2_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, double *tau, double *work,
                 lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGEQR2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqr2)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dgeqr2_work)(layout, *m, *n, a_r, lda_r, tau,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGEQR2", ret);
}

/******************************************************************************
 * DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGEQRF_TEST LAPACK_GLOBAL_SUFFIX(dgeqrf_test, DGEQRF_TEST)
void DGEQRF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, double *tau, double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgeqrf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("DGEQRF", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGEQRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqrf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dgeqrf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGEQRF", ret);
}

/******************************************************************************
 * DGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGEQRFP_TEST LAPACK_GLOBAL_SUFFIX(dgeqrfp_test, DGEQRFP_TEST)
void DGEQRFP_TEST(const lapack_int *m, const lapack_int *n, double *a,
                  const lapack_int *lda, double *tau, double *work,
                  const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgeqrfp_work)(LAPACK_COL_MAJOR, *m, *n, a,
                                               *lda, tau, work, *lwork);
        *info = lapacke_test_info("DGEQRFP", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGEQRFP", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqrfp)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dgeqrfp_work)(layout, *m, *n, a_r, lda_r, tau,
                                           work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGEQRFP", ret);
}

/******************************************************************************
 * DGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGERQF_TEST LAPACK_GLOBAL_SUFFIX(dgerqf_test, DGERQF_TEST)
void DGERQF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, double *tau, double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgerqf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("DGERQF", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGERQF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgerqf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dgerqf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGERQF", ret);
}

/******************************************************************************
 * DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DORGLQ_TEST LAPACK_GLOBAL_SUFFIX(dorglq_test, DORGLQ_TEST)
void DORGLQ_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 double *a, const lapack_int *lda, const double *tau,
                 double *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dorglq_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("DORGLQ", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DORGLQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dorglq)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dorglq_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DORGLQ", ret);
}

/******************************************************************************
 * DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DORGQL_TEST LAPACK_GLOBAL_SUFFIX(dorgql_test, DORGQL_TEST)
void DORGQL_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 double *a, const lapack_int *lda, const double *tau,
                 double *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dorgql_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("DORGQL", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DORGQL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dorgql)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dorgql_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DORGQL", ret);
}

/******************************************************************************
 * DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DORGQR_TEST LAPACK_GLOBAL_SUFFIX(dorgqr_test, DORGQR_TEST)
void DORGQR_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 double *a, const lapack_int *lda, const double *tau,
                 double *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dorgqr_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("DORGQR", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DORGQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dorgqr)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dorgqr_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DORGQR", ret);
}

/******************************************************************************
 * DORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DORGRQ_TEST LAPACK_GLOBAL_SUFFIX(dorgrq_test, DORGRQ_TEST)
void DORGRQ_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 double *a, const lapack_int *lda, const double *tau,
                 double *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dorgrq_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("DORGRQ", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DORGRQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dorgrq)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dorgrq_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DORGRQ", ret);
}

/******************************************************************************
 * DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define DORMLQ_TEST LAPACK_GLOBAL_SUFFIX(dormlq_test, DORMLQ_TEST)
void DORMLQ_TEST(const char *side, const char *trans, const lapack_int *m,
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
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dormlq_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("DORMLQ", ret);
        return;
    }

    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_dge_cm_to_rm(*k, nq, a, *lda, &lda_r);
    c_r = lapacke_test_dge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("DORMLQ", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dormlq)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_dormlq_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("DORMLQ", ret);
}

/******************************************************************************
 * DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define DORMQL_TEST LAPACK_GLOBAL_SUFFIX(dormql_test, DORMQL_TEST)
void DORMQL_TEST(const char *side, const char *trans, const lapack_int *m,
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
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dormql_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("DORMQL", ret);
        return;
    }

    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_dge_cm_to_rm(nq, *k, a, *lda, &lda_r);
    c_r = lapacke_test_dge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("DORMQL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dormql)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_dormql_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("DORMQL", ret);
}

/******************************************************************************
 * DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 ******************************************************************************/
#define DORMQR_TEST LAPACK_GLOBAL_SUFFIX(dormqr_test, DORMQR_TEST)
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
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dormqr_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("DORMQR", ret);
        return;
    }

    const double *a_r = a;
    lapack_int lda_r = *lda;
    double *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int nq = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_dge_cm_to_rm(nq, *k, a, *lda, &lda_r);
    c_r = lapacke_test_dge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("DORMQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dormqr)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_dormqr_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("DORMQR", ret);
}

// ========================================================================== //
//            Blocked and tall-skinny factorizations (QRT/TSQR/HR)            //
// ========================================================================== //

/******************************************************************************
 * DGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
 ******************************************************************************/
#define DGEQRT_TEST LAPACK_GLOBAL_SUFFIX(dgeqrt_test, DGEQRT_TEST)
void DGEQRT_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *nb,
                 double *a, const lapack_int *lda, double *t,
                 const lapack_int *ldt, double *work, lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
    double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_dge_cm_to_rm(*nb, MIN(*m, *n), t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("DGEQRT", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqrt)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                     ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_dgeqrt_work)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                          ldt_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*nb, MIN(*m, *n), t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("DGEQRT", ret);
}

/******************************************************************************
 * DGEQRT2( M, N, A, LDA, T, LDT, INFO )
 ******************************************************************************/
#define DGEQRT2_TEST LAPACK_GLOBAL_SUFFIX(dgeqrt2_test, DGEQRT2_TEST)
void DGEQRT2_TEST(const lapack_int *m, const lapack_int *n, double *a,
                  const lapack_int *lda, double *t, const lapack_int *ldt,
                  lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
    double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_dge_cm_to_rm(*n, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("DGEQRT2", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqrt2)(layout, *m, *n, a_r, lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_dgeqrt2_work)(layout, *m, *n, a_r, lda_r, t_r,
                                           ldt_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*n, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("DGEQRT2", ret);
}

/******************************************************************************
 * DGEQRT3( M, N, A, LDA, T, LDT, INFO )
 ******************************************************************************/
#define DGEQRT3_TEST LAPACK_GLOBAL_SUFFIX(dgeqrt3_test, DGEQRT3_TEST)
void DGEQRT3_TEST(const lapack_int *m, const lapack_int *n, double *a,
                  const lapack_int *lda, double *t, const lapack_int *ldt,
                  lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
    double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_dge_cm_to_rm(*n, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("DGEQRT3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqrt3)(layout, *m, *n, a_r, lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_dgeqrt3_work)(layout, *m, *n, a_r, lda_r, t_r,
                                           ldt_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*n, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("DGEQRT3", ret);
}

/******************************************************************************
 * DGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGETSQRHRT_TEST LAPACK_GLOBAL_SUFFIX(dgetsqrhrt_test, DGETSQRHRT_TEST)
void DGETSQRHRT_TEST(const lapack_int *m, const lapack_int *n,
                     const lapack_int *mb1, const lapack_int *nb1,
                     const lapack_int *nb2, double *a, const lapack_int *lda,
                     double *t, const lapack_int *ldt, double *work,
                     const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgetsqrhrt_work)(LAPACK_COL_MAJOR, *m, *n,
                                                  *mb1, *nb1, *nb2, a, *lda, t,
                                                  *ldt, work, *lwork);
        *info = lapacke_test_info("DGETSQRHRT", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_dge_cm_to_rm(*nb2, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("DGETSQRHRT", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgetsqrhrt)(layout, *m, *n, *mb1, *nb1, *nb2, a_r,
                                         lda_r, t_r, ldt_r);
#else
    ret = API_SUFFIX(LAPACKE_dgetsqrhrt_work)(
        layout, *m, *n, *mb1, *nb1, *nb2, a_r, lda_r, t_r, ldt_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*nb2, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("DGETSQRHRT", ret);
}

/******************************************************************************
 * DORHR_COL( M, N, NB, A, LDA, T, LDT, D, INFO )
 ******************************************************************************/
#define DORHR_COL_TEST LAPACK_GLOBAL_SUFFIX(dorhr_col_test, DORHR_COL_TEST)
void DORHR_COL_TEST(const lapack_int *m, const lapack_int *n,
                    const lapack_int *nb, double *a, const lapack_int *lda,
                    double *t, const lapack_int *ldt, double *d,
                    lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
    double *t_r = t;
    lapack_int ldt_r = *ldt;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    t_r = lapacke_test_dge_cm_to_rm(*ldt, *n, t, *ldt, &ldt_r);
    if (a_r == NULL || t_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(t_r);
        lapacke_test_report_alloc_failure("DORHR_COL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dorhr_col)(layout, *m, *n, *nb, a_r, lda_r, t_r,
                                        ldt_r, d);
#else
    ret = API_SUFFIX(LAPACKE_dorhr_col_work)(layout, *m, *n, *nb, a_r, lda_r,
                                             t_r, ldt_r, d);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(*ldt, *n, t_r, ldt_r, t, *ldt);
    LAPACKE_free(a_r);
    LAPACKE_free(t_r);
#endif
    *info = lapacke_test_info("DORHR_COL", ret);
}

// ========================================================================== //
//                    Rank-revealing factorizations (QP/TZ)                   //
// ========================================================================== //

/******************************************************************************
 * DGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGEQP3_TEST LAPACK_GLOBAL_SUFFIX(dgeqp3_test, DGEQP3_TEST)
void DGEQP3_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *jpvt, double *tau,
                 double *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgeqp3_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              jpvt, tau, work, *lwork);
        *info = lapacke_test_info("DGEQP3", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGEQP3", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqp3)(layout, *m, *n, a_r, lda_r, jpvt, tau);
#else
    ret = API_SUFFIX(LAPACKE_dgeqp3_work)(layout, *m, *n, a_r, lda_r, jpvt, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGEQP3", ret);
}

/******************************************************************************
 * DTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 ******************************************************************************/
#define DTZRZF_TEST LAPACK_GLOBAL_SUFFIX(dtzrzf_test, DTZRZF_TEST)
void DTZRZF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, double *tau, double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dtzrzf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("DTZRZF", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DTZRZF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dtzrzf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dtzrzf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DTZRZF", ret);
}

// ========================================================================== //
//                             Least squares (LS)                             //
// ========================================================================== //

/******************************************************************************
 * DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGELS_TEST LAPACK_GLOBAL_SUFFIX(dgels_test, DGELS_TEST)
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
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgels_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                             *nrhs, a, *lda, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("DGELS", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGELS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgels)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgels_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                         lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGELS", ret);
}

/******************************************************************************
 * DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO
 * )
 ******************************************************************************/
#define DGELSD_TEST LAPACK_GLOBAL_SUFFIX(dgelsd_test, DGELSD_TEST)
void DGELSD_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, double *a, const lapack_int *lda,
                 double *b, const lapack_int *ldb, double *s,
                 const double *rcond, lapack_int *rank, double *work,
                 const lapack_int *lwork, lapack_int *iwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgelsd_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, s, *rcond, rank,
                                              work, *lwork, iwork);
        *info = lapacke_test_info("DGELSD", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGELSD", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgelsd)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, s, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_dgelsd_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, s, *rcond, rank, work,
                                          *lwork, iwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGELSD", ret);
}

/******************************************************************************
 * DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGELSS_TEST LAPACK_GLOBAL_SUFFIX(dgelss_test, DGELSS_TEST)
void DGELSS_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, double *a, const lapack_int *lda,
                 double *b, const lapack_int *ldb, double *s,
                 const double *rcond, lapack_int *rank, double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgelss_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, s, *rcond, rank,
                                              work, *lwork);
        *info = lapacke_test_info("DGELSS", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGELSS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgelss)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, s, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_dgelss_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, s, *rcond, rank, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGELSS", ret);
}

/******************************************************************************
 * DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGELSY_TEST LAPACK_GLOBAL_SUFFIX(dgelsy_test, DGELSY_TEST)
void DGELSY_TEST(const lapack_int *m, const lapack_int *n,
                 const lapack_int *nrhs, double *a, const lapack_int *lda,
                 double *b, const lapack_int *ldb, lapack_int *jpvt,
                 const double *rcond, lapack_int *rank, double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgelsy_work)(LAPACK_COL_MAJOR, *m, *n, *nrhs,
                                              a, *lda, b, *ldb, jpvt, *rcond,
                                              rank, work, *lwork);
        *info = lapacke_test_info("DGELSY", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGELSY", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgelsy)(layout, *m, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r, jpvt, *rcond, rank);
#else
    ret = API_SUFFIX(LAPACKE_dgelsy_work)(layout, *m, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r, jpvt, *rcond, rank, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGELSY", ret);
}

/******************************************************************************
 * DGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGETSLS_TEST LAPACK_GLOBAL_SUFFIX(dgetsls_test, DGETSLS_TEST)
void DGETSLS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                  const lapack_int *nrhs, double *a, const lapack_int *lda,
                  double *b, const lapack_int *ldb, double *work,
                  const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgetsls_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                               *nrhs, a, *lda, b, *ldb, work,
                                               *lwork);
        *info = lapacke_test_info("DGETSLS", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *b_r = b;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(MAX(*m, *n), *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGETSLS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgetsls)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                      b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgetsls_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                           lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(MAX(*m, *n), *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGETSLS", ret);
}

// ========================================================================== //
//                             SVD and bidiagonal                             //
// ========================================================================== //

/******************************************************************************
 * DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )
 ******************************************************************************/
#define DBDSQR_TEST LAPACK_GLOBAL_SUFFIX(dbdsqr_test, DBDSQR_TEST)
void DBDSQR_TEST(const char *uplo, const lapack_int *n, const lapack_int *ncvt,
                 const lapack_int *nru, const lapack_int *ncc, double *d,
                 double *e, double *vt, const lapack_int *ldvt, double *u,
                 const lapack_int *ldu, double *c, const lapack_int *ldc,
                 double *work, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *vt_r = vt;
    lapack_int ldvt_r = *ldvt;
    double *u_r = u;
    lapack_int ldu_r = *ldu;
    double *c_r = c;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    vt_r = lapacke_test_dge_cm_to_rm(*n, *ncvt, vt, *ldvt, &ldvt_r);
    u_r = lapacke_test_dge_cm_to_rm(*nru, *n, u, *ldu, &ldu_r);
    c_r = lapacke_test_dge_cm_to_rm(*n, *ncc, c, *ldc, &ldc_r);
    if (vt_r == NULL || u_r == NULL || c_r == NULL) {
        LAPACKE_free(vt_r);
        LAPACKE_free(u_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("DBDSQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dbdsqr)(layout, *uplo, *n, *ncvt, *nru, *ncc, d, e,
                                     vt_r, ldvt_r, u_r, ldu_r, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_dbdsqr_work)(layout, *uplo, *n, *ncvt, *nru, *ncc,
                                          d, e, vt_r, ldvt_r, u_r, ldu_r, c_r,
                                          ldc_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *ncvt, vt_r, ldvt_r, vt, *ldvt);
    lapacke_test_dge_rm_to_cm(*nru, *n, u_r, ldu_r, u, *ldu);
    lapacke_test_dge_rm_to_cm(*n, *ncc, c_r, ldc_r, c, *ldc);
    LAPACKE_free(vt_r);
    LAPACKE_free(u_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("DBDSQR", ret);
}

/******************************************************************************
 * DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
 ******************************************************************************/
#define DGESVD_TEST LAPACK_GLOBAL_SUFFIX(dgesvd_test, DGESVD_TEST)
void DGESVD_TEST(const char *jobu, const char *jobvt, const lapack_int *m,
                 const lapack_int *n, double *a, const lapack_int *lda,
                 double *s, double *u, const lapack_int *ldu, double *vt,
                 const lapack_int *ldvt, double *work, const lapack_int *lwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN jobu_len, FORTRAN_STRLEN jobvt_len
#endif
)
{
    lapack_int ret = 0;
    if (*lwork == -1) { /* workspace query */
        ret = API_SUFFIX(LAPACKE_dgesvd_work)(LAPACK_COL_MAJOR, *jobu, *jobvt,
                                              *m, *n, a, *lda, s, u, *ldu, vt,
                                              *ldvt, work, *lwork);
        *info = lapacke_test_info("DGESVD", ret);
        return;
    }

    double *a_r = a;
    lapack_int lda_r = *lda;
    double *u_r = u;
    lapack_int ldu_r = *ldu;
    double *vt_r = vt;
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
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    u_r = lapacke_test_dge_cm_to_rm(urows, ucols, u, *ldu, &ldu_r);
    vt_r = lapacke_test_dge_cm_to_rm(vtrows, *n, vt, *ldvt, &ldvt_r);
    if (a_r == NULL || u_r == NULL || vt_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(u_r);
        LAPACKE_free(vt_r);
        lapacke_test_report_alloc_failure("DGESVD", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgesvd)(layout, *jobu, *jobvt, *m, *n, a_r, lda_r,
                                     s, u_r, ldu_r, vt_r, ldvt_r, work);
#else
    ret = API_SUFFIX(LAPACKE_dgesvd_work)(layout, *jobu, *jobvt, *m, *n, a_r,
                                          lda_r, s, u_r, ldu_r, vt_r, ldvt_r,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(urows, ucols, u_r, ldu_r, u, *ldu);
    lapacke_test_dge_rm_to_cm(vtrows, *n, vt_r, ldvt_r, vt, *ldvt);
    LAPACKE_free(a_r);
    LAPACKE_free(u_r);
    LAPACKE_free(vt_r);
#endif
    *info = lapacke_test_info("DGESVD", ret);
}

// ========================================================================== //
//                             Auxiliary routines                             //
// ========================================================================== //

/******************************************************************************
 * DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
 ******************************************************************************/
#define DLAMCH_TEST LAPACK_GLOBAL_SUFFIX(dlamch_test, DLAMCH_TEST)
double DLAMCH_TEST(const char *cmach
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN cmach_len
#endif
)
{
    double res = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_dlamch)(*cmach);
#else
    res = API_SUFFIX(LAPACKE_dlamch_work)(*cmach);
#endif

    return res;
}

/******************************************************************************
 * DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
 ******************************************************************************/
#define DLANGE_TEST LAPACK_GLOBAL_SUFFIX(dlange_test, DLANGE_TEST)
double DLANGE_TEST(const char *norm, const lapack_int *m, const lapack_int *n,
                   const double *a, const lapack_int *lda, double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len
#endif
)
{
    double res = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("DLANGE", &info);
        return (double)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_dlange)(layout, *norm, *m, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_dlange_work)(layout, *norm, *m, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
 ******************************************************************************/
#define DLANSY_TEST LAPACK_GLOBAL_SUFFIX(dlansy_test, DLANSY_TEST)
double DLANSY_TEST(const char *norm, const char *uplo, const lapack_int *n,
                   const double *a, const lapack_int *lda, double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len
#endif
)
{
    double res = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dsy_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("DLANSY", &info);
        return (double)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_dlansy)(layout, *norm, *uplo, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_dlansy_work)(layout, *norm, *uplo, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * DOUBLE PRECISION FUNCTION DLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
 ******************************************************************************/
#define DLANTR_TEST LAPACK_GLOBAL_SUFFIX(dlantr_test, DLANTR_TEST)
double DLANTR_TEST(const char *norm, const char *uplo, const char *diag,
                   const lapack_int *m, const lapack_int *n, const double *a,
                   const lapack_int *lda, double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len, FORTRAN_STRLEN uplo_len,
                   FORTRAN_STRLEN diag_len
#endif
)
{
    double res = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("DLANTR", &info);
        return (double)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_dlantr)(layout, *norm, *uplo, *diag, *m, *n, a_r,
                                     lda_r);
#else
    res = API_SUFFIX(LAPACKE_dlantr_work)(layout, *norm, *uplo, *diag, *m, *n,
                                          a_r, lda_r, work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}

/******************************************************************************
 * DLARNV( IDIST, ISEED, N, X )
 ******************************************************************************/
#define DLARNV_TEST LAPACK_GLOBAL_SUFFIX(dlarnv_test, DLARNV_TEST)
void DLARNV_TEST(const lapack_int *idist, lapack_int *iseed,
                 const lapack_int *n, double *x)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dlarnv)(*idist, iseed, *n, x);
#else
    ret = API_SUFFIX(LAPACKE_dlarnv_work)(*idist, iseed, *n, x);
#endif

    lapacke_test_info_unshifted("DLARNV", ret);
}

/******************************************************************************
 * DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
 ******************************************************************************/
#define DLASCL_TEST LAPACK_GLOBAL_SUFFIX(dlascl_test, DLASCL_TEST)
void DLASCL_TEST(const char *type, const lapack_int *kl, const lapack_int *ku,
                 const double *cfrom, const double *cto, const lapack_int *m,
                 const lapack_int *n, double *a, const lapack_int *lda,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN type_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int arows = API_SUFFIX(LAPACKE_lsame)(*type, 'b')   ? *kl + 1
                             : API_SUFFIX(LAPACKE_lsame)(*type, 'q') ? *ku + 1
                             : API_SUFFIX(LAPACKE_lsame)(*type, 'z')
                                 ? 2 * *kl + *ku + 1
                                 : *m;
    a_r = lapacke_test_dge_cm_to_rm(arows, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DLASCL", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dlascl)(layout, *type, *kl, *ku, *cfrom, *cto, *m,
                                     *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_dlascl_work)(layout, *type, *kl, *ku, *cfrom, *cto,
                                          *m, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(arows, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DLASCL", ret);
}

/******************************************************************************
 * DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
 ******************************************************************************/
#define DLASET_TEST LAPACK_GLOBAL_SUFFIX(dlaset_test, DLASET_TEST)
void DLASET_TEST(const char *uplo, const lapack_int *m, const lapack_int *n,
                 const double *alpha, const double *beta, double *a,
                 const lapack_int *lda
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("DLASET", &info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dlaset)(layout, *uplo, *m, *n, *alpha, *beta, a_r,
                                     lda_r);
#else
    ret = API_SUFFIX(LAPACKE_dlaset_work)(layout, *uplo, *m, *n, *alpha, *beta,
                                          a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    lapacke_test_info("DLASET", ret);
}

/******************************************************************************
 * DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
 ******************************************************************************/
#define DLASWP_TEST LAPACK_GLOBAL_SUFFIX(dlaswp_test, DLASWP_TEST)
void DLASWP_TEST(const lapack_int *n, double *a, const lapack_int *lda,
                 const lapack_int *k1, const lapack_int *k2,
                 const lapack_int *ipiv, const lapack_int *incx)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(
        lapacke_test_laswp_rows(*k1, *k2, ipiv, *incx), *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("DLASWP", &info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dlaswp)(layout, *n, a_r, lda_r, *k1, *k2, ipiv,
                                     *incx);
#else
    ret = API_SUFFIX(LAPACKE_dlaswp_work)(layout, *n, a_r, lda_r, *k1, *k2,
                                          ipiv, *incx);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(lapacke_test_laswp_rows(*k1, *k2, ipiv, *incx),
                              *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    lapacke_test_info("DLASWP", ret);
}
