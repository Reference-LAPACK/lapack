/******************************************************************************
 * Test matrix generation for the dedicated high-level LAPACKE tests
 * (double precision).
 * See lapacke_test.h.
 ******************************************************************************/

#include <stdio.h>

#include "lapacke_test.h"

/**
 * \brief Generate a random, well-conditioned rows-by-cols matrix.
 *
 * Uses ?LATMS (the generator the LAPACK test suite uses) to generate the
 * matrix in column-major order and scatters it into the target layout;
 * padding positions get a finite sentinel. The seed is fixed per fill type,
 * so every call produces the same matrix and test failures are
 * reproducible. A ?LATMS failure is recorded as a failed check.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  rows   Number of matrix rows (at most LAPACKE_TEST_LD).
 * \param[in]  cols   Number of matrix columns (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 * \param[in]  sym    ?LATMS SYM parameter: 'N' for a general matrix, 'P'
 *                    for symmetric positive definite, 'S' for symmetric
 *                    with random-sign eigenvalues (indefinite).
 * \param[in]  kl     Lower bandwidth (0 gives an upper triangular matrix).
 * \param[in]  ku     Upper bandwidth (0 gives a lower triangular matrix).
 * \param[in]  seed   The ?LATMS ISEED to start from (not modified; the
 *                    last entry must be odd).
 */
static void lapacke_test_dfill_latms(int layout, lapack_int rows,
                                     lapack_int cols, double *a, lapack_int ld,
                                     char sym, lapack_int kl, lapack_int ku,
                                     const lapack_int seed[4])
{
    double tmp[LAPACKE_TEST_LD * LAPACKE_TEST_LD];
    double d[LAPACKE_TEST_LD];
    double work[3 * LAPACKE_TEST_LD];
    lapack_int iseed[4];
    const char dist = 'U';
    const char pack = 'N';
    const lapack_int mode = 3;
    const double cond = 10.0;
    const double dmax = 1.0;
    const lapack_int ld_tmp = LAPACKE_TEST_LD;
    lapack_int info = 0;
    const size_t len = lapacke_test_alloc_len(layout, rows, cols, ld);
    size_t p;
    lapack_int i, j;

    iseed[0] = seed[0];
    iseed[1] = seed[1];
    iseed[2] = seed[2];
    iseed[3] = seed[3];
    LAPACK_dlatms(&rows, &cols, &dist, iseed, &sym, d, &mode, &cond, &dmax, &kl,
                  &ku, &pack, tmp, &ld_tmp, work, &info);
    if (info != 0) {
        lapacke_test_checks++;
        lapacke_test_failures++;
        printf("FAIL matrix generation: info = %d\n", (int)info);
    }

    for (p = 0; p < len; p++) {
        if (lapacke_test_map_position(layout, rows, cols, ld, p, &i, &j)) {
            a[p] = tmp[i + j * LAPACKE_TEST_LD];
        } else {
            a[p] = 0.25; /* padding, never referenced */
        }
    }
}

/**
 * \brief Fill a with a random general rows-by-cols matrix, double precision.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  rows   Number of matrix rows (at most LAPACKE_TEST_LD).
 * \param[in]  cols   Number of matrix columns (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill(int layout, lapack_int rows, lapack_int cols, double *a,
                        lapack_int ld)
{
    static const lapack_int seed[4] = {1988, 1989, 1990, 1991};
    lapacke_test_dfill_latms(layout, rows, cols, a, ld, 'N', rows - 1, cols - 1,
                             seed);
}

/**
 * \brief Fill a with a random symmetric positive definite n-by-n matrix, double
 * precision.
 *
 * The eigenvalues lie in [0.1, 1], so the matrix is safely positive
 * definite for the Cholesky-based tests.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill_spd(int layout, lapack_int n, double *a, lapack_int ld)
{
    static const lapack_int seed[4] = {1990, 1991, 1992, 1993};
    lapacke_test_dfill_latms(layout, n, n, a, ld, 'P', n - 1, n - 1, seed);
}

/**
 * \brief Fill a with a random symmetric indefinite n-by-n matrix, double
 * precision.
 *
 * The eigenvalues have magnitudes in [0.1, 1] and random (seed-fixed)
 * signs, so the matrix is well conditioned but not positive definite --
 * suitable for the symmetric indefinite (sy) factorizations.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill_sym(int layout, lapack_int n, double *a, lapack_int ld)
{
    static const lapack_int seed[4] = {2010, 2011, 2012, 2013};
    lapacke_test_dfill_latms(layout, n, n, a, ld, 'S', n - 1, n - 1, seed);
}

/**
 * \brief Fill a with a random nonsingular triangular n-by-n matrix, double
 * precision.
 *
 * Generated by restricting the bandwidth of a general ?LATMS matrix, so
 * the singular values stay in [0.1, 1]: the matrix is well conditioned
 * and nonsingular, and the entries beyond the uplo triangle are exactly
 * zero. Routines with a unit-diagonal option do not read the diagonal;
 * that is a question for their access region, not for the fill.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  uplo   'U' for upper, 'L' for lower triangular.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill_tri(int layout, char uplo, lapack_int n, double *a,
                            lapack_int ld)
{
    static const lapack_int seed[4] = {2020, 2021, 2022, 2023};
    const int upper = (uplo == 'U' || uplo == 'u');
    lapacke_test_dfill_latms(layout, n, n, a, ld, 'N', upper ? 0 : n - 1,
                             upper ? n - 1 : 0, seed);
}

/**
 * \brief Fill b with a random right-hand side, double precision.
 *
 * Same generator as lapacke_test_dfill with a different seed, so right-hand
 * sides differ from system matrices.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  rows   Number of matrix rows (at most LAPACKE_TEST_LD).
 * \param[in]  cols   Number of matrix columns (at most LAPACKE_TEST_LD).
 * \param[out] b      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of b in the given layout.
 */
void lapacke_test_dfill_rhs(int layout, lapack_int rows, lapack_int cols,
                            double *b, lapack_int ld)
{
    static const lapack_int seed[4] = {2000, 2001, 2002, 2003};
    lapacke_test_dfill_latms(layout, rows, cols, b, ld, 'N', rows - 1, cols - 1,
                             seed);
}

/**
 * \brief Fill every allocated position of a matrix buffer with NaN, double
 * precision.
 *
 * Used to verify that a routine with NaN checking disabled does not reject
 * NaNs anywhere in the buffer.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  rows   Number of matrix rows.
 * \param[in]  cols   Number of matrix columns.
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill_nan(int layout, lapack_int rows, lapack_int cols,
                            double *a, lapack_int ld)
{
    const size_t len = lapacke_test_alloc_len(layout, rows, cols, ld);
    for (size_t p = 0; p < len; p++) {
        a[p] = lapacke_create_nan();
    }
}

/**
 * \brief A quiet NaN, double precision.
 *
 * \return A quiet NaN.
 */
double lapacke_test_dnan(void)
{
    return lapacke_create_nan();
}

/**
 * \brief Fill ab with a random m-by-n general band matrix in band storage,
 * double precision.
 *
 * A ?LATMS matrix with kl subdiagonals and ku superdiagonals, scattered into
 * the last kl + ku + 1 rows of the rows-by-n band array (rows = 2 * kl + ku
 * + 1 leaves the kl fill-in rows a band LU factorization writes on top).
 * Fill-in rows, the unreferenced band corners and padding get a finite
 * sentinel.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  m      Number of matrix rows (at most LAPACKE_TEST_LD).
 * \param[in]  n      Number of matrix columns (at most LAPACKE_TEST_LD).
 * \param[in]  kl     Number of subdiagonals.
 * \param[in]  ku     Number of superdiagonals.
 * \param[in]  rows   Number of rows of the band array (at most
 *                    LAPACKE_TEST_LD).
 * \param[out] ab     Buffer of lapacke_test_alloc_len(layout, rows, n, ld)
 *                    elements.
 * \param[in]  ld     Leading dimension of ab in the given layout.
 */
void lapacke_test_dfill_gb(int layout, lapack_int m, lapack_int n,
                           lapack_int kl, lapack_int ku, lapack_int rows,
                           double *ab, lapack_int ld)
{
    static const lapack_int seed[4] = {2030, 2031, 2032, 2033};
    double tmp[LAPACKE_TEST_LD * LAPACKE_TEST_LD];
    const lapack_int top = rows - kl - ku - 1;
    const size_t len = lapacke_test_alloc_len(layout, rows, n, ld);
    size_t p;
    lapack_int i, j;

    lapacke_test_dfill_latms(LAPACK_COL_MAJOR, m, n, tmp, LAPACKE_TEST_LD, 'N',
                             kl, ku, seed);
    for (p = 0; p < len; p++) {
        if (lapacke_test_map_position(layout, rows, n, ld, p, &i, &j) &&
            i >= top && lapacke_test_region_band(i - top, j, m, ku)) {
            ab[p] = tmp[(i - top + j - ku) + j * LAPACKE_TEST_LD];
        } else {
            ab[p] = 0.25; /* never referenced */
        }
    }
}

/* Scatter the uplo band of the n-by-n column-major matrix tmp (leading
 * dimension LAPACKE_TEST_LD) into the (kd + 1)-by-n band array ab. */
static void lapacke_test_dscatter_band(int layout, char uplo, lapack_int n,
                                       lapack_int kd, const double *tmp,
                                       double *ab, lapack_int ld)
{
    const int upper = (uplo == 'U' || uplo == 'u');
    const lapack_int shift = upper ? kd : 0;
    const size_t len = lapacke_test_alloc_len(layout, kd + 1, n, ld);
    size_t p;
    lapack_int i, j;

    for (p = 0; p < len; p++) {
        if (lapacke_test_map_position(layout, kd + 1, n, ld, p, &i, &j) &&
            lapacke_test_region_band(i, j, n, shift)) {
            ab[p] = tmp[(i + j - shift) + j * LAPACKE_TEST_LD];
        } else {
            ab[p] = 0.25; /* never referenced */
        }
    }
}

/**
 * \brief Fill ab with a random symmetric positive definite n-by-n band
 * matrix in the uplo band storage, double precision.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  uplo   'U' or 'L'.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[in]  kd     Number of super/subdiagonals.
 * \param[out] ab     Buffer of lapacke_test_alloc_len(layout, kd + 1, n, ld)
 *                    elements.
 * \param[in]  ld     Leading dimension of ab in the given layout.
 */
void lapacke_test_dfill_pb(int layout, char uplo, lapack_int n, lapack_int kd,
                           double *ab, lapack_int ld)
{
    static const lapack_int seed[4] = {2040, 2041, 2042, 2043};
    double tmp[LAPACKE_TEST_LD * LAPACKE_TEST_LD];
    lapacke_test_dfill_latms(LAPACK_COL_MAJOR, n, n, tmp, LAPACKE_TEST_LD, 'P',
                             kd, kd, seed);
    lapacke_test_dscatter_band(layout, uplo, n, kd, tmp, ab, ld);
}

/**
 * \brief Fill ab with a random nonsingular triangular n-by-n band matrix in
 * the uplo band storage, double precision.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  uplo   'U' or 'L'.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[in]  kd     Number of super/subdiagonals.
 * \param[out] ab     Buffer of lapacke_test_alloc_len(layout, kd + 1, n, ld)
 *                    elements.
 * \param[in]  ld     Leading dimension of ab in the given layout.
 */
void lapacke_test_dfill_tb(int layout, char uplo, lapack_int n, lapack_int kd,
                           double *ab, lapack_int ld)
{
    static const lapack_int seed[4] = {2050, 2051, 2052, 2053};
    const int upper = (uplo == 'U' || uplo == 'u');
    double tmp[LAPACKE_TEST_LD * LAPACKE_TEST_LD];
    lapacke_test_dfill_latms(LAPACK_COL_MAJOR, n, n, tmp, LAPACKE_TEST_LD, 'N',
                             upper ? 0 : kd, upper ? kd : 0, seed);
    lapacke_test_dscatter_band(layout, uplo, n, kd, tmp, ab, ld);
}

/* Pack the uplo triangle of the n-by-n column-major matrix tmp (leading
 * dimension LAPACKE_TEST_LD) into ap in the packed storage of the layout,
 * as LAPACKE_?tp_trans lays it out. */
static void lapacke_test_dpack(int layout, char uplo, lapack_int n,
                               const double *tmp, double *ap)
{
    const int colmaj = (layout == LAPACK_COL_MAJOR);
    const int upper = (uplo == 'U' || uplo == 'u');
    lapack_int i, j;

    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            size_t p;
            if (upper ? i > j : i < j) {
                continue;
            }
            if (colmaj == upper) {
                p = colmaj ? i + (size_t)j * (j + 1) / 2
                           : j + (size_t)i * (i + 1) / 2;
            } else {
                p = colmaj ? (i - j) + (size_t)j * (2 * n - j + 1) / 2
                           : (j - i) + (size_t)i * (2 * n - i + 1) / 2;
            }
            ap[p] = tmp[i + j * LAPACKE_TEST_LD];
        }
    }
}

/**
 * \brief Fill ap with the uplo triangle of a random symmetric positive
 * definite n-by-n matrix in packed storage, double precision.
 *
 * The same matrix as lapacke_test_dfill_spd, packed.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  uplo   'U' or 'L'.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[out] ap     Buffer of n * (n + 1) / 2 elements.
 */
void lapacke_test_dfill_pp(int layout, char uplo, lapack_int n, double *ap)
{
    double tmp[LAPACKE_TEST_LD * LAPACKE_TEST_LD];
    lapacke_test_dfill_spd(LAPACK_COL_MAJOR, n, tmp, LAPACKE_TEST_LD);
    lapacke_test_dpack(layout, uplo, n, tmp, ap);
}

/**
 * \brief Fill ap with the uplo triangle of a random symmetric indefinite
 * n-by-n matrix in packed storage, double precision.
 *
 * The same matrix as lapacke_test_dfill_sym, packed.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  uplo   'U' or 'L'.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[out] ap     Buffer of n * (n + 1) / 2 elements.
 */
void lapacke_test_dfill_sp(int layout, char uplo, lapack_int n, double *ap)
{
    double tmp[LAPACKE_TEST_LD * LAPACKE_TEST_LD];
    lapacke_test_dfill_sym(LAPACK_COL_MAJOR, n, tmp, LAPACKE_TEST_LD);
    lapacke_test_dpack(layout, uplo, n, tmp, ap);
}

/**
 * \brief Fill ap with a random nonsingular uplo triangular n-by-n matrix in
 * packed storage, double precision.
 *
 * The same matrix as lapacke_test_dfill_tri, packed.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  uplo   'U' or 'L'.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[out] ap     Buffer of n * (n + 1) / 2 elements.
 */
void lapacke_test_dfill_tp(int layout, char uplo, lapack_int n, double *ap)
{
    double tmp[LAPACKE_TEST_LD * LAPACKE_TEST_LD];
    lapacke_test_dfill_tri(LAPACK_COL_MAJOR, uplo, n, tmp, LAPACKE_TEST_LD);
    lapacke_test_dpack(layout, uplo, n, tmp, ap);
}

/**
 * \brief Fill v with len deterministic nonzero values in (-0.5, 0.5),
 * double precision.
 *
 * For vectors whose content only has to be finite and nonzero: reflector
 * scalars, tridiagonal off-diagonals, right-hand sides.
 *
 * \param[in]  len Number of entries.
 * \param[out] v   The vector.
 */
void lapacke_test_dfill_vec(lapack_int len, double *v)
{
    lapack_int i;
    for (i = 0; i < len; i++) {
        const double re = (double)((i * 7 + 3) % 11) / 11.0 - 0.5;
        v[i] = re;
    }
}

/**
 * \brief Fill v with len deterministic values in [2, 3), double precision.
 *
 * For vectors that must be positive: the diagonal of a positive definite
 * tridiagonal matrix, scale factors, norms.
 *
 * \param[in]  len Number of entries.
 * \param[out] v   The vector.
 */
void lapacke_test_dfill_pos(lapack_int len, double *v)
{
    lapack_int i;
    for (i = 0; i < len; i++) {
        const double re = 2.0 + (double)((i * 7 + 3) % 11) / 11.0;
        v[i] = re;
    }
}
