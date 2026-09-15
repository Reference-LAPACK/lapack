#include <string.h>

#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Self-test for the matrix generation helpers: the routine tests rely on
 * the fills being deterministic and on the structural properties promised
 * by each fill type (also for the fills no routine test uses yet). */

/* Buffer index of logical element (i, j). */
static size_t idx(int layout, lapack_int i, lapack_int j)
{
    return layout == LAPACK_COL_MAJOR ? i + j * (size_t)LD : i * (size_t)LD + j;
}

/* The fills are deterministic, so identical runs must agree bit for bit. */
static int equal_buffers(const float *a, const float *b, size_t len)
{
    return memcmp(a, b, len * sizeof(*a)) == 0;
}

static int is_zero(float v)
{
    return v == 0.0f;
}

/* a[i][j] == conj(a[j][i]); plain equality in the real precisions. */
static int conj_equal(float x, float y)
{
    return x == y;
}

static int is_triangular(const float *a, int layout, char uplo)
{
    lapack_int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            const float v = a[idx(layout, i, j)];
            const int outside = (uplo == 'U') ? (i > j) : (i < j);
            if (outside && !is_zero(v)) return 0;
            if (i == j && is_zero(v)) return 0;
        }
    }
    return 1;
}

LAPACKE_TEST(smatgen)
{
    float a[LD * LD], b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];
        const char *lname = lapacke_test_layout_names[l];
        const size_t len_sq = lapacke_test_alloc_len(layout, N, N, LD);
        lapack_int i, j;
        int ok;

        /* The seeds are fixed, so repeated fills of every type must
         * produce bit-identical buffers. */
        ok = 1;
        lapacke_test_sfill(layout, M, N, a, LD);
        lapacke_test_sfill(layout, M, N, b, LD);
        ok = ok &&
             equal_buffers(a, b, lapacke_test_alloc_len(layout, M, N, LD));
        lapacke_test_sfill_rhs(layout, M, NRHS, a, LD);
        lapacke_test_sfill_rhs(layout, M, NRHS, b, LD);
        ok = ok &&
             equal_buffers(a, b, lapacke_test_alloc_len(layout, M, NRHS, LD));
        lapacke_test_sfill_spd(layout, N, a, LD);
        lapacke_test_sfill_spd(layout, N, b, LD);
        ok = ok && equal_buffers(a, b, len_sq);
        lapacke_test_sfill_sym(layout, N, a, LD);
        lapacke_test_sfill_sym(layout, N, b, LD);
        ok = ok && equal_buffers(a, b, len_sq);
        lapacke_test_sfill_tri(layout, 'U', N, a, LD);
        lapacke_test_sfill_tri(layout, 'U', N, b, LD);
        ok = ok && equal_buffers(a, b, len_sq);
        lapacke_test_check("smatgen fills deterministic", lname, ok ? 0 : -999,
                           0);

        lapacke_test_sfill_sym(layout, N, a, LD);
        ok = 1;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                ok = ok &&
                     conj_equal(a[idx(layout, i, j)], a[idx(layout, j, i)]);
            }
        }
        lapacke_test_check("smatgen sym is symmetric", lname, ok ? 0 : -999, 0);

        lapacke_test_sfill_tri(layout, 'U', N, a, LD);
        ok = is_triangular(a, layout, 'U');
        lapacke_test_sfill_tri(layout, 'L', N, a, LD);
        ok = ok && is_triangular(a, layout, 'L');
        lapacke_test_check("smatgen tri is triangular", lname, ok ? 0 : -999,
                           0);
    }

    /* The SPD fill must Cholesky-factorize; the indefinite fill must be
     * rejected with a positive info (not positive definite). */
    lapacke_test_sfill_spd(LAPACK_COL_MAJOR, N, a, LD);
    lapacke_test_check(
        "smatgen spd is positive definite", NULL,
        API_SUFFIX(LAPACKE_spotrf)(LAPACK_COL_MAJOR, 'U', N, a, LD), 0);

    lapacke_test_sfill_sym(LAPACK_COL_MAJOR, N, a, LD);
    lapacke_test_check(
        "smatgen sym is indefinite", NULL,
        API_SUFFIX(LAPACKE_spotrf)(LAPACK_COL_MAJOR, 'U', N, a, LD) > 0 ? 0
                                                                        : -999,
        0);

    /* The band, packed and vector fills: deterministic, and structured as
     * promised (positive definite band and packed matrices factorize,
     * triangular band and packed matrices are nonsingular). */
    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];
        const char *lname = lapacke_test_layout_names[l];
        float ab[LD * LD], ab2[LD * LD], ap[LD * (LD + 1) / 2], v[LD];
        lapack_int i;
        int ok;

        ok = 1;
        lapacke_test_sfill_gb(layout, M, N, 1, 1, 4, ab, LD);
        lapacke_test_sfill_gb(layout, M, N, 1, 1, 4, ab2, LD);
        ok = ok &&
             equal_buffers(ab, ab2, lapacke_test_alloc_len(layout, 4, N, LD));
        lapacke_test_sfill_pp(layout, 'U', N, ap);
        lapacke_test_sfill_pp(layout, 'U', N, ab2);
        ok = ok && equal_buffers(ap, ab2, (size_t)N * (N + 1) / 2);
        lapacke_test_check("smatgen band and packed fills deterministic", lname,
                           ok ? 0 : -999, 0);

        ok = 1;
        lapacke_test_sfill_pb(layout, 'U', N, 1, ab, LD);
        ok = ok && API_SUFFIX(LAPACKE_spbtrf)(layout, 'U', N, 1, ab, LD) == 0;
        lapacke_test_sfill_pb(layout, 'L', N, 1, ab, LD);
        ok = ok && API_SUFFIX(LAPACKE_spbtrf)(layout, 'L', N, 1, ab, LD) == 0;
        lapacke_test_sfill_pp(layout, 'U', N, ap);
        ok = ok && API_SUFFIX(LAPACKE_spptrf)(layout, 'U', N, ap) == 0;
        lapacke_test_sfill_pp(layout, 'L', N, ap);
        ok = ok && API_SUFFIX(LAPACKE_spptrf)(layout, 'L', N, ap) == 0;
        lapacke_test_check("smatgen pb and pp are positive definite", lname,
                           ok ? 0 : -999, 0);

        ok = 1;
        lapacke_test_sfill_tb(layout, 'U', N, 1, ab, LD);
        lapacke_test_sfill_rhs(layout, N, NRHS, ab2, LD);
        ok = ok && API_SUFFIX(LAPACKE_stbtrs)(layout, 'U', 'N', 'N', N, 1, NRHS,
                                              ab, LD, ab2, LD) == 0;
        lapacke_test_sfill_tb(layout, 'L', N, 1, ab, LD);
        ok = ok && API_SUFFIX(LAPACKE_stbtrs)(layout, 'L', 'N', 'N', N, 1, NRHS,
                                              ab, LD, ab2, LD) == 0;
        lapacke_test_sfill_tp(layout, 'U', N, ap);
        ok = ok && API_SUFFIX(LAPACKE_stptrs)(layout, 'U', 'N', 'N', N, NRHS,
                                              ap, ab2, LD) == 0;
        lapacke_test_sfill_tp(layout, 'L', N, ap);
        ok = ok && API_SUFFIX(LAPACKE_stptrs)(layout, 'L', 'N', 'N', N, NRHS,
                                              ap, ab2, LD) == 0;
        lapacke_test_check("smatgen tb and tp are nonsingular", lname,
                           ok ? 0 : -999, 0);

        ok = 1;
        lapacke_test_sfill_vec(LD, v);
        for (i = 0; i < LD; i++) {
            ok = ok && !is_zero(v[i]);
        }
        lapacke_test_sfill_pos(LD, v);
        for (i = 0; i < LD; i++) {
            ok = ok && (v[i]) >= 2;
        }
        lapacke_test_check("smatgen vec is nonzero and pos is positive", lname,
                           ok ? 0 : -999, 0);
    }
}
