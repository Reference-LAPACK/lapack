#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define K 2
#define NB 2
#define KBIG LAPACKE_TEST_N /* more reflectors than NB */

/* One call: refill V (nrows_v-by-K), T (NB-by-K) and C (rows-by-cols),
 * schedule the malloc failure countdown and check the info result of a
 * dgemqrt call applying Q from the given side in the indexed layout. */
#define LAPACKE_DGEMQRT_TEST(layout_index, side, rows, cols, countdown, name,  \
                             expected)                                         \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        const lapack_int nrows_v = (side) == 'L' ? (rows) : (cols);            \
        lapacke_test_dfill(layout, nrows_v, K, v, LD);                         \
        lapacke_test_dfill(layout, NB, K, t, LD);                              \
        lapacke_test_dfill(layout, rows, cols, c, LD);                         \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dgemqrt)(                        \
                               layout, side, 'N', rows, cols, K, NB, v, LD, t, \
                               LD, c, LD),                                     \
                           expected);                                          \
    } while (0)

/* Apply Q with KBIG > NB reflectors to the same V, T and C in both layouts
 * and return the largest difference between the two results, or -1 if a
 * call fails. V is nrows_v-by-KBIG and T is NB-by-KBIG in both layouts. */
static double lapacke_test_dgemqrt_layout_diff(char side, char trans,
                                               lapack_int rows,
                                               lapack_int cols)
{
    double v[2][LD * LD], t[2][LD * LD], c[2][LD * LD];
    const lapack_int nrows_v = side == 'L' ? rows : cols;
    double diff = 0.0;
    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];
        lapacke_test_dfill(layout, nrows_v, KBIG, v[l], LD);
        lapacke_test_dfill(layout, NB, KBIG, t[l], LD);
        lapacke_test_dfill(layout, rows, cols, c[l], LD);
        if (API_SUFFIX(LAPACKE_dgemqrt)(layout, side, trans, rows, cols, KBIG,
                                        NB, v[l], LD, t[l], LD, c[l],
                                        LD) != 0)
            return -1.0;
    }
    for (lapack_int i = 0; i < rows; i++)
        for (lapack_int j = 0; j < cols; j++) {
            double d = c[0][i + j * LD] - c[1][i * LD + j];
            d = d < 0 ? -d : d;
            if (d > diff)
                diff = d;
        }
    return diff;
}

LAPACKE_TEST(dgemqrt)
{
    double v[LD * LD], t[LD * LD], c[LD * LD];

    /* DGEMQRT needs a workspace of N*NB entries for side = 'L' and M*NB for
     * side = 'R', which LAPACKE allocates; C is N-by-M (more columns than
     * rows) for side = 'L' and M-by-N for side = 'R'. The test allocator
     * fails the test when a call writes past the end of a LAPACKE
     * allocation. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_DGEMQRT_TEST(l, 'L', N, M, -1, "dgemqrt side = L", 0);
        LAPACKE_DGEMQRT_TEST(l, 'R', M, N, -1, "dgemqrt side = R", 0);
    }

    /* In row-major layout V and T are copied as nrows_v-by-K and NB-by-K
     * matrices; with K > NB both layouts must give the same result. */
    for (size_t sd = 0; sd < 2; sd++) {
        for (size_t tt = 0; tt < 2; tt++) {
            const char side = sd ? 'R' : 'L', trans = tt ? 'T' : 'N';
            const double diff = lapacke_test_dgemqrt_layout_diff(
                side, trans, side == 'L' ? N : M, side == 'L' ? M : N);
            lapacke_test_check("dgemqrt row-major result (K > NB)",
                               side == 'L' ? "side L" : "side R",
                               diff >= 0.0 && diff <= 1e-12 ? 0 : 1, 0);
        }
    }

    /* Column-major: the workspace is the only allocation. */
    LAPACKE_DGEMQRT_TEST(0, 'L', N, M, 0, "dgemqrt work alloc failure",
                         LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_DGEMQRT_TEST(0, 'L', N, M, 1, "dgemqrt allocation count", 0);
    lapacke_test_check_alloc_count("dgemqrt col-major allocation count");

    /* Row-major allocates the workspace, then the transposed copies of V, T
     * and C. */
    LAPACKE_DGEMQRT_TEST(1, 'L', N, M, 0, "dgemqrt work alloc failure",
                         LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGEMQRT_TEST(1, 'L', N, M, 1, "dgemqrt transpose alloc failure (v)",
                         LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGEMQRT_TEST(1, 'L', N, M, 2, "dgemqrt transpose alloc failure (t)",
                         LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGEMQRT_TEST(1, 'L', N, M, 3, "dgemqrt transpose alloc failure (c)",
                         LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_DGEMQRT_TEST(1, 'L', N, M, 4, "dgemqrt allocation count", 0);
    lapacke_test_check_alloc_count("dgemqrt row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_DGEMQRT_TEST(2, 'L', N, M, 0, "dgemqrt invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgemqrt invalid layout allocation count");
}
