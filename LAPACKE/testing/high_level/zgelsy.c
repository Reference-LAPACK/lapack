#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zgelsy. */
#define LAPACKE_ZGELSY_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, M, N, a, LD);                               \
        lapacke_test_zfill_rhs(layout, M, NRHS, b, LD);                        \
        lapacke_test_fill_int(LD * LD, jpvt, 0);                               \
        rcond[0] = 0.0;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zgelsy)(layout, M, N, NRHS, a,   \
                                                      LD, b, LD, jpvt,         \
                                                      rcond[0], &rank),        \
                           expected);                                          \
    } while (0)

/* B has max(M, N) rows, all of them documented inputs for M > N. For
 * M < N only the first M rows are, while LAPACKE checks all max(M, N):
 * a divergence this fixture does not reach. */
LAPACKE_TEST(zgelsy)
{
    lapack_complex_double a[LD * LD];
    lapack_complex_double b[LD * LD];
    lapack_int jpvt[LD * LD];
    double rcond[1];
    lapack_int rank;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgelsy a", l, M, N, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_zfill(layout, M, N, a, LD),
             lapacke_test_zfill_rhs(layout, M, NRHS, b, LD),
             lapacke_test_fill_int(LD * LD, jpvt, 0), (rcond[0] = 0.0)),
            API_SUFFIX(LAPACKE_zgelsy)(layout, M, N, NRHS, a, LD, b, LD, jpvt,
                                       rcond[0], &rank));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgelsy b", l, M, NRHS, b, LD, lapacke_test_region_full, -7,
            (lapacke_test_zfill(layout, M, N, a, LD),
             lapacke_test_zfill_rhs(layout, M, NRHS, b, LD),
             lapacke_test_fill_int(LD * LD, jpvt, 0), (rcond[0] = 0.0)),
            API_SUFFIX(LAPACKE_zgelsy)(layout, M, N, NRHS, a, LD, b, LD, jpvt,
                                       rcond[0], &rank));

        LAPACKE_TEST_DNAN_SWEEP(
            "zgelsy rcond", l, 1, 1, rcond, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -10,
            (lapacke_test_zfill(layout, M, N, a, LD),
             lapacke_test_zfill_rhs(layout, M, NRHS, b, LD),
             lapacke_test_fill_int(LD * LD, jpvt, 0), (rcond[0] = 0.0)),
            API_SUFFIX(LAPACKE_zgelsy)(layout, M, N, NRHS, a, LD, b, LD, jpvt,
                                       rcond[0], &rank));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_int(LD * LD, jpvt, 0);
        lapacke_test_zfill_nan(layout, M, N, a, LD);
        lapacke_test_zfill_nan(layout, M, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, 1, 1, rcond,
                               LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_check(
            "zgelsy NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgelsy)(layout, M, N, NRHS, a, LD, b, LD, jpvt,
                                       rcond[0], &rank) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZGELSY_ALLOC_TEST(0, 0, "zgelsy work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGELSY_ALLOC_TEST(0, 1, "zgelsy work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGELSY_ALLOC_TEST(0, 2, "zgelsy allocation count", 0);
    lapacke_test_check_alloc_count("zgelsy col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZGELSY_ALLOC_TEST(1, 0, "zgelsy work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGELSY_ALLOC_TEST(1, 1, "zgelsy work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGELSY_ALLOC_TEST(1, 2, "zgelsy transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGELSY_ALLOC_TEST(1, 3, "zgelsy transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGELSY_ALLOC_TEST(1, 4, "zgelsy allocation count", 0);
    lapacke_test_check_alloc_count("zgelsy row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZGELSY_ALLOC_TEST(2, 0, "zgelsy invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgelsy invalid layout allocation count");
}
