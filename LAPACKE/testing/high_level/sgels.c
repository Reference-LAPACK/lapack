#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sgels. */
#define LAPACKE_SGELS_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, M, N, a, LD);                               \
        lapacke_test_sfill_rhs(layout, M, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_sgels)(layout, 'N', M, N, NRHS, a, LD, b, LD),  \
            expected);                                                         \
    } while (0)

/* B has max(M, N) rows, all of them documented inputs for M > N. For
 * M < N only the first M rows are, while LAPACKE checks all max(M, N):
 * a divergence this fixture does not reach. */
LAPACKE_TEST(sgels)
{
    float a[LD * LD];
    float b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgels a", l, M, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_sfill_rhs(layout, M, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_sgels)(layout, 'N', M, N, NRHS, a, LD, b, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgels b", l, M, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_sfill_rhs(layout, M, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_sgels)(layout, 'N', M, N, NRHS, a, LD, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, M, N, a, LD);
        lapacke_test_sfill_nan(layout, M, NRHS, b, LD);
        lapacke_test_check("sgels NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_sgels)(layout, 'N', M, N, NRHS, a,
                                                     LD, b, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SGELS_ALLOC_TEST(0, 0, "sgels work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGELS_ALLOC_TEST(0, 1, "sgels allocation count", 0);
    lapacke_test_check_alloc_count("sgels col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SGELS_ALLOC_TEST(1, 0, "sgels work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGELS_ALLOC_TEST(1, 1, "sgels transpose alloc failure (a_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGELS_ALLOC_TEST(1, 2, "sgels transpose alloc failure (b_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGELS_ALLOC_TEST(1, 3, "sgels allocation count", 0);
    lapacke_test_check_alloc_count("sgels row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SGELS_ALLOC_TEST(2, 0, "sgels invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgels invalid layout allocation count");
}
