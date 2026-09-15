#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dgels. */
#define LAPACKE_DGELS_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, M, N, a, LD);                               \
        lapacke_test_dfill_rhs(layout, M, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dgels)(layout, 'N', M, N, NRHS, a, LD, b, LD),  \
            expected);                                                         \
    } while (0)

/* B has max(M, N) rows, all of them documented inputs for M > N. For
 * M < N only the first M rows are, while LAPACKE checks all max(M, N):
 * a divergence this fixture does not reach. */
LAPACKE_TEST(dgels)
{
    double a[LD * LD];
    double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgels a", l, M, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_dfill(layout, M, N, a, LD),
             lapacke_test_dfill_rhs(layout, M, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgels)(layout, 'N', M, N, NRHS, a, LD, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgels b", l, M, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_dfill(layout, M, N, a, LD),
             lapacke_test_dfill_rhs(layout, M, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgels)(layout, 'N', M, N, NRHS, a, LD, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_nan(layout, M, N, a, LD);
        lapacke_test_dfill_nan(layout, M, NRHS, b, LD);
        lapacke_test_check("dgels NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_dgels)(layout, 'N', M, N, NRHS, a,
                                                     LD, b, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DGELS_ALLOC_TEST(0, 0, "dgels work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGELS_ALLOC_TEST(0, 1, "dgels allocation count", 0);
    lapacke_test_check_alloc_count("dgels col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGELS_ALLOC_TEST(1, 0, "dgels work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGELS_ALLOC_TEST(1, 1, "dgels transpose alloc failure (a_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGELS_ALLOC_TEST(1, 2, "dgels transpose alloc failure (b_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGELS_ALLOC_TEST(1, 3, "dgels allocation count", 0);
    lapacke_test_check_alloc_count("dgels row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGELS_ALLOC_TEST(2, 0, "dgels invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgels invalid layout allocation count");
}
