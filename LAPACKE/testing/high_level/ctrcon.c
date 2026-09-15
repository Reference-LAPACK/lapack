#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call ctrcon. */
#define LAPACKE_CTRCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_tri(layout, 'U', N, a, LD);                         \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_ctrcon)(layout, '1', 'U', 'N',   \
                                                      N, a, LD, &rcond),       \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(ctrcon)
{
    lapack_complex_float a[LD * LD];
    float rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP("ctrcon a uplo=U diag=N", l, N, N, a, LD,
                                lapacke_test_region_upper, -6,
                                (lapacke_test_cfill_tri(layout, 'U', N, a, LD)),
                                API_SUFFIX(LAPACKE_ctrcon)(
                                    layout, '1', 'U', 'N', N, a, LD, &rcond));

        LAPACKE_TEST_CNAN_SWEEP("ctrcon a uplo=U diag=U", l, N, N, a, LD,
                                lapacke_test_region_strict_upper, -6,
                                (lapacke_test_cfill_tri(layout, 'U', N, a, LD)),
                                API_SUFFIX(LAPACKE_ctrcon)(
                                    layout, '1', 'U', 'U', N, a, LD, &rcond));

        LAPACKE_TEST_CNAN_SWEEP("ctrcon a uplo=L diag=N", l, N, N, a, LD,
                                lapacke_test_region_lower, -6,
                                (lapacke_test_cfill_tri(layout, 'L', N, a, LD)),
                                API_SUFFIX(LAPACKE_ctrcon)(
                                    layout, '1', 'L', 'N', N, a, LD, &rcond));

        LAPACKE_TEST_CNAN_SWEEP("ctrcon a uplo=L diag=U", l, N, N, a, LD,
                                lapacke_test_region_strict_lower, -6,
                                (lapacke_test_cfill_tri(layout, 'L', N, a, LD)),
                                API_SUFFIX(LAPACKE_ctrcon)(
                                    layout, '1', 'L', 'U', N, a, LD, &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        lapacke_test_check("ctrcon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_ctrcon)(layout, '1', 'U', 'N', N,
                                                      a, LD, &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CTRCON_ALLOC_TEST(0, 0, "ctrcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTRCON_ALLOC_TEST(0, 1, "ctrcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTRCON_ALLOC_TEST(0, 2, "ctrcon allocation count", 0);
    lapacke_test_check_alloc_count("ctrcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CTRCON_ALLOC_TEST(1, 0, "ctrcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTRCON_ALLOC_TEST(1, 1, "ctrcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTRCON_ALLOC_TEST(1, 2, "ctrcon transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CTRCON_ALLOC_TEST(1, 3, "ctrcon allocation count", 0);
    lapacke_test_check_alloc_count("ctrcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CTRCON_ALLOC_TEST(2, 0, "ctrcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("ctrcon invalid layout allocation count");
}
