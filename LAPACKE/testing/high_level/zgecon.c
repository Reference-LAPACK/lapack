#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zgecon. */
#define LAPACKE_ZGECON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, N, N, a, LD);                               \
        anorm[0] = 1.0;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zgecon)(layout, '1', N, a, LD,   \
                                                      anorm[0], &rcond),       \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zgecon)
{
    lapack_complex_double a[LD * LD];
    double anorm[1];
    double rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgecon a", l, N, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_zfill(layout, N, N, a, LD), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zgecon)(layout, '1', N, a, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "zgecon anorm", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_zfill(layout, N, N, a, LD), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zgecon)(layout, '1', N, a, LD, anorm[0],
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        anorm[0] = 1.0;
        lapacke_test_zfill_nan(layout, N, N, a, LD);
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_check("zgecon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_zgecon)(layout, '1', N, a, LD,
                                                      anorm[0], &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZGECON_ALLOC_TEST(0, 0, "zgecon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGECON_ALLOC_TEST(0, 1, "zgecon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGECON_ALLOC_TEST(0, 2, "zgecon allocation count", 0);
    lapacke_test_check_alloc_count("zgecon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZGECON_ALLOC_TEST(1, 0, "zgecon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGECON_ALLOC_TEST(1, 1, "zgecon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGECON_ALLOC_TEST(1, 2, "zgecon transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGECON_ALLOC_TEST(1, 3, "zgecon allocation count", 0);
    lapacke_test_check_alloc_count("zgecon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZGECON_ALLOC_TEST(2, 0, "zgecon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgecon invalid layout allocation count");
}
