#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cpocon. */
#define LAPACKE_CPOCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_spd(layout, N, a, LD);                              \
        anorm[0] = 1.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cpocon)(layout, 'U', N, a, LD,   \
                                                      anorm[0], &rcond),       \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cpocon)
{
    lapack_complex_float a[LD * LD];
    float anorm[1];
    float rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cpocon a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -4,
            (lapacke_test_cfill_spd(layout, N, a, LD), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cpocon)(layout, 'U', N, a, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "cpocon anorm uplo=U", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_cfill_spd(layout, N, a, LD), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cpocon)(layout, 'U', N, a, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpocon a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -4,
            (lapacke_test_cfill_spd(layout, N, a, LD), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cpocon)(layout, 'L', N, a, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "cpocon anorm uplo=L", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_cfill_spd(layout, N, a, LD), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cpocon)(layout, 'L', N, a, LD, anorm[0],
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        anorm[0] = 1.0f;
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_check("cpocon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_cpocon)(layout, 'U', N, a, LD,
                                                      anorm[0], &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CPOCON_ALLOC_TEST(0, 0, "cpocon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPOCON_ALLOC_TEST(0, 1, "cpocon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPOCON_ALLOC_TEST(0, 2, "cpocon allocation count", 0);
    lapacke_test_check_alloc_count("cpocon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CPOCON_ALLOC_TEST(1, 0, "cpocon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPOCON_ALLOC_TEST(1, 1, "cpocon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPOCON_ALLOC_TEST(1, 2, "cpocon transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPOCON_ALLOC_TEST(1, 3, "cpocon allocation count", 0);
    lapacke_test_check_alloc_count("cpocon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CPOCON_ALLOC_TEST(2, 0, "cpocon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cpocon invalid layout allocation count");
}
