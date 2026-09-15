#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call spocon. */
#define LAPACKE_SPOCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_spd(layout, N, a, LD);                              \
        anorm[0] = 1.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_spocon)(layout, 'U', N, a, LD,   \
                                                      anorm[0], &rcond),       \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(spocon)
{
    float a[LD * LD];
    float anorm[1];
    float rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "spocon a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -4,
            (lapacke_test_sfill_spd(layout, N, a, LD), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_spocon)(layout, 'U', N, a, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "spocon anorm uplo=U", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_sfill_spd(layout, N, a, LD), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_spocon)(layout, 'U', N, a, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "spocon a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -4,
            (lapacke_test_sfill_spd(layout, N, a, LD), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_spocon)(layout, 'L', N, a, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "spocon anorm uplo=L", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_sfill_spd(layout, N, a, LD), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_spocon)(layout, 'L', N, a, LD, anorm[0],
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        anorm[0] = 1.0f;
        lapacke_test_sfill_nan(layout, N, N, a, LD);
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_check("spocon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_spocon)(layout, 'U', N, a, LD,
                                                      anorm[0], &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SPOCON_ALLOC_TEST(0, 0, "spocon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPOCON_ALLOC_TEST(0, 1, "spocon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPOCON_ALLOC_TEST(0, 2, "spocon allocation count", 0);
    lapacke_test_check_alloc_count("spocon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPOCON_ALLOC_TEST(1, 0, "spocon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPOCON_ALLOC_TEST(1, 1, "spocon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPOCON_ALLOC_TEST(1, 2, "spocon transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPOCON_ALLOC_TEST(1, 3, "spocon allocation count", 0);
    lapacke_test_check_alloc_count("spocon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPOCON_ALLOC_TEST(2, 0, "spocon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("spocon invalid layout allocation count");
}
