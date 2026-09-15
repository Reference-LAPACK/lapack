#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zppcon. */
#define LAPACKE_ZPPCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_pp(layout, 'U', N, ap);                             \
        anorm[0] = 1.0;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zppcon)(layout, 'U', N, ap, anorm[0], &rcond),  \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(zppcon)
{
    lapack_complex_double ap[LD * LD];
    double anorm[1];
    double rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "zppcon anorm uplo=U", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -5,
            (lapacke_test_zfill_pp(layout, 'U', N, ap), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zppcon)(layout, 'U', N, ap, anorm[0], &rcond));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppcon ap uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4, (lapacke_test_zfill_pp(layout, 'U', N, ap), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zppcon)(layout, 'U', N, ap, anorm[0], &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "zppcon anorm uplo=L", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -5,
            (lapacke_test_zfill_pp(layout, 'L', N, ap), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zppcon)(layout, 'L', N, ap, anorm[0], &rcond));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppcon ap uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4, (lapacke_test_zfill_pp(layout, 'L', N, ap), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zppcon)(layout, 'L', N, ap, anorm[0], &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        anorm[0] = 1.0;
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_zfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_check("zppcon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_zppcon)(layout, 'U', N, ap,
                                                      anorm[0], &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZPPCON_ALLOC_TEST(0, 0, "zppcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPPCON_ALLOC_TEST(0, 1, "zppcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPPCON_ALLOC_TEST(0, 2, "zppcon allocation count", 0);
    lapacke_test_check_alloc_count("zppcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPPCON_ALLOC_TEST(1, 0, "zppcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPPCON_ALLOC_TEST(1, 1, "zppcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPPCON_ALLOC_TEST(1, 2, "zppcon transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPPCON_ALLOC_TEST(1, 3, "zppcon allocation count", 0);
    lapacke_test_check_alloc_count("zppcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPPCON_ALLOC_TEST(2, 0, "zppcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zppcon invalid layout allocation count");
}
