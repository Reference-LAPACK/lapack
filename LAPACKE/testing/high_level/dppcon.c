#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dppcon. */
#define LAPACKE_DPPCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_pp(layout, 'U', N, ap);                             \
        anorm[0] = 1.0;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dppcon)(layout, 'U', N, ap, anorm[0], &rcond),  \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(dppcon)
{
    double ap[LD * LD];
    double anorm[1];
    double rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dppcon anorm uplo=U", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -5,
            (lapacke_test_dfill_pp(layout, 'U', N, ap), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_dppcon)(layout, 'U', N, ap, anorm[0], &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "dppcon ap uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4, (lapacke_test_dfill_pp(layout, 'U', N, ap), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_dppcon)(layout, 'U', N, ap, anorm[0], &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "dppcon anorm uplo=L", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -5,
            (lapacke_test_dfill_pp(layout, 'L', N, ap), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_dppcon)(layout, 'L', N, ap, anorm[0], &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "dppcon ap uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4, (lapacke_test_dfill_pp(layout, 'L', N, ap), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_dppcon)(layout, 'L', N, ap, anorm[0], &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        anorm[0] = 1.0;
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_dfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_check("dppcon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_dppcon)(layout, 'U', N, ap,
                                                      anorm[0], &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DPPCON_ALLOC_TEST(0, 0, "dppcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPPCON_ALLOC_TEST(0, 1, "dppcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPPCON_ALLOC_TEST(0, 2, "dppcon allocation count", 0);
    lapacke_test_check_alloc_count("dppcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DPPCON_ALLOC_TEST(1, 0, "dppcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPPCON_ALLOC_TEST(1, 1, "dppcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPPCON_ALLOC_TEST(1, 2, "dppcon transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPPCON_ALLOC_TEST(1, 3, "dppcon allocation count", 0);
    lapacke_test_check_alloc_count("dppcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DPPCON_ALLOC_TEST(2, 0, "dppcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dppcon invalid layout allocation count");
}
