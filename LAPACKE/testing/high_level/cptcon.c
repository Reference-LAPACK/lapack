#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cptcon. */
#define LAPACKE_CPTCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_sfill_pos(LD * LD, d);                                    \
        lapacke_test_cfill_vec(LD * LD, e);                                    \
        anorm[0] = 1.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_cptcon)(N, d, e, anorm[0], &rcond), expected);  \
    } while (0)

LAPACKE_TEST(cptcon)
{
    float d[LD * LD];
    lapack_complex_float e[LD * LD];
    float anorm[1];
    float rcond;

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "cptcon anorm", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -4,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, e), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cptcon)(N, d, e, anorm[0], &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "cptcon d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -2,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, e), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cptcon)(N, d, e, anorm[0], &rcond));

        LAPACKE_TEST_CNAN_SWEEP(
            "cptcon e", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -3,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, e), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cptcon)(N, d, e, anorm[0], &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        anorm[0] = 1.0f;
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_sfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_cfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check(
            "cptcon NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cptcon)(N, d, e, anorm[0], &rcond) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CPTCON_ALLOC_TEST(0, 0, "cptcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPTCON_ALLOC_TEST(0, 1, "cptcon allocation count", 0);
    lapacke_test_check_alloc_count("cptcon col-major allocation count");
}
