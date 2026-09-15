#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zptcon. */
#define LAPACKE_ZPTCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_dfill_pos(LD * LD, d);                                    \
        lapacke_test_zfill_vec(LD * LD, e);                                    \
        anorm[0] = 1.0;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zptcon)(N, d, e, anorm[0], &rcond), expected);  \
    } while (0)

LAPACKE_TEST(zptcon)
{
    double d[LD * LD];
    lapack_complex_double e[LD * LD];
    double anorm[1];
    double rcond;

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "zptcon anorm", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -4,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zptcon)(N, d, e, anorm[0], &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "zptcon d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -2,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zptcon)(N, d, e, anorm[0], &rcond));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptcon e", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -3,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zptcon)(N, d, e, anorm[0], &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        anorm[0] = 1.0;
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_dfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_zfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check(
            "zptcon NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zptcon)(N, d, e, anorm[0], &rcond) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZPTCON_ALLOC_TEST(0, 0, "zptcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPTCON_ALLOC_TEST(0, 1, "zptcon allocation count", 0);
    lapacke_test_check_alloc_count("zptcon col-major allocation count");
}
