#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cgtcon. */
#define LAPACKE_CGTCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_cfill_vec(LD * LD, dl);                                   \
        lapacke_test_cfill_pos(LD * LD, d);                                    \
        lapacke_test_cfill_vec(LD * LD, du);                                   \
        lapacke_test_cfill_vec(LD * LD, du2);                                  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        anorm[0] = 1.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cgtcon)('1', N, dl, d, du, du2,  \
                                                      ipiv, anorm[0], &rcond), \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cgtcon)
{
    lapack_complex_float dl[LD * LD];
    lapack_complex_float d[LD * LD];
    lapack_complex_float du[LD * LD];
    lapack_complex_float du2[LD * LD];
    lapack_int ipiv[LD * LD];
    float anorm[1];
    float rcond;

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "cgtcon anorm", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -8,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtcon d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -4,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtcon dl", l, N - 1, 1, dl, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -3,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtcon du", l, N - 1, 1, du, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -5,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtcon du2", l, N - 2, 1, du2, LAPACKE_TEST_VLD(layout, N - 2),
            lapacke_test_region_full, -6,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        anorm[0] = 1.0f;
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_cfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_cfill_nan(layout, N - 1, 1, dl,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_cfill_nan(layout, N - 1, 1, du,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_cfill_nan(layout, N - 2, 1, du2,
                               LAPACKE_TEST_VLD(layout, N - 2));
        lapacke_test_check(
            "cgtcon NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CGTCON_ALLOC_TEST(0, 0, "cgtcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGTCON_ALLOC_TEST(0, 1, "cgtcon allocation count", 0);
    lapacke_test_check_alloc_count("cgtcon col-major allocation count");
}
