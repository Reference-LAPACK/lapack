#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zgtcon. */
#define LAPACKE_ZGTCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_zfill_vec(LD * LD, dl);                                   \
        lapacke_test_zfill_pos(LD * LD, d);                                    \
        lapacke_test_zfill_vec(LD * LD, du);                                   \
        lapacke_test_zfill_vec(LD * LD, du2);                                  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        anorm[0] = 1.0;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zgtcon)('1', N, dl, d, du, du2,  \
                                                      ipiv, anorm[0], &rcond), \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zgtcon)
{
    lapack_complex_double dl[LD * LD];
    lapack_complex_double d[LD * LD];
    lapack_complex_double du[LD * LD];
    lapack_complex_double du2[LD * LD];
    lapack_int ipiv[LD * LD];
    double anorm[1];
    double rcond;

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "zgtcon anorm", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -8,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgtcon d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -4,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgtcon dl", l, N - 1, 1, dl, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -3,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgtcon du", l, N - 1, 1, du, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -5,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgtcon du2", l, N - 2, 1, du2, LAPACKE_TEST_VLD(layout, N - 2),
            lapacke_test_region_full, -6,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        anorm[0] = 1.0;
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_zfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_zfill_nan(layout, N - 1, 1, dl,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_zfill_nan(layout, N - 1, 1, du,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_zfill_nan(layout, N - 2, 1, du2,
                               LAPACKE_TEST_VLD(layout, N - 2));
        lapacke_test_check(
            "zgtcon NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgtcon)('1', N, dl, d, du, du2, ipiv, anorm[0],
                                       &rcond) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZGTCON_ALLOC_TEST(0, 0, "zgtcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGTCON_ALLOC_TEST(0, 1, "zgtcon allocation count", 0);
    lapacke_test_check_alloc_count("zgtcon col-major allocation count");
}
