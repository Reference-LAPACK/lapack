#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zgttrs. */
#define LAPACKE_ZGTTRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_vec(LD * LD, dl);                                   \
        lapacke_test_zfill_pos(LD * LD, d);                                    \
        lapacke_test_zfill_vec(LD * LD, du);                                   \
        lapacke_test_zfill_vec(LD * LD, du2);                                  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zgttrs)(layout, 'N', N, NRHS,    \
                                                      dl, d, du, du2, ipiv, b, \
                                                      LD),                     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zgttrs)
{
    lapack_complex_double dl[LD * LD];
    lapack_complex_double d[LD * LD];
    lapack_complex_double du[LD * LD];
    lapack_complex_double du2[LD * LD];
    lapack_int ipiv[LD * LD];
    lapack_complex_double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgttrs b", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgttrs d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -6,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgttrs dl", l, N - 1, 1, dl, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -5,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgttrs du", l, N - 1, 1, du, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -7,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgttrs du2", l, N - 2, 1, du2, LAPACKE_TEST_VLD(layout, N - 2),
            lapacke_test_region_full, -8,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_zfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_zfill_nan(layout, N - 1, 1, dl,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_zfill_nan(layout, N - 1, 1, du,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_zfill_nan(layout, N - 2, 1, du2,
                               LAPACKE_TEST_VLD(layout, N - 2));
        lapacke_test_check(
            "zgttrs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZGTTRS_ALLOC_TEST(0, 0, "zgttrs allocation count", 0);
    lapacke_test_check_alloc_count("zgttrs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZGTTRS_ALLOC_TEST(1, 0, "zgttrs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGTTRS_ALLOC_TEST(1, 1, "zgttrs allocation count", 0);
    lapacke_test_check_alloc_count("zgttrs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZGTTRS_ALLOC_TEST(2, 0, "zgttrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgttrs invalid layout allocation count");
}
