#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dgttrs. */
#define LAPACKE_DGTTRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_vec(LD * LD, dl);                                   \
        lapacke_test_dfill_pos(LD * LD, d);                                    \
        lapacke_test_dfill_vec(LD * LD, du);                                   \
        lapacke_test_dfill_vec(LD * LD, du2);                                  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dgttrs)(layout, 'N', N, NRHS,    \
                                                      dl, d, du, du2, ipiv, b, \
                                                      LD),                     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dgttrs)
{
    double dl[LD * LD];
    double d[LD * LD];
    double du[LD * LD];
    double du2[LD * LD];
    lapack_int ipiv[LD * LD];
    double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgttrs b", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgttrs d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -6,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgttrs dl", l, N - 1, 1, dl, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -5,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgttrs du", l, N - 1, 1, du, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -7,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgttrs du2", l, N - 2, 1, du2, LAPACKE_TEST_VLD(layout, N - 2),
            lapacke_test_region_full, -8,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_dfill_nan(layout, N - 1, 1, dl,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_dfill_nan(layout, N - 1, 1, du,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_dfill_nan(layout, N - 2, 1, du2,
                               LAPACKE_TEST_VLD(layout, N - 2));
        lapacke_test_check(
            "dgttrs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgttrs)(layout, 'N', N, NRHS, dl, d, du, du2,
                                       ipiv, b, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_DGTTRS_ALLOC_TEST(0, 0, "dgttrs allocation count", 0);
    lapacke_test_check_alloc_count("dgttrs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGTTRS_ALLOC_TEST(1, 0, "dgttrs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGTTRS_ALLOC_TEST(1, 1, "dgttrs allocation count", 0);
    lapacke_test_check_alloc_count("dgttrs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGTTRS_ALLOC_TEST(2, 0, "dgttrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgttrs invalid layout allocation count");
}
