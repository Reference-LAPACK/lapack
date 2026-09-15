#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zgtsv. */
#define LAPACKE_ZGTSV_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_vec(LD * LD, dl);                                   \
        lapacke_test_zfill_pos(LD * LD, d);                                    \
        lapacke_test_zfill_vec(LD * LD, du);                                   \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zgtsv)(layout, N, NRHS, dl, d, du, b, LD),      \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(zgtsv)
{
    lapack_complex_double dl[LD * LD];
    lapack_complex_double d[LD * LD];
    lapack_complex_double du[LD * LD];
    lapack_complex_double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgtsv b", l, N, NRHS, b, LD, lapacke_test_region_full, -7,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zgtsv)(layout, N, NRHS, dl, d, du, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgtsv d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -5,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zgtsv)(layout, N, NRHS, dl, d, du, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgtsv dl", l, N - 1, 1, dl, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -4,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zgtsv)(layout, N, NRHS, dl, d, du, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgtsv du", l, N - 1, 1, du, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -6,
            (lapacke_test_zfill_vec(LD * LD, dl),
             lapacke_test_zfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, du),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zgtsv)(layout, N, NRHS, dl, d, du, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_zfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_zfill_nan(layout, N - 1, 1, dl,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_zfill_nan(layout, N - 1, 1, du,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check(
            "zgtsv NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgtsv)(layout, N, NRHS, dl, d, du, b, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZGTSV_ALLOC_TEST(0, 0, "zgtsv allocation count", 0);
    lapacke_test_check_alloc_count("zgtsv col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZGTSV_ALLOC_TEST(1, 0, "zgtsv transpose alloc failure (b_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGTSV_ALLOC_TEST(1, 1, "zgtsv allocation count", 0);
    lapacke_test_check_alloc_count("zgtsv row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZGTSV_ALLOC_TEST(2, 0, "zgtsv invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgtsv invalid layout allocation count");
}
