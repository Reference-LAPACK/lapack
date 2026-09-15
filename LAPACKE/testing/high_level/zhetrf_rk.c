#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zhetrf_rk. */
#define LAPACKE_ZHETRF_RK_ALLOC_TEST(layout_index, countdown, name, expected)  \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_sym(layout, N, a, LD);                              \
        lapacke_test_zfill_vec(LD * LD, e);                                    \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zhetrf_rk)(layout, 'U', N, a, LD, e, ipiv),     \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(zhetrf_rk)
{
    lapack_complex_double a[LD * LD];
    lapack_complex_double e[LD * LD];
    lapack_int ipiv[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zhetrf_rk a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -4,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_zhetrf_rk)(layout, 'U', N, a, LD, e, ipiv));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zhetrf_rk a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -4,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_zhetrf_rk)(layout, 'L', N, a, LD, e, ipiv));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_vec(LD * LD, e);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_zfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "zhetrf_rk NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zhetrf_rk)(layout, 'U', N, a, LD, e, ipiv) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZHETRF_RK_ALLOC_TEST(0, 0, "zhetrf_rk work alloc failure (work)",
                                 LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZHETRF_RK_ALLOC_TEST(0, 1, "zhetrf_rk allocation count", 0);
    lapacke_test_check_alloc_count("zhetrf_rk col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZHETRF_RK_ALLOC_TEST(1, 0, "zhetrf_rk work alloc failure (work)",
                                 LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZHETRF_RK_ALLOC_TEST(1, 1,
                                 "zhetrf_rk transpose alloc failure (a_t)",
                                 LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZHETRF_RK_ALLOC_TEST(1, 2, "zhetrf_rk allocation count", 0);
    lapacke_test_check_alloc_count("zhetrf_rk row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZHETRF_RK_ALLOC_TEST(2, 0, "zhetrf_rk invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zhetrf_rk invalid layout allocation count");
}
