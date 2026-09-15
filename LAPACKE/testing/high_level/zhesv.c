#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zhesv. */
#define LAPACKE_ZHESV_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_sym(layout, N, a, LD);                              \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zhesv)(layout, 'U', N, NRHS, a,  \
                                                     LD, ipiv, b, LD),         \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zhesv)
{
    lapack_complex_double a[LD * LD];
    lapack_int ipiv[LD * LD];
    lapack_complex_double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zhesv a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -5,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zhesv)(layout, 'U', N, NRHS, a, LD, ipiv, b,
                                      LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zhesv b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zhesv)(layout, 'U', N, NRHS, a, LD, ipiv, b,
                                      LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zhesv a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -5,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zhesv)(layout, 'L', N, NRHS, a, LD, ipiv, b,
                                      LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zhesv b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zhesv)(layout, 'L', N, NRHS, a, LD, ipiv, b,
                                      LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_zfill_nan(layout, N, N, a, LD);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check("zhesv NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_zhesv)(layout, 'U', N, NRHS, a,
                                                     LD, ipiv, b, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZHESV_ALLOC_TEST(0, 0, "zhesv work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZHESV_ALLOC_TEST(0, 1, "zhesv allocation count", 0);
    lapacke_test_check_alloc_count("zhesv col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZHESV_ALLOC_TEST(1, 0, "zhesv work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZHESV_ALLOC_TEST(1, 1, "zhesv transpose alloc failure (a_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZHESV_ALLOC_TEST(1, 2, "zhesv transpose alloc failure (b_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZHESV_ALLOC_TEST(1, 3, "zhesv allocation count", 0);
    lapacke_test_check_alloc_count("zhesv row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZHESV_ALLOC_TEST(2, 0, "zhesv invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zhesv invalid layout allocation count");
}
