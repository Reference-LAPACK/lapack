#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

#define REGION_E_U(i, j) ((i) >= 1)
#define REGION_E_L(i, j) ((i) < N - 1)

/* Refill the inputs, schedule the malloc failure, call csycon_3. */
#define LAPACKE_CSYCON_3_ALLOC_TEST(layout_index, countdown, name, expected)   \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_sym(layout, N, a, LD);                              \
        lapacke_test_cfill_vec(LD * LD, e);                                    \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        anorm[0] = 1.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_csycon_3)(layout, 'U', N, a, LD, \
                                                        e, ipiv, anorm[0],     \
                                                        &rcond),               \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(csycon_3)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float e[LD * LD];
    lapack_int ipiv[LD * LD];
    float anorm[1];
    float rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "csycon_3 a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -4,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, e),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_csycon_3)(layout, 'U', N, a, LD, e, ipiv,
                                         anorm[0], &rcond));

        LAPACKE_TEST_CNAN_SWEEP(
            "csycon_3 e uplo=U", l, N, 1, e, LAPACKE_TEST_VLD(layout, N),
            REGION_E_U, -6,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, e),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_csycon_3)(layout, 'U', N, a, LD, e, ipiv,
                                         anorm[0], &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "csycon_3 anorm uplo=U", l, 1, 1, anorm,
            LAPACKE_TEST_VLD(layout, 1), lapacke_test_region_full, -8,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, e),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_csycon_3)(layout, 'U', N, a, LD, e, ipiv,
                                         anorm[0], &rcond));

        LAPACKE_TEST_CNAN_SWEEP(
            "csycon_3 a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -4,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, e),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_csycon_3)(layout, 'L', N, a, LD, e, ipiv,
                                         anorm[0], &rcond));

        LAPACKE_TEST_CNAN_SWEEP(
            "csycon_3 e uplo=L", l, N, 1, e, LAPACKE_TEST_VLD(layout, N),
            REGION_E_L, -6,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, e),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_csycon_3)(layout, 'L', N, a, LD, e, ipiv,
                                         anorm[0], &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "csycon_3 anorm uplo=L", l, 1, 1, anorm,
            LAPACKE_TEST_VLD(layout, 1), lapacke_test_region_full, -8,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, e),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_csycon_3)(layout, 'L', N, a, LD, e, ipiv,
                                         anorm[0], &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        anorm[0] = 1.0f;
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        lapacke_test_cfill_nan(layout, N, 1, e, LAPACKE_TEST_VLD(layout, N));
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_check(
            "csycon_3 NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_csycon_3)(layout, 'U', N, a, LD, e, ipiv,
                                         anorm[0], &rcond) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CSYCON_3_ALLOC_TEST(0, 0, "csycon_3 work alloc failure (work)",
                                LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CSYCON_3_ALLOC_TEST(0, 1, "csycon_3 allocation count", 0);
    lapacke_test_check_alloc_count("csycon_3 col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CSYCON_3_ALLOC_TEST(1, 0, "csycon_3 work alloc failure (work)",
                                LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CSYCON_3_ALLOC_TEST(1, 1, "csycon_3 transpose alloc failure (a_t)",
                                LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CSYCON_3_ALLOC_TEST(1, 2, "csycon_3 allocation count", 0);
    lapacke_test_check_alloc_count("csycon_3 row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CSYCON_3_ALLOC_TEST(2, 0, "csycon_3 invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("csycon_3 invalid layout allocation count");
}
