#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cgesvx. */
#define LAPACKE_CGESVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, N, N, a, LD);                               \
        lapacke_test_cfill(layout, N, N, af, LD);                              \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        equed = 'N';                                                           \
        lapacke_test_sfill_pos(LD * LD, r);                                    \
        lapacke_test_sfill_pos(LD * LD, c);                                    \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_cfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_cgesvx)(layout, 'N', 'N', N, NRHS, a, LD, af,   \
                                       LD, ipiv, &equed, r, c, b, LD, x, LD,   \
                                       &rcond, ferr, berr, &rpivot),           \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(cgesvx)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float af[LD * LD];
    lapack_int ipiv[LD * LD];
    char equed;
    float r[LD * LD];
    float c[LD * LD];
    lapack_complex_float b[LD * LD];
    lapack_complex_float x[LD * LD];
    float rcond;
    float ferr[LD * LD];
    float berr[LD * LD];
    float rpivot;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cgesvx a fact=N", l, N, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_cfill(layout, N, N, a, LD),
             lapacke_test_cfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgesvx)(layout, 'N', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgesvx b fact=N", l, N, NRHS, b, LD, lapacke_test_region_full, -14,
            (lapacke_test_cfill(layout, N, N, a, LD),
             lapacke_test_cfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgesvx)(layout, 'N', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgesvx a fact=F", l, N, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_cfill(layout, N, N, a, LD),
             lapacke_test_cfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_cgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_cgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgesvx af fact=F", l, N, N, af, LD, lapacke_test_region_full, -8,
            (lapacke_test_cfill(layout, N, N, a, LD),
             lapacke_test_cfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_cgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_cgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgesvx b fact=F", l, N, NRHS, b, LD, lapacke_test_region_full, -14,
            (lapacke_test_cfill(layout, N, N, a, LD),
             lapacke_test_cfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_cgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_cgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_SNAN_SWEEP(
            "cgesvx c fact=F", l, N, 1, c, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -13,
            (lapacke_test_cfill(layout, N, N, a, LD),
             lapacke_test_cfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_cgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_cgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_SNAN_SWEEP(
            "cgesvx r fact=F", l, N, 1, r, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -12,
            (lapacke_test_cfill(layout, N, N, a, LD),
             lapacke_test_cfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_cgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_cgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        equed = 'N';
        lapacke_test_cfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_cfill_nan(layout, N, N, af, LD);
        lapacke_test_sfill_nan(layout, N, 1, c, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_sfill_nan(layout, N, 1, r, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_check(
            "cgesvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgesvx)(layout, 'N', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CGESVX_ALLOC_TEST(0, 0, "cgesvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGESVX_ALLOC_TEST(0, 1, "cgesvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGESVX_ALLOC_TEST(0, 2, "cgesvx allocation count", 0);
    lapacke_test_check_alloc_count("cgesvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CGESVX_ALLOC_TEST(1, 0, "cgesvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGESVX_ALLOC_TEST(1, 1, "cgesvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGESVX_ALLOC_TEST(1, 2, "cgesvx transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGESVX_ALLOC_TEST(1, 3, "cgesvx transpose alloc failure (af_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGESVX_ALLOC_TEST(1, 4, "cgesvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGESVX_ALLOC_TEST(1, 5, "cgesvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGESVX_ALLOC_TEST(1, 6, "cgesvx allocation count", 0);
    lapacke_test_check_alloc_count("cgesvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CGESVX_ALLOC_TEST(2, 0, "cgesvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cgesvx invalid layout allocation count");
}
