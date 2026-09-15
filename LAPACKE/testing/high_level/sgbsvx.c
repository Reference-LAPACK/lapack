#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

#define REGION_AB(i, j) lapacke_test_region_band(i, j, N, KU)
#define REGION_AFB(i, j) lapacke_test_region_band(i, j, N, KL + KU)

/* Refill the inputs, schedule the malloc failure, call sgbsvx. */
#define LAPACKE_SGBSVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD);      \
        lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb, LD); \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        equed = 'N';                                                           \
        lapacke_test_sfill_pos(LD * LD, r);                                    \
        lapacke_test_sfill_pos(LD * LD, c);                                    \
        lapacke_test_sfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_sgbsvx)(layout, 'N', 'N', N, KL, KU, NRHS, ab,  \
                                       LD, afb, LD, ipiv, &equed, r, c, b, LD, \
                                       x, LD, &rcond, ferr, berr, &rpivot),    \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(sgbsvx)
{
    float ab[LD * LD];
    float afb[LD * LD];
    lapack_int ipiv[LD * LD];
    char equed;
    float r[LD * LD];
    float c[LD * LD];
    float b[LD * LD];
    float x[LD * LD];
    float rcond;
    float ferr[LD * LD];
    float berr[LD * LD];
    float rpivot;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgbsvx ab fact=N", l, KL + KU + 1, N, ab, LD, REGION_AB, -8,
            (lapacke_test_sfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sgbsvx)(layout, 'N', 'N', N, KL, KU, NRHS, ab,
                                       LD, afb, LD, ipiv, &equed, r, c, b, LD,
                                       x, LD, &rcond, ferr, berr, &rpivot));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgbsvx b fact=N", l, N, NRHS, b, LD, lapacke_test_region_full, -16,
            (lapacke_test_sfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sgbsvx)(layout, 'N', 'N', N, KL, KU, NRHS, ab,
                                       LD, afb, LD, ipiv, &equed, r, c, b, LD,
                                       x, LD, &rcond, ferr, berr, &rpivot));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgbsvx ab fact=F", l, KL + KU + 1, N, ab, LD, REGION_AB, -8,
            (lapacke_test_sfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_sgbtrf)(layout, N, N, KL, KU, afb, LD, ipiv)),
            API_SUFFIX(LAPACKE_sgbsvx)(layout, 'F', 'N', N, KL, KU, NRHS, ab,
                                       LD, afb, LD, ipiv, &equed, r, c, b, LD,
                                       x, LD, &rcond, ferr, berr, &rpivot));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgbsvx afb fact=F", l, 2 * KL + KU + 1, N, afb, LD, REGION_AFB,
            -10,
            (lapacke_test_sfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_sgbtrf)(layout, N, N, KL, KU, afb, LD, ipiv)),
            API_SUFFIX(LAPACKE_sgbsvx)(layout, 'F', 'N', N, KL, KU, NRHS, ab,
                                       LD, afb, LD, ipiv, &equed, r, c, b, LD,
                                       x, LD, &rcond, ferr, berr, &rpivot));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgbsvx b fact=F", l, N, NRHS, b, LD, lapacke_test_region_full, -16,
            (lapacke_test_sfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_sgbtrf)(layout, N, N, KL, KU, afb, LD, ipiv)),
            API_SUFFIX(LAPACKE_sgbsvx)(layout, 'F', 'N', N, KL, KU, NRHS, ab,
                                       LD, afb, LD, ipiv, &equed, r, c, b, LD,
                                       x, LD, &rcond, ferr, berr, &rpivot));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgbsvx c fact=F", l, N, 1, c, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -15,
            (lapacke_test_sfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_sgbtrf)(layout, N, N, KL, KU, afb, LD, ipiv)),
            API_SUFFIX(LAPACKE_sgbsvx)(layout, 'F', 'N', N, KL, KU, NRHS, ab,
                                       LD, afb, LD, ipiv, &equed, r, c, b, LD,
                                       x, LD, &rcond, ferr, berr, &rpivot));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgbsvx r fact=F", l, N, 1, r, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -14,
            (lapacke_test_sfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_sfill_pos(LD * LD, r),
             lapacke_test_sfill_pos(LD * LD, c),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_sgbtrf)(layout, N, N, KL, KU, afb, LD, ipiv)),
            API_SUFFIX(LAPACKE_sgbsvx)(layout, 'F', 'N', N, KL, KU, NRHS, ab,
                                       LD, afb, LD, ipiv, &equed, r, c, b, LD,
                                       x, LD, &rcond, ferr, berr, &rpivot));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        equed = 'N';
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_sfill_nan(layout, KL + KU + 1, N, ab, LD);
        lapacke_test_sfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_sfill_nan(layout, 2 * KL + KU + 1, N, afb, LD);
        lapacke_test_sfill_nan(layout, N, 1, c, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_sfill_nan(layout, N, 1, r, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_check(
            "sgbsvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sgbsvx)(layout, 'N', 'N', N, KL, KU, NRHS, ab,
                                       LD, afb, LD, ipiv, &equed, r, c, b, LD,
                                       x, LD, &rcond, ferr, berr, &rpivot) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SGBSVX_ALLOC_TEST(0, 0, "sgbsvx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGBSVX_ALLOC_TEST(0, 1, "sgbsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGBSVX_ALLOC_TEST(0, 2, "sgbsvx allocation count", 0);
    lapacke_test_check_alloc_count("sgbsvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SGBSVX_ALLOC_TEST(1, 0, "sgbsvx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGBSVX_ALLOC_TEST(1, 1, "sgbsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGBSVX_ALLOC_TEST(1, 2, "sgbsvx transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGBSVX_ALLOC_TEST(1, 3, "sgbsvx transpose alloc failure (afb_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGBSVX_ALLOC_TEST(1, 4, "sgbsvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGBSVX_ALLOC_TEST(1, 5, "sgbsvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGBSVX_ALLOC_TEST(1, 6, "sgbsvx allocation count", 0);
    lapacke_test_check_alloc_count("sgbsvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SGBSVX_ALLOC_TEST(2, 0, "sgbsvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgbsvx invalid layout allocation count");
}
