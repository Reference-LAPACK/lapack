#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dgesvx. */
#define LAPACKE_DGESVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, N, N, a, LD);                               \
        lapacke_test_dfill(layout, N, N, af, LD);                              \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        equed = 'N';                                                           \
        lapacke_test_dfill_pos(LD * LD, r);                                    \
        lapacke_test_dfill_pos(LD * LD, c);                                    \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_dfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dgesvx)(layout, 'N', 'N', N, NRHS, a, LD, af,   \
                                       LD, ipiv, &equed, r, c, b, LD, x, LD,   \
                                       &rcond, ferr, berr, &rpivot),           \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(dgesvx)
{
    double a[LD * LD];
    double af[LD * LD];
    lapack_int ipiv[LD * LD];
    char equed;
    double r[LD * LD];
    double c[LD * LD];
    double b[LD * LD];
    double x[LD * LD];
    double rcond;
    double ferr[LD * LD];
    double berr[LD * LD];
    double rpivot;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgesvx a fact=N", l, N, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_dfill(layout, N, N, a, LD),
             lapacke_test_dfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgesvx)(layout, 'N', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgesvx b fact=N", l, N, NRHS, b, LD, lapacke_test_region_full, -14,
            (lapacke_test_dfill(layout, N, N, a, LD),
             lapacke_test_dfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgesvx)(layout, 'N', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgesvx a fact=F", l, N, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_dfill(layout, N, N, a, LD),
             lapacke_test_dfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_dgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_dgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgesvx af fact=F", l, N, N, af, LD, lapacke_test_region_full, -8,
            (lapacke_test_dfill(layout, N, N, a, LD),
             lapacke_test_dfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_dgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_dgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgesvx b fact=F", l, N, NRHS, b, LD, lapacke_test_region_full, -14,
            (lapacke_test_dfill(layout, N, N, a, LD),
             lapacke_test_dfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_dgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_dgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgesvx c fact=F", l, N, 1, c, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -13,
            (lapacke_test_dfill(layout, N, N, a, LD),
             lapacke_test_dfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_dgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_dgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgesvx r fact=F", l, N, 1, r, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -12,
            (lapacke_test_dfill(layout, N, N, a, LD),
             lapacke_test_dfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (equed = 'B'),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_dgetrf)(layout, N, N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_dgesvx)(layout, 'F', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        equed = 'N';
        lapacke_test_dfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_dfill_nan(layout, N, N, a, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, N, af, LD);
        lapacke_test_dfill_nan(layout, N, 1, c, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_dfill_nan(layout, N, 1, r, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_check(
            "dgesvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgesvx)(layout, 'N', 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, &equed, r, c, b, LD, x, LD, &rcond,
                                       ferr, berr, &rpivot) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DGESVX_ALLOC_TEST(0, 0, "dgesvx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGESVX_ALLOC_TEST(0, 1, "dgesvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGESVX_ALLOC_TEST(0, 2, "dgesvx allocation count", 0);
    lapacke_test_check_alloc_count("dgesvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGESVX_ALLOC_TEST(1, 0, "dgesvx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGESVX_ALLOC_TEST(1, 1, "dgesvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGESVX_ALLOC_TEST(1, 2, "dgesvx transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGESVX_ALLOC_TEST(1, 3, "dgesvx transpose alloc failure (af_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGESVX_ALLOC_TEST(1, 4, "dgesvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGESVX_ALLOC_TEST(1, 5, "dgesvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGESVX_ALLOC_TEST(1, 6, "dgesvx allocation count", 0);
    lapacke_test_check_alloc_count("dgesvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGESVX_ALLOC_TEST(2, 0, "dgesvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgesvx invalid layout allocation count");
}
