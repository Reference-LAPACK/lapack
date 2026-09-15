#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sgecxx. */
#define LAPACKE_SGECXX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_fill_int(LD * LD, desel_rows, 0);                         \
        lapacke_test_fill_int(LD * LD, sel_desel_cols, 0);                     \
        abstol[0] = -1.0f;                                                     \
        reltol[0] = -1.0f;                                                     \
        lapacke_test_sfill(layout, M, N, a, LD);                               \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_fill_int(LD * LD, jpiv, 0);                               \
        lapacke_test_sfill_vec(LD * LD, tau);                                  \
        lapacke_test_sfill_rhs(layout, M, N, c, LD);                           \
        lapacke_test_sfill_rhs(layout, M, N, qrc, LD);                         \
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sgecxx)(                         \
                               layout, 'X', 'N', M, N, desel_rows,             \
                               sel_desel_cols, N, abstol[0], reltol[0], a, LD, \
                               &k, &maxc2nrmk, &relmaxc2nrmk, &fnrmk, ipiv,    \
                               jpiv, tau, c, LD, qrc, LD, x, LD),              \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sgecxx)
{
    lapack_int desel_rows[LD * LD];
    lapack_int sel_desel_cols[LD * LD];
    float abstol[1];
    float reltol[1];
    float a[LD * LD];
    lapack_int k;
    float maxc2nrmk;
    float relmaxc2nrmk;
    float fnrmk;
    lapack_int ipiv[LD * LD];
    lapack_int jpiv[LD * LD];
    float tau[LD * LD];
    float c[LD * LD];
    float qrc[LD * LD];
    float x[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgecxx a", l, M, N, a, LD, lapacke_test_region_full, -11,
            (lapacke_test_fill_int(LD * LD, desel_rows, 0),
             lapacke_test_fill_int(LD * LD, sel_desel_cols, 0),
             (abstol[0] = -1.0f), (reltol[0] = -1.0f),
             lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_fill_int(LD * LD, jpiv, 0),
             lapacke_test_sfill_vec(LD * LD, tau),
             lapacke_test_sfill_rhs(layout, M, N, c, LD),
             lapacke_test_sfill_rhs(layout, M, N, qrc, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD)),
            API_SUFFIX(LAPACKE_sgecxx)(
                layout, 'X', 'N', M, N, desel_rows, sel_desel_cols, N,
                abstol[0], reltol[0], a, LD, &k, &maxc2nrmk, &relmaxc2nrmk,
                &fnrmk, ipiv, jpiv, tau, c, LD, qrc, LD, x, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_int(LD * LD, desel_rows, 0);
        lapacke_test_fill_int(LD * LD, sel_desel_cols, 0);
        abstol[0] = -1.0f;
        reltol[0] = -1.0f;
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_fill_int(LD * LD, jpiv, 0);
        lapacke_test_sfill_vec(LD * LD, tau);
        lapacke_test_sfill_rhs(layout, M, N, c, LD);
        lapacke_test_sfill_rhs(layout, M, N, qrc, LD);
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_sfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "sgecxx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sgecxx)(
                layout, 'X', 'N', M, N, desel_rows, sel_desel_cols, N,
                abstol[0], reltol[0], a, LD, &k, &maxc2nrmk, &relmaxc2nrmk,
                &fnrmk, ipiv, jpiv, tau, c, LD, qrc, LD, x, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SGECXX_ALLOC_TEST(0, 0, "sgecxx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGECXX_ALLOC_TEST(0, 1, "sgecxx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGECXX_ALLOC_TEST(0, 2, "sgecxx allocation count", 0);
    lapacke_test_check_alloc_count("sgecxx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SGECXX_ALLOC_TEST(1, 0, "sgecxx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGECXX_ALLOC_TEST(1, 1, "sgecxx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGECXX_ALLOC_TEST(1, 2, "sgecxx transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGECXX_ALLOC_TEST(1, 3, "sgecxx transpose alloc failure (c_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGECXX_ALLOC_TEST(1, 4, "sgecxx transpose alloc failure (qrc_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGECXX_ALLOC_TEST(1, 5, "sgecxx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGECXX_ALLOC_TEST(1, 6, "sgecxx allocation count", 0);
    lapacke_test_check_alloc_count("sgecxx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SGECXX_ALLOC_TEST(2, 0, "sgecxx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgecxx invalid layout allocation count");
}
