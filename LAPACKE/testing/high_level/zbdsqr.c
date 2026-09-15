#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zbdsqr. */
#define LAPACKE_ZBDSQR_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_pos(LD * LD, d);                                    \
        lapacke_test_dfill_vec(LD * LD, e);                                    \
        lapacke_test_zfill_rhs(layout, N, N, vt, LD);                          \
        lapacke_test_zfill_rhs(layout, M, N, u, LD);                           \
        lapacke_test_zfill_rhs(layout, N, NRHS, c, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zbdsqr)(layout, 'U', N, N, M,    \
                                                      NRHS, d, e, vt, LD, u,   \
                                                      LD, c, LD),              \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zbdsqr)
{
    double d[LD * LD];
    double e[LD * LD];
    lapack_complex_double vt[LD * LD];
    lapack_complex_double u[LD * LD];
    lapack_complex_double c[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zbdsqr c", l, N, NRHS, c, LD, lapacke_test_region_full, -13,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, e),
             lapacke_test_zfill_rhs(layout, N, N, vt, LD),
             lapacke_test_zfill_rhs(layout, M, N, u, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_zbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "zbdsqr d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -7,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, e),
             lapacke_test_zfill_rhs(layout, N, N, vt, LD),
             lapacke_test_zfill_rhs(layout, M, N, u, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_zbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "zbdsqr e", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -8,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, e),
             lapacke_test_zfill_rhs(layout, N, N, vt, LD),
             lapacke_test_zfill_rhs(layout, M, N, u, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_zbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zbdsqr u", l, M, N, u, LD, lapacke_test_region_full, -11,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, e),
             lapacke_test_zfill_rhs(layout, N, N, vt, LD),
             lapacke_test_zfill_rhs(layout, M, N, u, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_zbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zbdsqr vt", l, N, N, vt, LD, lapacke_test_region_full, -9,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, e),
             lapacke_test_zfill_rhs(layout, N, N, vt, LD),
             lapacke_test_zfill_rhs(layout, M, N, u, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_zbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, N, NRHS, c, LD);
        lapacke_test_dfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_dfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_zfill_nan(layout, M, N, u, LD);
        lapacke_test_zfill_nan(layout, N, N, vt, LD);
        lapacke_test_check(
            "zbdsqr NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZBDSQR_ALLOC_TEST(0, 0, "zbdsqr work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZBDSQR_ALLOC_TEST(0, 1, "zbdsqr allocation count", 0);
    lapacke_test_check_alloc_count("zbdsqr col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZBDSQR_ALLOC_TEST(1, 0, "zbdsqr work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZBDSQR_ALLOC_TEST(1, 1, "zbdsqr transpose alloc failure (vt_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZBDSQR_ALLOC_TEST(1, 2, "zbdsqr transpose alloc failure (u_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZBDSQR_ALLOC_TEST(1, 3, "zbdsqr transpose alloc failure (c_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZBDSQR_ALLOC_TEST(1, 4, "zbdsqr allocation count", 0);
    lapacke_test_check_alloc_count("zbdsqr row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZBDSQR_ALLOC_TEST(2, 0, "zbdsqr invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zbdsqr invalid layout allocation count");
}
