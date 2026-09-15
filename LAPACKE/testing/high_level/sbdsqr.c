#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sbdsqr. */
#define LAPACKE_SBDSQR_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_pos(LD * LD, d);                                    \
        lapacke_test_sfill_vec(LD * LD, e);                                    \
        lapacke_test_sfill_rhs(layout, N, N, vt, LD);                          \
        lapacke_test_sfill_rhs(layout, M, N, u, LD);                           \
        lapacke_test_sfill_rhs(layout, N, NRHS, c, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sbdsqr)(layout, 'U', N, N, M,    \
                                                      NRHS, d, e, vt, LD, u,   \
                                                      LD, c, LD),              \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sbdsqr)
{
    float d[LD * LD];
    float e[LD * LD];
    float vt[LD * LD];
    float u[LD * LD];
    float c[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sbdsqr c", l, N, NRHS, c, LD, lapacke_test_region_full, -13,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_rhs(layout, N, N, vt, LD),
             lapacke_test_sfill_rhs(layout, M, N, u, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_sbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sbdsqr d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -7,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_rhs(layout, N, N, vt, LD),
             lapacke_test_sfill_rhs(layout, M, N, u, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_sbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sbdsqr e", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -8,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_rhs(layout, N, N, vt, LD),
             lapacke_test_sfill_rhs(layout, M, N, u, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_sbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sbdsqr u", l, M, N, u, LD, lapacke_test_region_full, -11,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_rhs(layout, N, N, vt, LD),
             lapacke_test_sfill_rhs(layout, M, N, u, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_sbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sbdsqr vt", l, N, N, vt, LD, lapacke_test_region_full, -9,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_rhs(layout, N, N, vt, LD),
             lapacke_test_sfill_rhs(layout, M, N, u, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, c, LD)),
            API_SUFFIX(LAPACKE_sbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, N, NRHS, c, LD);
        lapacke_test_sfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_sfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_sfill_nan(layout, M, N, u, LD);
        lapacke_test_sfill_nan(layout, N, N, vt, LD);
        lapacke_test_check(
            "sbdsqr NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sbdsqr)(layout, 'U', N, N, M, NRHS, d, e, vt, LD,
                                       u, LD, c, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SBDSQR_ALLOC_TEST(0, 0, "sbdsqr work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SBDSQR_ALLOC_TEST(0, 1, "sbdsqr allocation count", 0);
    lapacke_test_check_alloc_count("sbdsqr col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SBDSQR_ALLOC_TEST(1, 0, "sbdsqr work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SBDSQR_ALLOC_TEST(1, 1, "sbdsqr transpose alloc failure (vt_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SBDSQR_ALLOC_TEST(1, 2, "sbdsqr transpose alloc failure (u_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SBDSQR_ALLOC_TEST(1, 3, "sbdsqr transpose alloc failure (c_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SBDSQR_ALLOC_TEST(1, 4, "sbdsqr allocation count", 0);
    lapacke_test_check_alloc_count("sbdsqr row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SBDSQR_ALLOC_TEST(2, 0, "sbdsqr invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sbdsqr invalid layout allocation count");
}
