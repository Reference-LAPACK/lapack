#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cgesvd. */
#define LAPACKE_CGESVD_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, M, N, a, LD);                               \
        lapacke_test_sfill_pos(LD * LD, s);                                    \
        lapacke_test_cfill_rhs(layout, M, M, u, LD);                           \
        lapacke_test_cfill_rhs(layout, N, N, vt, LD);                          \
        lapacke_test_sfill_vec(LD * LD, superb);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cgesvd)(layout, 'A', 'A', M, N,  \
                                                      a, LD, s, u, LD, vt, LD, \
                                                      superb),                 \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cgesvd)
{
    lapack_complex_float a[LD * LD];
    float s[LD * LD];
    lapack_complex_float u[LD * LD];
    lapack_complex_float vt[LD * LD];
    float superb[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cgesvd a", l, M, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_cfill(layout, M, N, a, LD),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_cfill_rhs(layout, M, M, u, LD),
             lapacke_test_cfill_rhs(layout, N, N, vt, LD),
             lapacke_test_sfill_vec(LD * LD, superb)),
            API_SUFFIX(LAPACKE_cgesvd)(layout, 'A', 'A', M, N, a, LD, s, u, LD,
                                       vt, LD, superb));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_pos(LD * LD, s);
        lapacke_test_cfill_rhs(layout, M, M, u, LD);
        lapacke_test_cfill_rhs(layout, N, N, vt, LD);
        lapacke_test_sfill_vec(LD * LD, superb);
        lapacke_test_cfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "cgesvd NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgesvd)(layout, 'A', 'A', M, N, a, LD, s, u, LD,
                                       vt, LD, superb) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CGESVD_ALLOC_TEST(0, 0, "cgesvd work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGESVD_ALLOC_TEST(0, 1, "cgesvd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGESVD_ALLOC_TEST(0, 2, "cgesvd allocation count", 0);
    lapacke_test_check_alloc_count("cgesvd col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CGESVD_ALLOC_TEST(1, 0, "cgesvd work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGESVD_ALLOC_TEST(1, 1, "cgesvd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGESVD_ALLOC_TEST(1, 2, "cgesvd transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGESVD_ALLOC_TEST(1, 3, "cgesvd transpose alloc failure (u_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGESVD_ALLOC_TEST(1, 4, "cgesvd transpose alloc failure (vt_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGESVD_ALLOC_TEST(1, 5, "cgesvd allocation count", 0);
    lapacke_test_check_alloc_count("cgesvd row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CGESVD_ALLOC_TEST(2, 0, "cgesvd invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cgesvd invalid layout allocation count");
}
