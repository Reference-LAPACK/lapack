#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sgelsd. */
#define LAPACKE_SGELSD_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, M, N, a, LD);                               \
        lapacke_test_sfill_rhs(layout, M, NRHS, b, LD);                        \
        lapacke_test_sfill_pos(LD * LD, s);                                    \
        rcond[0] = 0.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sgelsd)(layout, M, N, NRHS, a,   \
                                                      LD, b, LD, s, rcond[0],  \
                                                      &rank),                  \
                           expected);                                          \
    } while (0)

/* B has max(M, N) rows, all of them documented inputs for M > N. For
 * M < N only the first M rows are, while LAPACKE checks all max(M, N):
 * a divergence this fixture does not reach. */
LAPACKE_TEST(sgelsd)
{
    float a[LD * LD];
    float b[LD * LD];
    float s[LD * LD];
    float rcond[1];
    lapack_int rank;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgelsd a", l, M, N, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_sfill_rhs(layout, M, NRHS, b, LD),
             lapacke_test_sfill_pos(LD * LD, s), (rcond[0] = 0.0f)),
            API_SUFFIX(LAPACKE_sgelsd)(layout, M, N, NRHS, a, LD, b, LD, s,
                                       rcond[0], &rank));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgelsd b", l, M, NRHS, b, LD, lapacke_test_region_full, -7,
            (lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_sfill_rhs(layout, M, NRHS, b, LD),
             lapacke_test_sfill_pos(LD * LD, s), (rcond[0] = 0.0f)),
            API_SUFFIX(LAPACKE_sgelsd)(layout, M, N, NRHS, a, LD, b, LD, s,
                                       rcond[0], &rank));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgelsd rcond", l, 1, 1, rcond, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -10,
            (lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_sfill_rhs(layout, M, NRHS, b, LD),
             lapacke_test_sfill_pos(LD * LD, s), (rcond[0] = 0.0f)),
            API_SUFFIX(LAPACKE_sgelsd)(layout, M, N, NRHS, a, LD, b, LD, s,
                                       rcond[0], &rank));

        /* NaN checks off: no test, xGELSD scales the bidiagonal problem by its
         * norm, and xLASCL rejects a NaN scale factor through XERBLA, which
         * stops the process. */
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SGELSD_ALLOC_TEST(0, 0, "sgelsd work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGELSD_ALLOC_TEST(0, 1, "sgelsd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGELSD_ALLOC_TEST(0, 2, "sgelsd allocation count", 0);
    lapacke_test_check_alloc_count("sgelsd col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SGELSD_ALLOC_TEST(1, 0, "sgelsd work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGELSD_ALLOC_TEST(1, 1, "sgelsd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGELSD_ALLOC_TEST(1, 2, "sgelsd transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGELSD_ALLOC_TEST(1, 3, "sgelsd transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGELSD_ALLOC_TEST(1, 4, "sgelsd allocation count", 0);
    lapacke_test_check_alloc_count("sgelsd row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SGELSD_ALLOC_TEST(2, 0, "sgelsd invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgelsd invalid layout allocation count");
}
