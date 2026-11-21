#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One call computing only the singular values: refill the rows-by-cols A,
 * schedule the malloc failure countdown and check the info result of a
 * cgejsv call in the indexed layout. */
#define LAPACKE_CGEJSV_TEST(layout_index, rows, cols, countdown, name,         \
                            expected)                                          \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, rows, cols, a, LD);                         \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cgejsv)(                         \
                               layout, 'C', 'N', 'N', 'N', 'N', 'N', rows,     \
                               cols, a, LD, sva, u, LD, v, LD, stat, istat),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cgejsv)
{
    lapack_complex_float a[LD * LD], u[LD * LD], v[LD * LD];
    float sva[LD], stat[7];
    lapack_int istat[3];

    /* CGEJSV always sets IWORK(1:4), so LAPACKE must allocate at least four
     * entries of IWORK even when M+2*N is smaller, as for a 1-by-1 matrix.
     * The test allocator fails the test when a call writes past the end of
     * a LAPACKE allocation. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_CGEJSV_TEST(l, 1, 1, -1, "cgejsv 1-by-1", 0);
        LAPACKE_CGEJSV_TEST(l, M, N, -1, "cgejsv M-by-N", 0);
    }

    /* Column-major: IWORK, the complex and the real workspace are the only
     * allocations. */
    LAPACKE_CGEJSV_TEST(0, M, N, 0, "cgejsv iwork alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGEJSV_TEST(0, M, N, 1, "cgejsv cwork alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGEJSV_TEST(0, M, N, 2, "cgejsv rwork alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_CGEJSV_TEST(0, M, N, 3, "cgejsv allocation count", 0);
    lapacke_test_check_alloc_count("cgejsv col-major allocation count");

    /* Row-major allocates the same workspace, then the transposed copy of A
     * (U and V are not referenced). */
    LAPACKE_CGEJSV_TEST(1, M, N, 3, "cgejsv transpose alloc failure (a)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_CGEJSV_TEST(1, M, N, 4, "cgejsv allocation count", 0);
    lapacke_test_check_alloc_count("cgejsv row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_CGEJSV_TEST(2, M, N, 0, "cgejsv invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cgejsv invalid layout allocation count");
}
