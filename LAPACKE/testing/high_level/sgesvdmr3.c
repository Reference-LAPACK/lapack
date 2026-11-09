#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One call with singular vectors: refill the rows-by-cols A, schedule the
 * malloc failure countdown and check the info result of a sgesvdmr3 call in
 * the indexed layout. */
#define LAPACKE_SGESVDMR3_TEST(layout_index, rows, cols, countdown, name,      \
                               expected)                                       \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, rows, cols, a, LD);                         \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sgesvdmr3)(                      \
                               layout, 'V', 'V', 'A', rows, cols, a, LD, 0.0f, \
                               0.0f, 0, 0, &ns, s, u, LD, vt, LD, superb),     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sgesvdmr3)
{
    float a[LD * LD], s[LD], u[LD * LD], vt[LD * LD];
    lapack_int ns, superb[12 * LD];

    /* With singular vectors, SGESVDMR3 needs an IWORK of 24*MIN(M,N)
     * entries, which LAPACKE allocates; the test allocator fails the test
     * when a call writes past the end of a LAPACKE allocation. M > N and
     * M < N take the QR and the LQ path. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_SGESVDMR3_TEST(l, M, N, -1, "sgesvdmr3 M > N", 0);
        lapacke_test_check("sgesvdmr3 M > N singular values found",
                           lapacke_test_layout_names[l], ns, N);
        LAPACKE_SGESVDMR3_TEST(l, N, M, -1, "sgesvdmr3 M < N", 0);
        lapacke_test_check("sgesvdmr3 M < N singular values found",
                           lapacke_test_layout_names[l], ns, N);
    }

    /* Column-major: the workspace and IWORK are the only allocations. */
    LAPACKE_SGESVDMR3_TEST(0, M, N, 0, "sgesvdmr3 work alloc failure",
                           LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGESVDMR3_TEST(0, M, N, 1, "sgesvdmr3 iwork alloc failure",
                           LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_SGESVDMR3_TEST(0, M, N, 2, "sgesvdmr3 allocation count", 0);
    lapacke_test_check_alloc_count("sgesvdmr3 col-major allocation count");

    /* Row-major: allocation order is the workspace and IWORK, then the
     * transposed copies of A, U and VT inside the middle level. */
    LAPACKE_SGESVDMR3_TEST(1, M, N, 0, "sgesvdmr3 work alloc failure",
                           LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGESVDMR3_TEST(1, M, N, 1, "sgesvdmr3 iwork alloc failure",
                           LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGESVDMR3_TEST(1, M, N, 2, "sgesvdmr3 transpose alloc failure (a)",
                           LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGESVDMR3_TEST(1, M, N, 3, "sgesvdmr3 transpose alloc failure (u)",
                           LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGESVDMR3_TEST(1, M, N, 4, "sgesvdmr3 transpose alloc failure (vt)",
                           LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_SGESVDMR3_TEST(1, M, N, 5, "sgesvdmr3 allocation count", 0);
    lapacke_test_check_alloc_count("sgesvdmr3 row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_SGESVDMR3_TEST(2, M, N, 0, "sgesvdmr3 invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgesvdmr3 invalid layout allocation count");
}
