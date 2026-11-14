#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One call: refill the positive definite tridiagonal (d, e), schedule the
 * malloc failure countdown and check the info result of a cpteqr call with
 * the given compz in the indexed layout. */
#define LAPACKE_CPTEQR_TEST(layout_index, compz, countdown, name, expected)    \
    do {                                                                       \
        for (lapack_int i = 0; i < N; i++) {                                   \
            d[i] = 4.0 + i;                                                    \
            e[i] = 1.0;                                                        \
        }                                                                      \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cpteqr)(                         \
                               lapacke_test_layouts[layout_index], compz, N,   \
                               d, e, z, LD),                                   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cpteqr)
{
    float d[N], e[N];
    lapack_complex_float z[LD * LD];

    /* CPTEQR needs a real WORK of 4*N entries for every compz: with
     * compz = 'N' it computes the eigenvalues with SLASQ1 (through CBDSQR),
     * which uses all of them. The test allocator fails the test when a call
     * writes past the end of a LAPACKE allocation. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_CPTEQR_TEST(l, 'N', -1, "cpteqr compz = N", 0);
        LAPACKE_CPTEQR_TEST(l, 'I', -1, "cpteqr compz = I", 0);
    }

    /* Column-major: the workspace is the only allocation. */
    LAPACKE_CPTEQR_TEST(0, 'I', 0, "cpteqr work alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_CPTEQR_TEST(0, 'I', 1, "cpteqr allocation count", 0);
    lapacke_test_check_alloc_count("cpteqr col-major allocation count");

    /* Row-major allocates the workspace, then the transposed copy of Z. */
    LAPACKE_CPTEQR_TEST(1, 'I', 0, "cpteqr work alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPTEQR_TEST(1, 'I', 1, "cpteqr transpose alloc failure (z)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_CPTEQR_TEST(1, 'I', 2, "cpteqr allocation count", 0);
    lapacke_test_check_alloc_count("cpteqr row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_CPTEQR_TEST(2, 'I', 0, "cpteqr invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cpteqr invalid layout allocation count");
}
