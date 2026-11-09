#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One call: refill the positive definite tridiagonal (d, e), schedule the
 * malloc failure countdown and check the info result of a spteqr call with
 * the given compz in the indexed layout. */
#define LAPACKE_SPTEQR_TEST(layout_index, compz, countdown, name, expected)    \
    do {                                                                       \
        for (lapack_int i = 0; i < N; i++) {                                   \
            d[i] = 4.0 + i;                                                    \
            e[i] = 1.0;                                                        \
        }                                                                      \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_spteqr)(                         \
                               lapacke_test_layouts[layout_index], compz, N,   \
                               d, e, z, LD),                                   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(spteqr)
{
    float d[N], e[N], z[LD * LD];

    /* SPTEQR needs a WORK of 4*N entries for every compz: with compz = 'N'
     * it computes the eigenvalues with SLASQ1 (through SBDSQR), which uses
     * all of them. The test allocator fails the test when a call writes past
     * the end of a LAPACKE allocation. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_SPTEQR_TEST(l, 'N', -1, "spteqr compz = N", 0);
        LAPACKE_SPTEQR_TEST(l, 'I', -1, "spteqr compz = I", 0);
    }

    /* Column-major: the workspace is the only allocation. */
    LAPACKE_SPTEQR_TEST(0, 'I', 0, "spteqr work alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_SPTEQR_TEST(0, 'I', 1, "spteqr allocation count", 0);
    lapacke_test_check_alloc_count("spteqr col-major allocation count");

    /* Row-major allocates the workspace, then the transposed copy of Z. */
    LAPACKE_SPTEQR_TEST(1, 'I', 0, "spteqr work alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPTEQR_TEST(1, 'I', 1, "spteqr transpose alloc failure (z)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_SPTEQR_TEST(1, 'I', 2, "spteqr allocation count", 0);
    lapacke_test_check_alloc_count("spteqr row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_SPTEQR_TEST(2, 'I', 0, "spteqr invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("spteqr invalid layout allocation count");
}
