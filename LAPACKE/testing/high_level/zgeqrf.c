#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One allocation test: refill the input, schedule the malloc failure
 * countdown and check the info result of a zgeqrf call in the indexed
 * layout. */
#define LAPACKE_ZGEQRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, M, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zgeqrf)(layout, M, N, a, LD, tau), expected);   \
    } while (0)

LAPACKE_TEST(zgeqrf)
{
    lapack_complex_double a[LD * LD], tau[N];

    /* The whole M-by-N matrix is a documented input (tau is output only). */
    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];
        LAPACKE_TEST_ZNAN_SWEEP(
            "zgeqrf a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_zfill(layout, M, N, a, LD)),
            API_SUFFIX(LAPACKE_zgeqrf)(layout, M, N, a, LD, tau));

        /* With NaN checking disabled even all-NaN input must go through to
         * the Fortran routine (valid arguments, so info must not be
         * negative). */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "zgeqrf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgeqrf)(layout, M, N, a, LD, tau) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* The high level allocates the workspace (column-major: the only
     * allocation). */
    LAPACKE_ZGEQRF_ALLOC_TEST(0, 0, "zgeqrf work alloc failure",
                              LAPACK_WORK_MEMORY_ERROR);

    /* Row-major allocates the workspace, then the transposed copy of A. */
    LAPACKE_ZGEQRF_ALLOC_TEST(1, 0, "zgeqrf work alloc failure",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGEQRF_ALLOC_TEST(1, 1, "zgeqrf transpose alloc failure (a)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Recovery, scheduled one past the last allocation. */
    LAPACKE_ZGEQRF_ALLOC_TEST(1, 2, "zgeqrf recovers after alloc failure", 0);
    lapacke_test_check_alloc_count("zgeqrf allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_ZGEQRF_ALLOC_TEST(2, 0, "zgeqrf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgeqrf invalid layout allocation count");
}
