#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* dlange is a value-returning function; on NaN detection it returns the
 * (negative) info code as its value. Map a plausible norm (finite,
 * nonnegative) to 0 so the sweep can compare info codes. A NaN result
 * means the routine consumed the swept NaN instead of rejecting it; map
 * it to a sentinel instead of casting (casting NaN to an integer is
 * undefined behavior and yields 0 on AArch64, masking the failure). */
static lapack_int dlange_info(double value)
{
    if (value != value) {
        return -999;
    }
    if (value >= 0.0) {
        return 0;
    }
    return (lapack_int)value;
}

/* One allocation test: refill the input, schedule the malloc failure
 * countdown and check the mapped result of a dlange call with the given
 * norm in the indexed layout. */
#define LAPACKE_DLANGE_ALLOC_TEST(layout_index, countdown, norm, name,         \
                                  expected)                                    \
    do {                                                                       \
        lapacke_test_fill(lapacke_test_layouts[layout_index], M, N, a, LD);    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            dlange_info(LAPACKE_dlange(lapacke_test_layouts[layout_index],     \
                                       norm, M, N, a, LD)),                    \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(dlange)
{
    double a[LD * LD];
    double res;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];
        LAPACKE_TEST_NAN_SWEEP(
            "dlange a", l, M, N, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_fill(layout, M, N, a, LD)),
            dlange_info(LAPACKE_dlange(layout, '1', M, N, a, LD)));

        /* With NaN checking disabled the NaN must not be rejected with
         * -5. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_nan(layout, M, N, a, LD);
        res = LAPACKE_dlange(layout, '1', M, N, a, LD);
        lapacke_test_check("dlange NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           dlange_info(res) == -5 ? -5 : 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* The infinity norm allocates a work array in the high level. On
     * allocation failure LAPACKE_dlange returns 0.0 -- indistinguishable
     * from a genuine zero norm; checked exactly here (the macro's mapping
     * would accept any plausible norm). Pinned as the current behavior;
     * questionable enough to raise upstream. */
    lapacke_test_fill(LAPACK_COL_MAJOR, M, N, a, LD);
    lapacke_test_schedule_malloc_failure(0);
    res = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', M, N, a, LD);
    lapacke_test_check("dlange work alloc failure returns 0.0", "col-major",
                       res == 0.0 ? 0 : -999, 0);

    /* Recovery, scheduled one past the last allocation. */
    LAPACKE_DLANGE_ALLOC_TEST(0, 1, 'I', "dlange recovers after alloc failure",
                              0);
    lapacke_test_check_alloc_count("dlange allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_DLANGE_ALLOC_TEST(2, 0, '1', "dlange invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dlange invalid layout allocation count");
}
