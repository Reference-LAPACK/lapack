#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* zlange is a value-returning function; on NaN detection it returns the
 * (negative) info code as its value. Map a plausible norm (finite,
 * nonnegative) to 0 so the sweep can compare info codes. A NaN result
 * means the routine consumed the swept NaN instead of rejecting it; map
 * it to a sentinel instead of casting, which is undefined behavior for a
 * NaN and yields 0 on AArch64, masking the failure. */
static lapack_int zlange_info(double value)
{
    if (value != value) {
        return -999;
    }
    if (value >= 0.0) {
        return 0;
    }
    return (lapack_int)value;
}

/* Refill the inputs, schedule the malloc failure, call zlange. */
#define LAPACKE_ZLANGE_ALLOC_TEST(layout_index, countdown, norm, name,         \
                                  expected)                                    \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, M, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           zlange_info(API_SUFFIX(LAPACKE_zlange)(             \
                               layout, norm, M, N, a, LD)),                    \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zlange)
{
    lapack_complex_double a[LD * LD];
    double res;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zlange a", l, M, N, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_zfill(layout, M, N, a, LD)),
            zlange_info(API_SUFFIX(LAPACKE_zlange)(layout, '1', M, N, a, LD)));

        /* NaN checks off: the NaN must not be rejected as an error. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, M, N, a, LD);
        res = API_SUFFIX(LAPACKE_zlange)(layout, '1', M, N, a, LD);
        lapacke_test_check("zlange NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           zlange_info(res) == -5 ? -5 : 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* On allocation failure the norm functions return 0.0, indistinguishable
     * from a zero norm; pinned as the current behavior. */
    lapacke_test_zfill(LAPACK_COL_MAJOR, M, N, a, LD);
    lapacke_test_schedule_malloc_failure(0);
    res = API_SUFFIX(LAPACKE_zlange)(LAPACK_COL_MAJOR, 'I', M, N, a, LD);
    lapacke_test_check("zlange alloc failure (work) returns 0.0", "col-major",
                       res == 0 ? 0 : -999, 0);
    LAPACKE_ZLANGE_ALLOC_TEST(0, 1, 'I', "zlange allocation count", 0);
    lapacke_test_check_alloc_count("zlange col-major allocation count");

    lapacke_test_zfill(LAPACK_ROW_MAJOR, M, N, a, LD);
    lapacke_test_schedule_malloc_failure(0);
    res = API_SUFFIX(LAPACKE_zlange)(LAPACK_ROW_MAJOR, 'I', M, N, a, LD);
    lapacke_test_check("zlange alloc failure (work) returns 0.0", "row-major",
                       res == 0 ? 0 : -999, 0);
    LAPACKE_ZLANGE_ALLOC_TEST(1, 1, 'I', "zlange allocation count", 0);
    lapacke_test_check_alloc_count("zlange row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZLANGE_ALLOC_TEST(2, 0, '1', "zlange invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zlange invalid layout allocation count");
}
