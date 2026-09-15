#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* clanhe is a value-returning function; on NaN detection it returns the
 * (negative) info code as its value. Map a plausible norm (finite,
 * nonnegative) to 0 so the sweep can compare info codes. A NaN result
 * means the routine consumed the swept NaN instead of rejecting it; map
 * it to a sentinel instead of casting, which is undefined behavior for a
 * NaN and yields 0 on AArch64, masking the failure. */
static lapack_int clanhe_info(float value)
{
    if (value != value) {
        return -999;
    }
    if (value >= 0.0) {
        return 0;
    }
    return (lapack_int)value;
}

/* Refill the inputs, schedule the malloc failure, call clanhe. */
#define LAPACKE_CLANHE_ALLOC_TEST(layout_index, countdown, norm, name,         \
                                  expected)                                    \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_sym(layout, N, a, LD);                              \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           clanhe_info(API_SUFFIX(LAPACKE_clanhe)(             \
                               layout, norm, 'U', N, a, LD)),                  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(clanhe)
{
    lapack_complex_float a[LD * LD];
    float res;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP("clanhe a uplo=U", l, N, N, a, LD,
                                lapacke_test_region_upper, -5,
                                (lapacke_test_cfill_sym(layout, N, a, LD)),
                                clanhe_info(API_SUFFIX(LAPACKE_clanhe)(
                                    layout, '1', 'U', N, a, LD)));

        LAPACKE_TEST_CNAN_SWEEP("clanhe a uplo=L", l, N, N, a, LD,
                                lapacke_test_region_lower, -5,
                                (lapacke_test_cfill_sym(layout, N, a, LD)),
                                clanhe_info(API_SUFFIX(LAPACKE_clanhe)(
                                    layout, '1', 'L', N, a, LD)));

        /* NaN checks off: the NaN must not be rejected as an error. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        res = API_SUFFIX(LAPACKE_clanhe)(layout, '1', 'U', N, a, LD);
        lapacke_test_check("clanhe NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           clanhe_info(res) == -5 ? -5 : 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* On allocation failure the norm functions return 0.0, indistinguishable
     * from a zero norm; pinned as the current behavior. */
    lapacke_test_cfill_sym(LAPACK_COL_MAJOR, N, a, LD);
    lapacke_test_schedule_malloc_failure(0);
    res = API_SUFFIX(LAPACKE_clanhe)(LAPACK_COL_MAJOR, 'I', 'U', N, a, LD);
    lapacke_test_check("clanhe alloc failure (work) returns 0.0", "col-major",
                       res == 0 ? 0 : -999, 0);
    LAPACKE_CLANHE_ALLOC_TEST(0, 1, 'I', "clanhe allocation count", 0);
    lapacke_test_check_alloc_count("clanhe col-major allocation count");

    lapacke_test_cfill_sym(LAPACK_ROW_MAJOR, N, a, LD);
    lapacke_test_schedule_malloc_failure(0);
    res = API_SUFFIX(LAPACKE_clanhe)(LAPACK_ROW_MAJOR, 'I', 'U', N, a, LD);
    lapacke_test_check("clanhe alloc failure (work) returns 0.0", "row-major",
                       res == 0 ? 0 : -999, 0);
    lapacke_test_cfill_sym(LAPACK_ROW_MAJOR, N, a, LD);
    lapacke_test_schedule_malloc_failure(1);
    res = API_SUFFIX(LAPACKE_clanhe)(LAPACK_ROW_MAJOR, 'I', 'U', N, a, LD);
    lapacke_test_check("clanhe alloc failure (a_t) returns 0.0", "row-major",
                       res == 0 ? 0 : -999, 0);
    LAPACKE_CLANHE_ALLOC_TEST(1, 2, 'I', "clanhe allocation count", 0);
    lapacke_test_check_alloc_count("clanhe row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CLANHE_ALLOC_TEST(2, 0, '1', "clanhe invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("clanhe invalid layout allocation count");
}
