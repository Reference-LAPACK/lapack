#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)
#define REGION_AFB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AFB_L(i, j) lapacke_test_region_band(i, j, N, 0)

/* Refill the inputs, schedule the malloc failure, call spbsvx. */
#define LAPACKE_SPBSVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_sfill_pb(layout, 'U', N, KD, afb, LD);                    \
        equed = 'N';                                                           \
        lapacke_test_sfill_pos(LD * LD, s);                                    \
        lapacke_test_sfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_spbsvx)(                         \
                               layout, 'N', 'U', N, KD, NRHS, ab, LD, afb, LD, \
                               &equed, s, b, LD, x, LD, &rcond, ferr, berr),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(spbsvx)
{
    float ab[LD * LD];
    float afb[LD * LD];
    char equed;
    float s[LD * LD];
    float b[LD * LD];
    float x[LD * LD];
    float rcond;
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx ab fact=N uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -7,
            (lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'N', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx b fact=N uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -13,
            (lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'N', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx ab fact=N uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -7,
            (lapacke_test_sfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'N', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx b fact=N uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -13,
            (lapacke_test_sfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'N', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx ab fact=F uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -7,
            (lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spbtrf)(layout, 'U', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'F', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx afb fact=F uplo=U", l, KD + 1, N, afb, LD, REGION_AFB_U, -9,
            (lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spbtrf)(layout, 'U', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'F', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx b fact=F uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -13,
            (lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spbtrf)(layout, 'U', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'F', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx s fact=F uplo=U", l, N, 1, s, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -12,
            (lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spbtrf)(layout, 'U', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'F', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx ab fact=F uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -7,
            (lapacke_test_sfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spbtrf)(layout, 'L', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'F', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx afb fact=F uplo=L", l, KD + 1, N, afb, LD, REGION_AFB_L, -9,
            (lapacke_test_sfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spbtrf)(layout, 'L', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'F', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx b fact=F uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -13,
            (lapacke_test_sfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spbtrf)(layout, 'L', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'F', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbsvx s fact=F uplo=L", l, N, 1, s, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -12,
            (lapacke_test_sfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spbtrf)(layout, 'L', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'F', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        equed = 'N';
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_sfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_sfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_sfill_nan(layout, KD + 1, N, afb, LD);
        lapacke_test_sfill_nan(layout, N, 1, s, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_check(
            "spbsvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_spbsvx)(layout, 'N', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SPBSVX_ALLOC_TEST(0, 0, "spbsvx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPBSVX_ALLOC_TEST(0, 1, "spbsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPBSVX_ALLOC_TEST(0, 2, "spbsvx allocation count", 0);
    lapacke_test_check_alloc_count("spbsvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPBSVX_ALLOC_TEST(1, 0, "spbsvx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPBSVX_ALLOC_TEST(1, 1, "spbsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPBSVX_ALLOC_TEST(1, 2, "spbsvx transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPBSVX_ALLOC_TEST(1, 3, "spbsvx transpose alloc failure (afb_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPBSVX_ALLOC_TEST(1, 4, "spbsvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPBSVX_ALLOC_TEST(1, 5, "spbsvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPBSVX_ALLOC_TEST(1, 6, "spbsvx allocation count", 0);
    lapacke_test_check_alloc_count("spbsvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPBSVX_ALLOC_TEST(2, 0, "spbsvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("spbsvx invalid layout allocation count");
}
