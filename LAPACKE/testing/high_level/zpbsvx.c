#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)
#define REGION_AFB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AFB_L(i, j) lapacke_test_region_band(i, j, N, 0)

/* Refill the inputs, schedule the malloc failure, call zpbsvx. */
#define LAPACKE_ZPBSVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_zfill_pb(layout, 'U', N, KD, afb, LD);                    \
        equed = 'N';                                                           \
        lapacke_test_dfill_pos(LD * LD, s);                                    \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_zfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zpbsvx)(                         \
                               layout, 'N', 'U', N, KD, NRHS, ab, LD, afb, LD, \
                               &equed, s, b, LD, x, LD, &rcond, ferr, berr),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zpbsvx)
{
    lapack_complex_double ab[LD * LD];
    lapack_complex_double afb[LD * LD];
    char equed;
    double s[LD * LD];
    lapack_complex_double b[LD * LD];
    lapack_complex_double x[LD * LD];
    double rcond;
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx ab fact=N uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -7,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'N', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx b fact=N uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -13,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'N', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx ab fact=N uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -7,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'N', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx b fact=N uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -13,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'N', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx ab fact=F uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -7,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpbtrf)(layout, 'U', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'F', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx afb fact=F uplo=U", l, KD + 1, N, afb, LD, REGION_AFB_U, -9,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpbtrf)(layout, 'U', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'F', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx b fact=F uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -13,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpbtrf)(layout, 'U', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'F', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "zpbsvx s fact=F uplo=U", l, N, 1, s, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -12,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'U', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpbtrf)(layout, 'U', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'F', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx ab fact=F uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -7,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpbtrf)(layout, 'L', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'F', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx afb fact=F uplo=L", l, KD + 1, N, afb, LD, REGION_AFB_L, -9,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpbtrf)(layout, 'L', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'F', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbsvx b fact=F uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -13,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpbtrf)(layout, 'L', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'F', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "zpbsvx s fact=F uplo=L", l, N, 1, s, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -12,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_zfill_pb(layout, 'L', N, KD, afb, LD), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpbtrf)(layout, 'L', N, KD, afb, LD)),
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'F', 'L', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        equed = 'N';
        lapacke_test_zfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_zfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_zfill_nan(layout, KD + 1, N, afb, LD);
        lapacke_test_dfill_nan(layout, N, 1, s, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_check(
            "zpbsvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zpbsvx)(layout, 'N', 'U', N, KD, NRHS, ab, LD,
                                       afb, LD, &equed, s, b, LD, x, LD, &rcond,
                                       ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZPBSVX_ALLOC_TEST(0, 0, "zpbsvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPBSVX_ALLOC_TEST(0, 1, "zpbsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPBSVX_ALLOC_TEST(0, 2, "zpbsvx allocation count", 0);
    lapacke_test_check_alloc_count("zpbsvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPBSVX_ALLOC_TEST(1, 0, "zpbsvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPBSVX_ALLOC_TEST(1, 1, "zpbsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPBSVX_ALLOC_TEST(1, 2, "zpbsvx transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPBSVX_ALLOC_TEST(1, 3, "zpbsvx transpose alloc failure (afb_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPBSVX_ALLOC_TEST(1, 4, "zpbsvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPBSVX_ALLOC_TEST(1, 5, "zpbsvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPBSVX_ALLOC_TEST(1, 6, "zpbsvx allocation count", 0);
    lapacke_test_check_alloc_count("zpbsvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPBSVX_ALLOC_TEST(2, 0, "zpbsvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zpbsvx invalid layout allocation count");
}
