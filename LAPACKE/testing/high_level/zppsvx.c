#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zppsvx. */
#define LAPACKE_ZPPSVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_pp(layout, 'U', N, ap);                             \
        lapacke_test_zfill_pp(layout, 'U', N, afp);                            \
        equed = 'N';                                                           \
        lapacke_test_dfill_pos(LD * LD, s);                                    \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_zfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zppsvx)(                         \
                               layout, 'N', 'U', N, NRHS, ap, afp, &equed, s,  \
                               b, LD, x, LD, &rcond, ferr, berr),              \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zppsvx)
{
    lapack_complex_double ap[LD * LD];
    lapack_complex_double afp[LD * LD];
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
            "zppsvx ap fact=N uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_zfill_pp(layout, 'U', N, ap),
             lapacke_test_zfill_pp(layout, 'U', N, afp), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'N', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppsvx b fact=N uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_zfill_pp(layout, 'U', N, ap),
             lapacke_test_zfill_pp(layout, 'U', N, afp), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'N', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppsvx ap fact=N uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_zfill_pp(layout, 'L', N, ap),
             lapacke_test_zfill_pp(layout, 'L', N, afp), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'N', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppsvx b fact=N uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_zfill_pp(layout, 'L', N, ap),
             lapacke_test_zfill_pp(layout, 'L', N, afp), (equed = 'N'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'N', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppsvx afp fact=F uplo=U", l, N * (N + 1) / 2, 1, afp,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -7,
            (lapacke_test_zfill_pp(layout, 'U', N, ap),
             lapacke_test_zfill_pp(layout, 'U', N, afp), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpptrf)(layout, 'U', N, afp)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'F', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppsvx ap fact=F uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_zfill_pp(layout, 'U', N, ap),
             lapacke_test_zfill_pp(layout, 'U', N, afp), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpptrf)(layout, 'U', N, afp)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'F', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppsvx b fact=F uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_zfill_pp(layout, 'U', N, ap),
             lapacke_test_zfill_pp(layout, 'U', N, afp), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpptrf)(layout, 'U', N, afp)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'F', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "zppsvx s fact=F uplo=U", l, N, 1, s, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -9,
            (lapacke_test_zfill_pp(layout, 'U', N, ap),
             lapacke_test_zfill_pp(layout, 'U', N, afp), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpptrf)(layout, 'U', N, afp)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'F', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppsvx afp fact=F uplo=L", l, N * (N + 1) / 2, 1, afp,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -7,
            (lapacke_test_zfill_pp(layout, 'L', N, ap),
             lapacke_test_zfill_pp(layout, 'L', N, afp), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpptrf)(layout, 'L', N, afp)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'F', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppsvx ap fact=F uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_zfill_pp(layout, 'L', N, ap),
             lapacke_test_zfill_pp(layout, 'L', N, afp), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpptrf)(layout, 'L', N, afp)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'F', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zppsvx b fact=F uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_zfill_pp(layout, 'L', N, ap),
             lapacke_test_zfill_pp(layout, 'L', N, afp), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpptrf)(layout, 'L', N, afp)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'F', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "zppsvx s fact=F uplo=L", l, N, 1, s, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -9,
            (lapacke_test_zfill_pp(layout, 'L', N, ap),
             lapacke_test_zfill_pp(layout, 'L', N, afp), (equed = 'Y'),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpptrf)(layout, 'L', N, afp)),
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'F', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        equed = 'N';
        lapacke_test_zfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_zfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_zfill_nan(layout, N * (N + 1) / 2, 1, afp,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_dfill_nan(layout, N, 1, s, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_check(
            "zppsvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zppsvx)(layout, 'N', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZPPSVX_ALLOC_TEST(0, 0, "zppsvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPPSVX_ALLOC_TEST(0, 1, "zppsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPPSVX_ALLOC_TEST(0, 2, "zppsvx allocation count", 0);
    lapacke_test_check_alloc_count("zppsvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPPSVX_ALLOC_TEST(1, 0, "zppsvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPPSVX_ALLOC_TEST(1, 1, "zppsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPPSVX_ALLOC_TEST(1, 2, "zppsvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPPSVX_ALLOC_TEST(1, 3, "zppsvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPPSVX_ALLOC_TEST(1, 4, "zppsvx transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPPSVX_ALLOC_TEST(1, 5, "zppsvx transpose alloc failure (afp_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPPSVX_ALLOC_TEST(1, 6, "zppsvx allocation count", 0);
    lapacke_test_check_alloc_count("zppsvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPPSVX_ALLOC_TEST(2, 0, "zppsvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zppsvx invalid layout allocation count");
}
