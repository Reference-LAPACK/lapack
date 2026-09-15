#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sppsvx. */
#define LAPACKE_SPPSVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_pp(layout, 'U', N, ap);                             \
        lapacke_test_sfill_pp(layout, 'U', N, afp);                            \
        equed = 'N';                                                           \
        lapacke_test_sfill_pos(LD * LD, s);                                    \
        lapacke_test_sfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sppsvx)(                         \
                               layout, 'N', 'U', N, NRHS, ap, afp, &equed, s,  \
                               b, LD, x, LD, &rcond, ferr, berr),              \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sppsvx)
{
    float ap[LD * LD];
    float afp[LD * LD];
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
            "sppsvx ap fact=N uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_sfill_pp(layout, 'U', N, ap),
             lapacke_test_sfill_pp(layout, 'U', N, afp), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'N', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx b fact=N uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_sfill_pp(layout, 'U', N, ap),
             lapacke_test_sfill_pp(layout, 'U', N, afp), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'N', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx ap fact=N uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_sfill_pp(layout, 'L', N, ap),
             lapacke_test_sfill_pp(layout, 'L', N, afp), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'N', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx b fact=N uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_sfill_pp(layout, 'L', N, ap),
             lapacke_test_sfill_pp(layout, 'L', N, afp), (equed = 'N'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'N', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx afp fact=F uplo=U", l, N * (N + 1) / 2, 1, afp,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -7,
            (lapacke_test_sfill_pp(layout, 'U', N, ap),
             lapacke_test_sfill_pp(layout, 'U', N, afp), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spptrf)(layout, 'U', N, afp)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'F', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx ap fact=F uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_sfill_pp(layout, 'U', N, ap),
             lapacke_test_sfill_pp(layout, 'U', N, afp), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spptrf)(layout, 'U', N, afp)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'F', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx b fact=F uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_sfill_pp(layout, 'U', N, ap),
             lapacke_test_sfill_pp(layout, 'U', N, afp), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spptrf)(layout, 'U', N, afp)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'F', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx s fact=F uplo=U", l, N, 1, s, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -9,
            (lapacke_test_sfill_pp(layout, 'U', N, ap),
             lapacke_test_sfill_pp(layout, 'U', N, afp), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spptrf)(layout, 'U', N, afp)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'F', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx afp fact=F uplo=L", l, N * (N + 1) / 2, 1, afp,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -7,
            (lapacke_test_sfill_pp(layout, 'L', N, ap),
             lapacke_test_sfill_pp(layout, 'L', N, afp), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spptrf)(layout, 'L', N, afp)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'F', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx ap fact=F uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_sfill_pp(layout, 'L', N, ap),
             lapacke_test_sfill_pp(layout, 'L', N, afp), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spptrf)(layout, 'L', N, afp)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'F', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx b fact=F uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_sfill_pp(layout, 'L', N, ap),
             lapacke_test_sfill_pp(layout, 'L', N, afp), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spptrf)(layout, 'L', N, afp)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'F', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppsvx s fact=F uplo=L", l, N, 1, s, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -9,
            (lapacke_test_sfill_pp(layout, 'L', N, ap),
             lapacke_test_sfill_pp(layout, 'L', N, afp), (equed = 'Y'),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spptrf)(layout, 'L', N, afp)),
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'F', 'L', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        equed = 'N';
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_sfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_sfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_sfill_nan(layout, N * (N + 1) / 2, 1, afp,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_sfill_nan(layout, N, 1, s, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_check(
            "sppsvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sppsvx)(layout, 'N', 'U', N, NRHS, ap, afp,
                                       &equed, s, b, LD, x, LD, &rcond, ferr,
                                       berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SPPSVX_ALLOC_TEST(0, 0, "sppsvx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPPSVX_ALLOC_TEST(0, 1, "sppsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPPSVX_ALLOC_TEST(0, 2, "sppsvx allocation count", 0);
    lapacke_test_check_alloc_count("sppsvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPPSVX_ALLOC_TEST(1, 0, "sppsvx work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPPSVX_ALLOC_TEST(1, 1, "sppsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPPSVX_ALLOC_TEST(1, 2, "sppsvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPPSVX_ALLOC_TEST(1, 3, "sppsvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPPSVX_ALLOC_TEST(1, 4, "sppsvx transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPPSVX_ALLOC_TEST(1, 5, "sppsvx transpose alloc failure (afp_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPPSVX_ALLOC_TEST(1, 6, "sppsvx allocation count", 0);
    lapacke_test_check_alloc_count("sppsvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPPSVX_ALLOC_TEST(2, 0, "sppsvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sppsvx invalid layout allocation count");
}
