#ifndef CBLAS_64_H
#define CBLAS_64_H
#include <stddef.h>
#include <stdint.h>
#include <inttypes.h>

#include "cblas.h"

#ifdef __cplusplus
extern "C" {            /* Assume C declarations for C++ */
#endif /* __cplusplus */

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */

double cblas_dcabs1_64(const void  *z);
float  cblas_scabs1_64(const void  *c);

float  cblas_sdsdot_64(const int64_t N, const float alpha, const float *X,
                    const int64_t incX, const float *Y, const int64_t incY);
double cblas_dsdot_64(const int64_t N, const float *X, const int64_t incX, const float *Y,
                   const int64_t incY);
float  cblas_sdot_64(const int64_t N, const float  *X, const int64_t incX,
                  const float  *Y, const int64_t incY);
double cblas_ddot_64(const int64_t N, const double *X, const int64_t incX,
                  const double *Y, const int64_t incY);

/*
 * Functions having prefixes Z and C only
 */
void   cblas_cdotu_sub_64(const int64_t N, const void *X, const int64_t incX,
                       const void *Y, const int64_t incY, void *dotu);
void   cblas_cdotc_sub_64(const int64_t N, const void *X, const int64_t incX,
                       const void *Y, const int64_t incY, void *dotc);

void   cblas_zdotu_sub_64(const int64_t N, const void *X, const int64_t incX,
                       const void *Y, const int64_t incY, void *dotu);
void   cblas_zdotc_sub_64(const int64_t N, const void *X, const int64_t incX,
                       const void *Y, const int64_t incY, void *dotc);


/*
 * Functions having prefixes S D SC DZ
 */
float  cblas_snrm2_64(const int64_t N, const float *X, const int64_t incX);
float  cblas_sasum_64(const int64_t N, const float *X, const int64_t incX);

double cblas_dnrm2_64(const int64_t N, const double *X, const int64_t incX);
double cblas_dasum_64(const int64_t N, const double *X, const int64_t incX);

float  cblas_scnrm2_64(const int64_t N, const void *X, const int64_t incX);
float  cblas_scasum_64(const int64_t N, const void *X, const int64_t incX);

double cblas_dznrm2_64(const int64_t N, const void *X, const int64_t incX);
double cblas_dzasum_64(const int64_t N, const void *X, const int64_t incX);


/*
 * Functions having standard 4 prefixes (S D C Z)
 */
CBLAS_INDEX cblas_isamax_64(const int64_t N, const float  *X, const int64_t incX);
CBLAS_INDEX cblas_idamax_64(const int64_t N, const double *X, const int64_t incX);
CBLAS_INDEX cblas_icamax_64(const int64_t N, const void   *X, const int64_t incX);
CBLAS_INDEX cblas_izamax_64(const int64_t N, const void   *X, const int64_t incX);

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void cblas_sswap_64(const int64_t N, float *X, const int64_t incX,
                 float *Y, const int64_t incY);
void cblas_scopy_64(const int64_t N, const float *X, const int64_t incX,
                 float *Y, const int64_t incY);
void cblas_saxpy_64(const int64_t N, const float alpha, const float *X,
                 const int64_t incX, float *Y, const int64_t incY);

void cblas_dswap_64(const int64_t N, double *X, const int64_t incX,
                 double *Y, const int64_t incY);
void cblas_dcopy_64(const int64_t N, const double *X, const int64_t incX,
                 double *Y, const int64_t incY);
void cblas_daxpy_64(const int64_t N, const double alpha, const double *X,
                 const int64_t incX, double *Y, const int64_t incY);

void cblas_cswap_64(const int64_t N, void *X, const int64_t incX,
                 void *Y, const int64_t incY);
void cblas_ccopy_64(const int64_t N, const void *X, const int64_t incX,
                 void *Y, const int64_t incY);
void cblas_caxpy_64(const int64_t N, const void *alpha, const void *X,
                 const int64_t incX, void *Y, const int64_t incY);

void cblas_zswap_64(const int64_t N, void *X, const int64_t incX,
                 void *Y, const int64_t incY);
void cblas_zcopy_64(const int64_t N, const void *X, const int64_t incX,
                 void *Y, const int64_t incY);
void cblas_zaxpy_64(const int64_t N, const void *alpha, const void *X,
                 const int64_t incX, void *Y, const int64_t incY);


/*
 * Routines with S and D prefix only
 */
void cblas_srotmg_64(float *d1, float *d2, float *b1, const float b2, float *P);
void cblas_srotm_64(const int64_t N, float *X, const int64_t incX,
                 float *Y, const int64_t incY, const float *P);
void cblas_drotmg_64(double *d1, double *d2, double *b1, const double b2, double *P);
void cblas_drotm_64(const int64_t N, double *X, const int64_t incX,
                 double *Y, const int64_t incY, const double *P);



/*
 * Routines with S D C Z CS and ZD prefixes
 */
void cblas_sscal_64(const int64_t N, const float alpha, float *X, const int64_t incX);
void cblas_dscal_64(const int64_t N, const double alpha, double *X, const int64_t incX);
void cblas_cscal_64(const int64_t N, const void *alpha, void *X, const int64_t incX);
void cblas_zscal_64(const int64_t N, const void *alpha, void *X, const int64_t incX);
void cblas_csscal_64(const int64_t N, const float alpha, void *X, const int64_t incX);
void cblas_zdscal_64(const int64_t N, const double alpha, void *X, const int64_t incX);

void cblas_srotg_64(float *a, float *b, float *c, float *s);
void cblas_drotg_64(double *a, double *b, double *c, double *s);
void cblas_crotg_64(void *a, void *b, float *c, void *s);
void cblas_zrotg_64(void *a, void *b, double *c, void *s);

void cblas_srot_64(const int64_t N, float *X, const int64_t incX,
                float *Y, const int64_t incY, const float c, const float s);
void cblas_drot_64(const int64_t N, double *X, const int64_t incX,
                double *Y, const int64_t incY, const double c, const double  s);
void cblas_csrot_64(const int64_t N, void *X, const int64_t incX,
                 void *Y, const int64_t incY, const float c, const float s);
void cblas_zdrot_64(const int64_t N, void *X, const int64_t incX,
                 void *Y, const int64_t incY, const double c, const double s);

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemv_64(const CBLAS_LAYOUT layout,
                 const CBLAS_TRANSPOSE TransA, const int64_t M, const int64_t N,
                 const float alpha, const float *A, const int64_t lda,
                 const float *X, const int64_t incX, const float beta,
                 float *Y, const int64_t incY);
void cblas_sgbmv_64(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int64_t M, const int64_t N,
                 const int64_t KL, const int64_t KU, const float alpha,
                 const float *A, const int64_t lda, const float *X,
                 const int64_t incX, const float beta, float *Y, const int64_t incY);
void cblas_strmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const float *A, const int64_t lda,
                 float *X, const int64_t incX);
void cblas_stbmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const int64_t K, const float *A, const int64_t lda,
                 float *X, const int64_t incX);
void cblas_stpmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const float *Ap, float *X, const int64_t incX);
void cblas_strsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const float *A, const int64_t lda, float *X,
                 const int64_t incX);
void cblas_stbsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const int64_t K, const float *A, const int64_t lda,
                 float *X, const int64_t incX);
void cblas_stpsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const float *Ap, float *X, const int64_t incX);

void cblas_dgemv_64(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int64_t M, const int64_t N,
                 const double alpha, const double *A, const int64_t lda,
                 const double *X, const int64_t incX, const double beta,
                 double *Y, const int64_t incY);
void cblas_dgbmv_64(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int64_t M, const int64_t N,
                 const int64_t KL, const int64_t KU, const double alpha,
                 const double *A, const int64_t lda, const double *X,
                 const int64_t incX, const double beta, double *Y, const int64_t incY);
void cblas_dtrmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const double *A, const int64_t lda,
                 double *X, const int64_t incX);
void cblas_dtbmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const int64_t K, const double *A, const int64_t lda,
                 double *X, const int64_t incX);
void cblas_dtpmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const double *Ap, double *X, const int64_t incX);
void cblas_dtrsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const double *A, const int64_t lda, double *X,
                 const int64_t incX);
void cblas_dtbsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const int64_t K, const double *A, const int64_t lda,
                 double *X, const int64_t incX);
void cblas_dtpsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const double *Ap, double *X, const int64_t incX);

void cblas_cgemv_64(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 const void *X, const int64_t incX, const void *beta,
                 void *Y, const int64_t incY);
void cblas_cgbmv_64(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int64_t M, const int64_t N,
                 const int64_t KL, const int64_t KU, const void *alpha,
                 const void *A, const int64_t lda, const void *X,
                 const int64_t incX, const void *beta, void *Y, const int64_t incY);
void cblas_ctrmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const void *A, const int64_t lda,
                 void *X, const int64_t incX);
void cblas_ctbmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const int64_t K, const void *A, const int64_t lda,
                 void *X, const int64_t incX);
void cblas_ctpmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const void *Ap, void *X, const int64_t incX);
void cblas_ctrsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const void *A, const int64_t lda, void *X,
                 const int64_t incX);
void cblas_ctbsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const int64_t K, const void *A, const int64_t lda,
                 void *X, const int64_t incX);
void cblas_ctpsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const void *Ap, void *X, const int64_t incX);

void cblas_zgemv_64(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 const void *X, const int64_t incX, const void *beta,
                 void *Y, const int64_t incY);
void cblas_zgbmv_64(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int64_t M, const int64_t N,
                 const int64_t KL, const int64_t KU, const void *alpha,
                 const void *A, const int64_t lda, const void *X,
                 const int64_t incX, const void *beta, void *Y, const int64_t incY);
void cblas_ztrmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const void *A, const int64_t lda,
                 void *X, const int64_t incX);
void cblas_ztbmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const int64_t K, const void *A, const int64_t lda,
                 void *X, const int64_t incX);
void cblas_ztpmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const void *Ap, void *X, const int64_t incX);
void cblas_ztrsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const void *A, const int64_t lda, void *X,
                 const int64_t incX);
void cblas_ztbsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const int64_t K, const void *A, const int64_t lda,
                 void *X, const int64_t incX);
void cblas_ztpsv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int64_t N, const void *Ap, void *X, const int64_t incX);


/*
 * Routines with S and D prefixes only
 */
void cblas_ssymv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const float alpha, const float *A,
                 const int64_t lda, const float *X, const int64_t incX,
                 const float beta, float *Y, const int64_t incY);
void cblas_ssbmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const int64_t K, const float alpha, const float *A,
                 const int64_t lda, const float *X, const int64_t incX,
                 const float beta, float *Y, const int64_t incY);
void cblas_sspmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const float alpha, const float *Ap,
                 const float *X, const int64_t incX,
                 const float beta, float *Y, const int64_t incY);
void cblas_sger_64(CBLAS_LAYOUT layout, const int64_t M, const int64_t N,
                const float alpha, const float *X, const int64_t incX,
                const float *Y, const int64_t incY, float *A, const int64_t lda);
void cblas_ssyr_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const float alpha, const float *X,
                const int64_t incX, float *A, const int64_t lda);
void cblas_sspr_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const float alpha, const float *X,
                const int64_t incX, float *Ap);
void cblas_ssyr2_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const float alpha, const float *X,
                const int64_t incX, const float *Y, const int64_t incY, float *A,
                const int64_t lda);
void cblas_sspr2_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const float alpha, const float *X,
                const int64_t incX, const float *Y, const int64_t incY, float *A);

void cblas_dsymv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const double alpha, const double *A,
                 const int64_t lda, const double *X, const int64_t incX,
                 const double beta, double *Y, const int64_t incY);
void cblas_dsbmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const int64_t K, const double alpha, const double *A,
                 const int64_t lda, const double *X, const int64_t incX,
                 const double beta, double *Y, const int64_t incY);
void cblas_dspmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const double alpha, const double *Ap,
                 const double *X, const int64_t incX,
                 const double beta, double *Y, const int64_t incY);
void cblas_dger_64(CBLAS_LAYOUT layout, const int64_t M, const int64_t N,
                const double alpha, const double *X, const int64_t incX,
                const double *Y, const int64_t incY, double *A, const int64_t lda);
void cblas_dsyr_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const double alpha, const double *X,
                const int64_t incX, double *A, const int64_t lda);
void cblas_dspr_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const double alpha, const double *X,
                const int64_t incX, double *Ap);
void cblas_dsyr2_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const double alpha, const double *X,
                const int64_t incX, const double *Y, const int64_t incY, double *A,
                const int64_t lda);
void cblas_dspr2_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const double alpha, const double *X,
                const int64_t incX, const double *Y, const int64_t incY, double *A);


/*
 * Routines with C and Z prefixes only
 */
void cblas_chemv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const void *alpha, const void *A,
                 const int64_t lda, const void *X, const int64_t incX,
                 const void *beta, void *Y, const int64_t incY);
void cblas_chbmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const int64_t K, const void *alpha, const void *A,
                 const int64_t lda, const void *X, const int64_t incX,
                 const void *beta, void *Y, const int64_t incY);
void cblas_chpmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const void *alpha, const void *Ap,
                 const void *X, const int64_t incX,
                 const void *beta, void *Y, const int64_t incY);
void cblas_cgeru_64(CBLAS_LAYOUT layout, const int64_t M, const int64_t N,
                 const void *alpha, const void *X, const int64_t incX,
                 const void *Y, const int64_t incY, void *A, const int64_t lda);
void cblas_cgerc_64(CBLAS_LAYOUT layout, const int64_t M, const int64_t N,
                 const void *alpha, const void *X, const int64_t incX,
                 const void *Y, const int64_t incY, void *A, const int64_t lda);
void cblas_cher_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const float alpha, const void *X, const int64_t incX,
                void *A, const int64_t lda);
void cblas_chpr_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const float alpha, const void *X,
                const int64_t incX, void *A);
void cblas_cher2_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, const int64_t N,
                const void *alpha, const void *X, const int64_t incX,
                const void *Y, const int64_t incY, void *A, const int64_t lda);
void cblas_chpr2_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, const int64_t N,
                const void *alpha, const void *X, const int64_t incX,
                const void *Y, const int64_t incY, void *Ap);

void cblas_zhemv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const void *alpha, const void *A,
                 const int64_t lda, const void *X, const int64_t incX,
                 const void *beta, void *Y, const int64_t incY);
void cblas_zhbmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const int64_t K, const void *alpha, const void *A,
                 const int64_t lda, const void *X, const int64_t incX,
                 const void *beta, void *Y, const int64_t incY);
void cblas_zhpmv_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int64_t N, const void *alpha, const void *Ap,
                 const void *X, const int64_t incX,
                 const void *beta, void *Y, const int64_t incY);
void cblas_zgeru_64(CBLAS_LAYOUT layout, const int64_t M, const int64_t N,
                 const void *alpha, const void *X, const int64_t incX,
                 const void *Y, const int64_t incY, void *A, const int64_t lda);
void cblas_zgerc_64(CBLAS_LAYOUT layout, const int64_t M, const int64_t N,
                 const void *alpha, const void *X, const int64_t incX,
                 const void *Y, const int64_t incY, void *A, const int64_t lda);
void cblas_zher_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const double alpha, const void *X, const int64_t incX,
                void *A, const int64_t lda);
void cblas_zhpr_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int64_t N, const double alpha, const void *X,
                const int64_t incX, void *A);
void cblas_zher2_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, const int64_t N,
                const void *alpha, const void *X, const int64_t incX,
                const void *Y, const int64_t incY, void *A, const int64_t lda);
void cblas_zhpr2_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, const int64_t N,
                const void *alpha, const void *X, const int64_t incX,
                const void *Y, const int64_t incY, void *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemm_64(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int64_t M, const int64_t N,
                 const int64_t K, const float alpha, const float *A,
                 const int64_t lda, const float *B, const int64_t ldb,
                 const float beta, float *C, const int64_t ldc);
void cblas_ssymm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int64_t M, const int64_t N,
                 const float alpha, const float *A, const int64_t lda,
                 const float *B, const int64_t ldb, const float beta,
                 float *C, const int64_t ldc);
void cblas_ssyrk_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                 const float alpha, const float *A, const int64_t lda,
                 const float beta, float *C, const int64_t ldc);
void cblas_ssyr2k_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                  const float alpha, const float *A, const int64_t lda,
                  const float *B, const int64_t ldb, const float beta,
                  float *C, const int64_t ldc);
void cblas_strmm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int64_t M, const int64_t N,
                 const float alpha, const float *A, const int64_t lda,
                 float *B, const int64_t ldb);
void cblas_strsm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int64_t M, const int64_t N,
                 const float alpha, const float *A, const int64_t lda,
                 float *B, const int64_t ldb);

void cblas_dgemm_64(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int64_t M, const int64_t N,
                 const int64_t K, const double alpha, const double *A,
                 const int64_t lda, const double *B, const int64_t ldb,
                 const double beta, double *C, const int64_t ldc);
void cblas_dsymm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int64_t M, const int64_t N,
                 const double alpha, const double *A, const int64_t lda,
                 const double *B, const int64_t ldb, const double beta,
                 double *C, const int64_t ldc);
void cblas_dsyrk_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                 const double alpha, const double *A, const int64_t lda,
                 const double beta, double *C, const int64_t ldc);
void cblas_dsyr2k_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                  const double alpha, const double *A, const int64_t lda,
                  const double *B, const int64_t ldb, const double beta,
                  double *C, const int64_t ldc);
void cblas_dtrmm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int64_t M, const int64_t N,
                 const double alpha, const double *A, const int64_t lda,
                 double *B, const int64_t ldb);
void cblas_dtrsm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int64_t M, const int64_t N,
                 const double alpha, const double *A, const int64_t lda,
                 double *B, const int64_t ldb);

void cblas_cgemm_64(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int64_t M, const int64_t N,
                 const int64_t K, const void *alpha, const void *A,
                 const int64_t lda, const void *B, const int64_t ldb,
                 const void *beta, void *C, const int64_t ldc);
void cblas_csymm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 const void *B, const int64_t ldb, const void *beta,
                 void *C, const int64_t ldc);
void cblas_csyrk_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                 const void *alpha, const void *A, const int64_t lda,
                 const void *beta, void *C, const int64_t ldc);
void cblas_csyr2k_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                  const void *alpha, const void *A, const int64_t lda,
                  const void *B, const int64_t ldb, const void *beta,
                  void *C, const int64_t ldc);
void cblas_ctrmm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 void *B, const int64_t ldb);
void cblas_ctrsm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 void *B, const int64_t ldb);

void cblas_zgemm_64(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int64_t M, const int64_t N,
                 const int64_t K, const void *alpha, const void *A,
                 const int64_t lda, const void *B, const int64_t ldb,
                 const void *beta, void *C, const int64_t ldc);
void cblas_zsymm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 const void *B, const int64_t ldb, const void *beta,
                 void *C, const int64_t ldc);
void cblas_zsyrk_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                 const void *alpha, const void *A, const int64_t lda,
                 const void *beta, void *C, const int64_t ldc);
void cblas_zsyr2k_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                  const void *alpha, const void *A, const int64_t lda,
                  const void *B, const int64_t ldb, const void *beta,
                  void *C, const int64_t ldc);
void cblas_ztrmm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 void *B, const int64_t ldb);
void cblas_ztrsm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 void *B, const int64_t ldb);


/*
 * Routines with prefixes C and Z only
 */
void cblas_chemm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 const void *B, const int64_t ldb, const void *beta,
                 void *C, const int64_t ldc);
void cblas_cherk_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                 const float alpha, const void *A, const int64_t lda,
                 const float beta, void *C, const int64_t ldc);
void cblas_cher2k_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                  const void *alpha, const void *A, const int64_t lda,
                  const void *B, const int64_t ldb, const float beta,
                  void *C, const int64_t ldc);

void cblas_zhemm_64(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int64_t M, const int64_t N,
                 const void *alpha, const void *A, const int64_t lda,
                 const void *B, const int64_t ldb, const void *beta,
                 void *C, const int64_t ldc);
void cblas_zherk_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                 const double alpha, const void *A, const int64_t lda,
                 const double beta, void *C, const int64_t ldc);
void cblas_zher2k_64(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int64_t N, const int64_t K,
                  const void *alpha, const void *A, const int64_t lda,
                  const void *B, const int64_t ldb, const double beta,
                  void *C, const int64_t ldc);

void
#ifdef HAS_ATTRIBUTE_WEAK_SUPPORT
__attribute__((weak))
#endif
cblas_xerbla_64(int64_t p, const char *rout, const char *form, ...);

#ifdef __cplusplus
}
#endif
#endif
