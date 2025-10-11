/*****************************************************************************
  Copyright (c) 2014, Intel Corp.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************
* Contents: Native C interface to LAPACK
* Author: Intel Corporation
*****************************************************************************/

#ifndef _LAPACKE_64_H_
#define _LAPACKE_64_H_

#include "lapacke.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* C-LAPACK function prototypes */

int64_t LAPACKE_sbdsdc_64( int matrix_layout, char uplo, char compq,
                           int64_t n, float* d, float* e, float* u,
                           int64_t ldu, float* vt, int64_t ldvt, float* q,
                           int64_t* iq );
int64_t LAPACKE_dbdsdc_64( int matrix_layout, char uplo, char compq,
                           int64_t n, double* d, double* e, double* u,
                           int64_t ldu, double* vt, int64_t ldvt,
                           double* q, int64_t* iq );

int64_t LAPACKE_sbdsqr_64( int matrix_layout, char uplo, int64_t n,
                           int64_t ncvt, int64_t nru, int64_t ncc,
                           float* d, float* e, float* vt, int64_t ldvt,
                           float* u, int64_t ldu, float* c, int64_t ldc );
int64_t LAPACKE_dbdsqr_64( int matrix_layout, char uplo, int64_t n,
                           int64_t ncvt, int64_t nru, int64_t ncc,
                           double* d, double* e, double* vt, int64_t ldvt,
                           double* u, int64_t ldu, double* c,
                           int64_t ldc );
int64_t LAPACKE_cbdsqr_64( int matrix_layout, char uplo, int64_t n,
                           int64_t ncvt, int64_t nru, int64_t ncc,
                           float* d, float* e, lapack_complex_float* vt,
                           int64_t ldvt, lapack_complex_float* u,
                           int64_t ldu, lapack_complex_float* c,
                           int64_t ldc );
int64_t LAPACKE_zbdsqr_64( int matrix_layout, char uplo, int64_t n,
                           int64_t ncvt, int64_t nru, int64_t ncc,
                           double* d, double* e, lapack_complex_double* vt,
                           int64_t ldvt, lapack_complex_double* u,
                           int64_t ldu, lapack_complex_double* c,
                           int64_t ldc );
int64_t LAPACKE_sbdsvdx_64( int matrix_layout, char uplo, char jobz, char range,
                           int64_t n, float* d, float* e,
                           float vl, float vu,
                           int64_t il, int64_t iu, int64_t* ns,
                           float* s, float* z, int64_t ldz,
                           int64_t* superb );
int64_t LAPACKE_dbdsvdx_64( int matrix_layout, char uplo, char jobz, char range,
                           int64_t n, double* d, double* e,
                           double vl, double vu,
                           int64_t il, int64_t iu, int64_t* ns,
                           double* s, double* z, int64_t ldz,
                           int64_t* superb );
int64_t LAPACKE_sdisna_64( char job, int64_t m, int64_t n, const float* d,
                           float* sep );
int64_t LAPACKE_ddisna_64( char job, int64_t m, int64_t n,
                           const double* d, double* sep );

int64_t LAPACKE_sgbbrd_64( int matrix_layout, char vect, int64_t m,
                           int64_t n, int64_t ncc, int64_t kl,
                           int64_t ku, float* ab, int64_t ldab, float* d,
                           float* e, float* q, int64_t ldq, float* pt,
                           int64_t ldpt, float* c, int64_t ldc );
int64_t LAPACKE_dgbbrd_64( int matrix_layout, char vect, int64_t m,
                           int64_t n, int64_t ncc, int64_t kl,
                           int64_t ku, double* ab, int64_t ldab,
                           double* d, double* e, double* q, int64_t ldq,
                           double* pt, int64_t ldpt, double* c,
                           int64_t ldc );
int64_t LAPACKE_cgbbrd_64( int matrix_layout, char vect, int64_t m,
                           int64_t n, int64_t ncc, int64_t kl,
                           int64_t ku, lapack_complex_float* ab,
                           int64_t ldab, float* d, float* e,
                           lapack_complex_float* q, int64_t ldq,
                           lapack_complex_float* pt, int64_t ldpt,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zgbbrd_64( int matrix_layout, char vect, int64_t m,
                           int64_t n, int64_t ncc, int64_t kl,
                           int64_t ku, lapack_complex_double* ab,
                           int64_t ldab, double* d, double* e,
                           lapack_complex_double* q, int64_t ldq,
                           lapack_complex_double* pt, int64_t ldpt,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_sgbcon_64( int matrix_layout, char norm, int64_t n,
                           int64_t kl, int64_t ku, const float* ab,
                           int64_t ldab, const int64_t* ipiv, float anorm,
                           float* rcond );
int64_t LAPACKE_dgbcon_64( int matrix_layout, char norm, int64_t n,
                           int64_t kl, int64_t ku, const double* ab,
                           int64_t ldab, const int64_t* ipiv,
                           double anorm, double* rcond );
int64_t LAPACKE_cgbcon_64( int matrix_layout, char norm, int64_t n,
                           int64_t kl, int64_t ku,
                           const lapack_complex_float* ab, int64_t ldab,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_zgbcon_64( int matrix_layout, char norm, int64_t n,
                           int64_t kl, int64_t ku,
                           const lapack_complex_double* ab, int64_t ldab,
                           const int64_t* ipiv, double anorm,
                           double* rcond );

int64_t LAPACKE_sgbequ_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku, const float* ab,
                           int64_t ldab, float* r, float* c, float* rowcnd,
                           float* colcnd, float* amax );
int64_t LAPACKE_dgbequ_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku, const double* ab,
                           int64_t ldab, double* r, double* c,
                           double* rowcnd, double* colcnd, double* amax );
int64_t LAPACKE_cgbequ_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku,
                           const lapack_complex_float* ab, int64_t ldab,
                           float* r, float* c, float* rowcnd, float* colcnd,
                           float* amax );
int64_t LAPACKE_zgbequ_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku,
                           const lapack_complex_double* ab, int64_t ldab,
                           double* r, double* c, double* rowcnd, double* colcnd,
                           double* amax );

int64_t LAPACKE_sgbequb_64( int matrix_layout, int64_t m, int64_t n,
                            int64_t kl, int64_t ku, const float* ab,
                            int64_t ldab, float* r, float* c, float* rowcnd,
                            float* colcnd, float* amax );
int64_t LAPACKE_dgbequb_64( int matrix_layout, int64_t m, int64_t n,
                            int64_t kl, int64_t ku, const double* ab,
                            int64_t ldab, double* r, double* c,
                            double* rowcnd, double* colcnd, double* amax );
int64_t LAPACKE_cgbequb_64( int matrix_layout, int64_t m, int64_t n,
                            int64_t kl, int64_t ku,
                            const lapack_complex_float* ab, int64_t ldab,
                            float* r, float* c, float* rowcnd, float* colcnd,
                            float* amax );
int64_t LAPACKE_zgbequb_64( int matrix_layout, int64_t m, int64_t n,
                            int64_t kl, int64_t ku,
                            const lapack_complex_double* ab, int64_t ldab,
                            double* r, double* c, double* rowcnd,
                            double* colcnd, double* amax );

int64_t LAPACKE_sgbrfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t kl, int64_t ku, int64_t nrhs,
                           const float* ab, int64_t ldab, const float* afb,
                           int64_t ldafb, const int64_t* ipiv,
                           const float* b, int64_t ldb, float* x,
                           int64_t ldx, float* ferr, float* berr );
int64_t LAPACKE_dgbrfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t kl, int64_t ku, int64_t nrhs,
                           const double* ab, int64_t ldab, const double* afb,
                           int64_t ldafb, const int64_t* ipiv,
                           const double* b, int64_t ldb, double* x,
                           int64_t ldx, double* ferr, double* berr );
int64_t LAPACKE_cgbrfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t kl, int64_t ku, int64_t nrhs,
                           const lapack_complex_float* ab, int64_t ldab,
                           const lapack_complex_float* afb, int64_t ldafb,
                           const int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zgbrfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t kl, int64_t ku, int64_t nrhs,
                           const lapack_complex_double* ab, int64_t ldab,
                           const lapack_complex_double* afb, int64_t ldafb,
                           const int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_sgbrfsx_64( int matrix_layout, char trans, char equed,
                            int64_t n, int64_t kl, int64_t ku,
                            int64_t nrhs, const float* ab, int64_t ldab,
                            const float* afb, int64_t ldafb,
                            const int64_t* ipiv, const float* r,
                            const float* c, const float* b, int64_t ldb,
                            float* x, int64_t ldx, float* rcond, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_dgbrfsx_64( int matrix_layout, char trans, char equed,
                            int64_t n, int64_t kl, int64_t ku,
                            int64_t nrhs, const double* ab, int64_t ldab,
                            const double* afb, int64_t ldafb,
                            const int64_t* ipiv, const double* r,
                            const double* c, const double* b, int64_t ldb,
                            double* x, int64_t ldx, double* rcond,
                            double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );
int64_t LAPACKE_cgbrfsx_64( int matrix_layout, char trans, char equed,
                            int64_t n, int64_t kl, int64_t ku,
                            int64_t nrhs, const lapack_complex_float* ab,
                            int64_t ldab, const lapack_complex_float* afb,
                            int64_t ldafb, const int64_t* ipiv,
                            const float* r, const float* c,
                            const lapack_complex_float* b, int64_t ldb,
                            lapack_complex_float* x, int64_t ldx,
                            float* rcond, float* berr, int64_t n_err_bnds,
                            float* err_bnds_norm, float* err_bnds_comp,
                            int64_t nparams, float* params );
int64_t LAPACKE_zgbrfsx_64( int matrix_layout, char trans, char equed,
                            int64_t n, int64_t kl, int64_t ku,
                            int64_t nrhs, const lapack_complex_double* ab,
                            int64_t ldab, const lapack_complex_double* afb,
                            int64_t ldafb, const int64_t* ipiv,
                            const double* r, const double* c,
                            const lapack_complex_double* b, int64_t ldb,
                            lapack_complex_double* x, int64_t ldx,
                            double* rcond, double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );

int64_t LAPACKE_sgbsv_64( int matrix_layout, int64_t n, int64_t kl,
                          int64_t ku, int64_t nrhs, float* ab,
                          int64_t ldab, int64_t* ipiv, float* b,
                          int64_t ldb );
int64_t LAPACKE_dgbsv_64( int matrix_layout, int64_t n, int64_t kl,
                          int64_t ku, int64_t nrhs, double* ab,
                          int64_t ldab, int64_t* ipiv, double* b,
                          int64_t ldb );
int64_t LAPACKE_cgbsv_64( int matrix_layout, int64_t n, int64_t kl,
                          int64_t ku, int64_t nrhs,
                          lapack_complex_float* ab, int64_t ldab,
                          int64_t* ipiv, lapack_complex_float* b,
                          int64_t ldb );
int64_t LAPACKE_zgbsv_64( int matrix_layout, int64_t n, int64_t kl,
                          int64_t ku, int64_t nrhs,
                          lapack_complex_double* ab, int64_t ldab,
                          int64_t* ipiv, lapack_complex_double* b,
                          int64_t ldb );

int64_t LAPACKE_sgbsvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t kl, int64_t ku,
                           int64_t nrhs, float* ab, int64_t ldab,
                           float* afb, int64_t ldafb, int64_t* ipiv,
                           char* equed, float* r, float* c, float* b,
                           int64_t ldb, float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr,
                           float* rpivot );
int64_t LAPACKE_dgbsvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t kl, int64_t ku,
                           int64_t nrhs, double* ab, int64_t ldab,
                           double* afb, int64_t ldafb, int64_t* ipiv,
                           char* equed, double* r, double* c, double* b,
                           int64_t ldb, double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr,
                           double* rpivot );
int64_t LAPACKE_cgbsvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t kl, int64_t ku,
                           int64_t nrhs, lapack_complex_float* ab,
                           int64_t ldab, lapack_complex_float* afb,
                           int64_t ldafb, int64_t* ipiv, char* equed,
                           float* r, float* c, lapack_complex_float* b,
                           int64_t ldb, lapack_complex_float* x,
                           int64_t ldx, float* rcond, float* ferr,
                           float* berr, float* rpivot );
int64_t LAPACKE_zgbsvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t kl, int64_t ku,
                           int64_t nrhs, lapack_complex_double* ab,
                           int64_t ldab, lapack_complex_double* afb,
                           int64_t ldafb, int64_t* ipiv, char* equed,
                           double* r, double* c, lapack_complex_double* b,
                           int64_t ldb, lapack_complex_double* x,
                           int64_t ldx, double* rcond, double* ferr,
                           double* berr, double* rpivot );

int64_t LAPACKE_sgbsvxx_64( int matrix_layout, char fact, char trans,
                            int64_t n, int64_t kl, int64_t ku,
                            int64_t nrhs, float* ab, int64_t ldab,
                            float* afb, int64_t ldafb, int64_t* ipiv,
                            char* equed, float* r, float* c, float* b,
                            int64_t ldb, float* x, int64_t ldx,
                            float* rcond, float* rpvgrw, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_dgbsvxx_64( int matrix_layout, char fact, char trans,
                            int64_t n, int64_t kl, int64_t ku,
                            int64_t nrhs, double* ab, int64_t ldab,
                            double* afb, int64_t ldafb, int64_t* ipiv,
                            char* equed, double* r, double* c, double* b,
                            int64_t ldb, double* x, int64_t ldx,
                            double* rcond, double* rpvgrw, double* berr,
                            int64_t n_err_bnds, double* err_bnds_norm,
                            double* err_bnds_comp, int64_t nparams,
                            double* params );
int64_t LAPACKE_cgbsvxx_64( int matrix_layout, char fact, char trans,
                            int64_t n, int64_t kl, int64_t ku,
                            int64_t nrhs, lapack_complex_float* ab,
                            int64_t ldab, lapack_complex_float* afb,
                            int64_t ldafb, int64_t* ipiv, char* equed,
                            float* r, float* c, lapack_complex_float* b,
                            int64_t ldb, lapack_complex_float* x,
                            int64_t ldx, float* rcond, float* rpvgrw,
                            float* berr, int64_t n_err_bnds,
                            float* err_bnds_norm, float* err_bnds_comp,
                            int64_t nparams, float* params );
int64_t LAPACKE_zgbsvxx_64( int matrix_layout, char fact, char trans,
                            int64_t n, int64_t kl, int64_t ku,
                            int64_t nrhs, lapack_complex_double* ab,
                            int64_t ldab, lapack_complex_double* afb,
                            int64_t ldafb, int64_t* ipiv, char* equed,
                            double* r, double* c, lapack_complex_double* b,
                            int64_t ldb, lapack_complex_double* x,
                            int64_t ldx, double* rcond, double* rpvgrw,
                            double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );

int64_t LAPACKE_sgbtrf_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku, float* ab,
                           int64_t ldab, int64_t* ipiv );
int64_t LAPACKE_dgbtrf_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku, double* ab,
                           int64_t ldab, int64_t* ipiv );
int64_t LAPACKE_cgbtrf_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku,
                           lapack_complex_float* ab, int64_t ldab,
                           int64_t* ipiv );
int64_t LAPACKE_zgbtrf_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku,
                           lapack_complex_double* ab, int64_t ldab,
                           int64_t* ipiv );

int64_t LAPACKE_sgbtrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t kl, int64_t ku, int64_t nrhs,
                           const float* ab, int64_t ldab,
                           const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dgbtrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t kl, int64_t ku, int64_t nrhs,
                           const double* ab, int64_t ldab,
                           const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_cgbtrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t kl, int64_t ku, int64_t nrhs,
                           const lapack_complex_float* ab, int64_t ldab,
                           const int64_t* ipiv, lapack_complex_float* b,
                           int64_t ldb );
int64_t LAPACKE_zgbtrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t kl, int64_t ku, int64_t nrhs,
                           const lapack_complex_double* ab, int64_t ldab,
                           const int64_t* ipiv, lapack_complex_double* b,
                           int64_t ldb );

int64_t LAPACKE_sgebak_64( int matrix_layout, char job, char side, int64_t n,
                           int64_t ilo, int64_t ihi, const float* scale,
                           int64_t m, float* v, int64_t ldv );
int64_t LAPACKE_dgebak_64( int matrix_layout, char job, char side, int64_t n,
                           int64_t ilo, int64_t ihi, const double* scale,
                           int64_t m, double* v, int64_t ldv );
int64_t LAPACKE_cgebak_64( int matrix_layout, char job, char side, int64_t n,
                           int64_t ilo, int64_t ihi, const float* scale,
                           int64_t m, lapack_complex_float* v,
                           int64_t ldv );
int64_t LAPACKE_zgebak_64( int matrix_layout, char job, char side, int64_t n,
                           int64_t ilo, int64_t ihi, const double* scale,
                           int64_t m, lapack_complex_double* v,
                           int64_t ldv );

int64_t LAPACKE_sgebal_64( int matrix_layout, char job, int64_t n, float* a,
                           int64_t lda, int64_t* ilo, int64_t* ihi,
                           float* scale );
int64_t LAPACKE_dgebal_64( int matrix_layout, char job, int64_t n, double* a,
                           int64_t lda, int64_t* ilo, int64_t* ihi,
                           double* scale );
int64_t LAPACKE_cgebal_64( int matrix_layout, char job, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ilo, int64_t* ihi, float* scale );
int64_t LAPACKE_zgebal_64( int matrix_layout, char job, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ilo, int64_t* ihi, double* scale );

int64_t LAPACKE_sgebrd_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, float* d, float* e,
                           float* tauq, float* taup );
int64_t LAPACKE_dgebrd_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, double* d, double* e,
                           double* tauq, double* taup );
int64_t LAPACKE_cgebrd_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda, float* d,
                           float* e, lapack_complex_float* tauq,
                           lapack_complex_float* taup );
int64_t LAPACKE_zgebrd_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda, double* d,
                           double* e, lapack_complex_double* tauq,
                           lapack_complex_double* taup );

int64_t LAPACKE_sgecon_64( int matrix_layout, char norm, int64_t n,
                           const float* a, int64_t lda, float anorm,
                           float* rcond );
int64_t LAPACKE_dgecon_64( int matrix_layout, char norm, int64_t n,
                           const double* a, int64_t lda, double anorm,
                           double* rcond );
int64_t LAPACKE_cgecon_64( int matrix_layout, char norm, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           float anorm, float* rcond );
int64_t LAPACKE_zgecon_64( int matrix_layout, char norm, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           double anorm, double* rcond );

int64_t LAPACKE_sgeequ_64( int matrix_layout, int64_t m, int64_t n,
                           const float* a, int64_t lda, float* r, float* c,
                           float* rowcnd, float* colcnd, float* amax );
int64_t LAPACKE_dgeequ_64( int matrix_layout, int64_t m, int64_t n,
                           const double* a, int64_t lda, double* r,
                           double* c, double* rowcnd, double* colcnd,
                           double* amax );
int64_t LAPACKE_cgeequ_64( int matrix_layout, int64_t m, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           float* r, float* c, float* rowcnd, float* colcnd,
                           float* amax );
int64_t LAPACKE_zgeequ_64( int matrix_layout, int64_t m, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           double* r, double* c, double* rowcnd, double* colcnd,
                           double* amax );

int64_t LAPACKE_sgeequb_64( int matrix_layout, int64_t m, int64_t n,
                            const float* a, int64_t lda, float* r, float* c,
                            float* rowcnd, float* colcnd, float* amax );
int64_t LAPACKE_dgeequb_64( int matrix_layout, int64_t m, int64_t n,
                            const double* a, int64_t lda, double* r,
                            double* c, double* rowcnd, double* colcnd,
                            double* amax );
int64_t LAPACKE_cgeequb_64( int matrix_layout, int64_t m, int64_t n,
                            const lapack_complex_float* a, int64_t lda,
                            float* r, float* c, float* rowcnd, float* colcnd,
                            float* amax );
int64_t LAPACKE_zgeequb_64( int matrix_layout, int64_t m, int64_t n,
                            const lapack_complex_double* a, int64_t lda,
                            double* r, double* c, double* rowcnd,
                            double* colcnd, double* amax );

int64_t LAPACKE_sgees_64( int matrix_layout, char jobvs, char sort,
                          LAPACK_S_SELECT2 select, int64_t n, float* a,
                          int64_t lda, int64_t* sdim, float* wr,
                          float* wi, float* vs, int64_t ldvs );
int64_t LAPACKE_dgees_64( int matrix_layout, char jobvs, char sort,
                          LAPACK_D_SELECT2 select, int64_t n, double* a,
                          int64_t lda, int64_t* sdim, double* wr,
                          double* wi, double* vs, int64_t ldvs );
int64_t LAPACKE_cgees_64( int matrix_layout, char jobvs, char sort,
                          LAPACK_C_SELECT1 select, int64_t n,
                          lapack_complex_float* a, int64_t lda,
                          int64_t* sdim, lapack_complex_float* w,
                          lapack_complex_float* vs, int64_t ldvs );
int64_t LAPACKE_zgees_64( int matrix_layout, char jobvs, char sort,
                          LAPACK_Z_SELECT1 select, int64_t n,
                          lapack_complex_double* a, int64_t lda,
                          int64_t* sdim, lapack_complex_double* w,
                          lapack_complex_double* vs, int64_t ldvs );

int64_t LAPACKE_sgeesx_64( int matrix_layout, char jobvs, char sort,
                           LAPACK_S_SELECT2 select, char sense, int64_t n,
                           float* a, int64_t lda, int64_t* sdim,
                           float* wr, float* wi, float* vs, int64_t ldvs,
                           float* rconde, float* rcondv );
int64_t LAPACKE_dgeesx_64( int matrix_layout, char jobvs, char sort,
                           LAPACK_D_SELECT2 select, char sense, int64_t n,
                           double* a, int64_t lda, int64_t* sdim,
                           double* wr, double* wi, double* vs, int64_t ldvs,
                           double* rconde, double* rcondv );
int64_t LAPACKE_cgeesx_64( int matrix_layout, char jobvs, char sort,
                           LAPACK_C_SELECT1 select, char sense, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* sdim, lapack_complex_float* w,
                           lapack_complex_float* vs, int64_t ldvs,
                           float* rconde, float* rcondv );
int64_t LAPACKE_zgeesx_64( int matrix_layout, char jobvs, char sort,
                           LAPACK_Z_SELECT1 select, char sense, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* sdim, lapack_complex_double* w,
                           lapack_complex_double* vs, int64_t ldvs,
                           double* rconde, double* rcondv );

int64_t LAPACKE_sgeev_64( int matrix_layout, char jobvl, char jobvr,
                          int64_t n, float* a, int64_t lda, float* wr,
                          float* wi, float* vl, int64_t ldvl, float* vr,
                          int64_t ldvr );
int64_t LAPACKE_dgeev_64( int matrix_layout, char jobvl, char jobvr,
                          int64_t n, double* a, int64_t lda, double* wr,
                          double* wi, double* vl, int64_t ldvl, double* vr,
                          int64_t ldvr );
int64_t LAPACKE_cgeev_64( int matrix_layout, char jobvl, char jobvr,
                          int64_t n, lapack_complex_float* a, int64_t lda,
                          lapack_complex_float* w, lapack_complex_float* vl,
                          int64_t ldvl, lapack_complex_float* vr,
                          int64_t ldvr );
int64_t LAPACKE_zgeev_64( int matrix_layout, char jobvl, char jobvr,
                          int64_t n, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* w,
                          lapack_complex_double* vl, int64_t ldvl,
                          lapack_complex_double* vr, int64_t ldvr );

int64_t LAPACKE_sgeevx_64( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, int64_t n, float* a,
                           int64_t lda, float* wr, float* wi, float* vl,
                           int64_t ldvl, float* vr, int64_t ldvr,
                           int64_t* ilo, int64_t* ihi, float* scale,
                           float* abnrm, float* rconde, float* rcondv );
int64_t LAPACKE_dgeevx_64( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, int64_t n, double* a,
                           int64_t lda, double* wr, double* wi, double* vl,
                           int64_t ldvl, double* vr, int64_t ldvr,
                           int64_t* ilo, int64_t* ihi, double* scale,
                           double* abnrm, double* rconde, double* rcondv );
int64_t LAPACKE_cgeevx_64( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* w, lapack_complex_float* vl,
                           int64_t ldvl, lapack_complex_float* vr,
                           int64_t ldvr, int64_t* ilo, int64_t* ihi,
                           float* scale, float* abnrm, float* rconde,
                           float* rcondv );
int64_t LAPACKE_zgeevx_64( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* w, lapack_complex_double* vl,
                           int64_t ldvl, lapack_complex_double* vr,
                           int64_t ldvr, int64_t* ilo, int64_t* ihi,
                           double* scale, double* abnrm, double* rconde,
                           double* rcondv );

int64_t LAPACKE_sgehrd_64( int matrix_layout, int64_t n, int64_t ilo,
                           int64_t ihi, float* a, int64_t lda,
                           float* tau );
int64_t LAPACKE_dgehrd_64( int matrix_layout, int64_t n, int64_t ilo,
                           int64_t ihi, double* a, int64_t lda,
                           double* tau );
int64_t LAPACKE_cgehrd_64( int matrix_layout, int64_t n, int64_t ilo,
                           int64_t ihi, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* tau );
int64_t LAPACKE_zgehrd_64( int matrix_layout, int64_t n, int64_t ilo,
                           int64_t ihi, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* tau );

int64_t LAPACKE_sgejsv_64( int matrix_layout, char joba, char jobu, char jobv,
                           char jobr, char jobt, char jobp, int64_t m,
                           int64_t n, float* a, int64_t lda, float* sva,
                           float* u, int64_t ldu, float* v, int64_t ldv,
                           float* stat, int64_t* istat );
int64_t LAPACKE_dgejsv_64( int matrix_layout, char joba, char jobu, char jobv,
                           char jobr, char jobt, char jobp, int64_t m,
                           int64_t n, double* a, int64_t lda, double* sva,
                           double* u, int64_t ldu, double* v, int64_t ldv,
                           double* stat, int64_t* istat );
int64_t LAPACKE_cgejsv_64( int matrix_layout, char joba, char jobu, char jobv,
                           char jobr, char jobt, char jobp, int64_t m,
                           int64_t n, lapack_complex_float* a, int64_t lda, float* sva,
                           lapack_complex_float* u, int64_t ldu, lapack_complex_float* v, int64_t ldv,
                           float* stat, int64_t* istat );
int64_t LAPACKE_zgejsv_64( int matrix_layout, char joba, char jobu, char jobv,
                           char jobr, char jobt, char jobp, int64_t m,
                           int64_t n, lapack_complex_double* a, int64_t lda, double* sva,
                           lapack_complex_double* u, int64_t ldu, lapack_complex_double* v, int64_t ldv,
                           double* stat, int64_t* istat );

int64_t LAPACKE_sgelq2_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, float* tau );
int64_t LAPACKE_dgelq2_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, double* tau );
int64_t LAPACKE_cgelq2_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* tau );
int64_t LAPACKE_zgelq2_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* tau );

int64_t LAPACKE_sgelqf_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, float* tau );
int64_t LAPACKE_dgelqf_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, double* tau );
int64_t LAPACKE_cgelqf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* tau );
int64_t LAPACKE_zgelqf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* tau );

int64_t LAPACKE_sgels_64( int matrix_layout, char trans, int64_t m,
                          int64_t n, int64_t nrhs, float* a,
                          int64_t lda, float* b, int64_t ldb );
int64_t LAPACKE_dgels_64( int matrix_layout, char trans, int64_t m,
                          int64_t n, int64_t nrhs, double* a,
                          int64_t lda, double* b, int64_t ldb );
int64_t LAPACKE_cgels_64( int matrix_layout, char trans, int64_t m,
                          int64_t n, int64_t nrhs,
                          lapack_complex_float* a, int64_t lda,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zgels_64( int matrix_layout, char trans, int64_t m,
                          int64_t n, int64_t nrhs,
                          lapack_complex_double* a, int64_t lda,
                          lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sgelsd_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, float* a, int64_t lda, float* b,
                           int64_t ldb, float* s, float rcond,
                           int64_t* rank );
int64_t LAPACKE_dgelsd_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, double* a, int64_t lda,
                           double* b, int64_t ldb, double* s, double rcond,
                           int64_t* rank );
int64_t LAPACKE_cgelsd_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb, float* s, float rcond,
                           int64_t* rank );
int64_t LAPACKE_zgelsd_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb, double* s, double rcond,
                           int64_t* rank );

int64_t LAPACKE_sgelss_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, float* a, int64_t lda, float* b,
                           int64_t ldb, float* s, float rcond,
                           int64_t* rank );
int64_t LAPACKE_dgelss_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, double* a, int64_t lda,
                           double* b, int64_t ldb, double* s, double rcond,
                           int64_t* rank );
int64_t LAPACKE_cgelss_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb, float* s, float rcond,
                           int64_t* rank );
int64_t LAPACKE_zgelss_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb, double* s, double rcond,
                           int64_t* rank );

int64_t LAPACKE_sgelsy_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, float* a, int64_t lda, float* b,
                           int64_t ldb, int64_t* jpvt, float rcond,
                           int64_t* rank );
int64_t LAPACKE_dgelsy_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, double* a, int64_t lda,
                           double* b, int64_t ldb, int64_t* jpvt,
                           double rcond, int64_t* rank );
int64_t LAPACKE_cgelsy_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb, int64_t* jpvt, float rcond,
                           int64_t* rank );
int64_t LAPACKE_zgelsy_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nrhs, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb, int64_t* jpvt, double rcond,
                           int64_t* rank );

int64_t LAPACKE_sgeqlf_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, float* tau );
int64_t LAPACKE_dgeqlf_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, double* tau );
int64_t LAPACKE_cgeqlf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* tau );
int64_t LAPACKE_zgeqlf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* tau );

int64_t LAPACKE_sgeqp3_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, int64_t* jpvt,
                           float* tau );
int64_t LAPACKE_dgeqp3_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, int64_t* jpvt,
                           double* tau );
int64_t LAPACKE_cgeqp3_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* jpvt, lapack_complex_float* tau );
int64_t LAPACKE_zgeqp3_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* jpvt, lapack_complex_double* tau );

int64_t LAPACKE_sgeqpf_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, int64_t* jpvt,
                           float* tau );
int64_t LAPACKE_dgeqpf_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, int64_t* jpvt,
                           double* tau );
int64_t LAPACKE_cgeqpf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* jpvt, lapack_complex_float* tau );
int64_t LAPACKE_zgeqpf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* jpvt, lapack_complex_double* tau );

int64_t LAPACKE_sgeqr2_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, float* tau );
int64_t LAPACKE_dgeqr2_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, double* tau );
int64_t LAPACKE_cgeqr2_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* tau );
int64_t LAPACKE_zgeqr2_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* tau );

int64_t LAPACKE_sgeqrf_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, float* tau );
int64_t LAPACKE_dgeqrf_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, double* tau );
int64_t LAPACKE_cgeqrf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* tau );
int64_t LAPACKE_zgeqrf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* tau );

int64_t LAPACKE_sgeqrfp_64( int matrix_layout, int64_t m, int64_t n,
                            float* a, int64_t lda, float* tau );
int64_t LAPACKE_dgeqrfp_64( int matrix_layout, int64_t m, int64_t n,
                            double* a, int64_t lda, double* tau );
int64_t LAPACKE_cgeqrfp_64( int matrix_layout, int64_t m, int64_t n,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* tau );
int64_t LAPACKE_zgeqrfp_64( int matrix_layout, int64_t m, int64_t n,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* tau );

int64_t LAPACKE_sgerfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const float* af, int64_t ldaf,
                           const int64_t* ipiv, const float* b,
                           int64_t ldb, float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_dgerfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           const double* af, int64_t ldaf,
                           const int64_t* ipiv, const double* b,
                           int64_t ldb, double* x, int64_t ldx,
                           double* ferr, double* berr );
int64_t LAPACKE_cgerfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* af,
                           int64_t ldaf, const int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zgerfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* af,
                           int64_t ldaf, const int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_sgerfsx_64( int matrix_layout, char trans, char equed,
                            int64_t n, int64_t nrhs, const float* a,
                            int64_t lda, const float* af, int64_t ldaf,
                            const int64_t* ipiv, const float* r,
                            const float* c, const float* b, int64_t ldb,
                            float* x, int64_t ldx, float* rcond, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_dgerfsx_64( int matrix_layout, char trans, char equed,
                            int64_t n, int64_t nrhs, const double* a,
                            int64_t lda, const double* af, int64_t ldaf,
                            const int64_t* ipiv, const double* r,
                            const double* c, const double* b, int64_t ldb,
                            double* x, int64_t ldx, double* rcond,
                            double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );
int64_t LAPACKE_cgerfsx_64( int matrix_layout, char trans, char equed,
                            int64_t n, int64_t nrhs,
                            const lapack_complex_float* a, int64_t lda,
                            const lapack_complex_float* af, int64_t ldaf,
                            const int64_t* ipiv, const float* r,
                            const float* c, const lapack_complex_float* b,
                            int64_t ldb, lapack_complex_float* x,
                            int64_t ldx, float* rcond, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_zgerfsx_64( int matrix_layout, char trans, char equed,
                            int64_t n, int64_t nrhs,
                            const lapack_complex_double* a, int64_t lda,
                            const lapack_complex_double* af, int64_t ldaf,
                            const int64_t* ipiv, const double* r,
                            const double* c, const lapack_complex_double* b,
                            int64_t ldb, lapack_complex_double* x,
                            int64_t ldx, double* rcond, double* berr,
                            int64_t n_err_bnds, double* err_bnds_norm,
                            double* err_bnds_comp, int64_t nparams,
                            double* params );

int64_t LAPACKE_sgerqf_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, float* tau );
int64_t LAPACKE_dgerqf_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, double* tau );
int64_t LAPACKE_cgerqf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* tau );
int64_t LAPACKE_zgerqf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* tau );

int64_t LAPACKE_sgesdd_64( int matrix_layout, char jobz, int64_t m,
                           int64_t n, float* a, int64_t lda, float* s,
                           float* u, int64_t ldu, float* vt,
                           int64_t ldvt );
int64_t LAPACKE_dgesdd_64( int matrix_layout, char jobz, int64_t m,
                           int64_t n, double* a, int64_t lda, double* s,
                           double* u, int64_t ldu, double* vt,
                           int64_t ldvt );
int64_t LAPACKE_cgesdd_64( int matrix_layout, char jobz, int64_t m,
                           int64_t n, lapack_complex_float* a,
                           int64_t lda, float* s, lapack_complex_float* u,
                           int64_t ldu, lapack_complex_float* vt,
                           int64_t ldvt );
int64_t LAPACKE_zgesdd_64( int matrix_layout, char jobz, int64_t m,
                           int64_t n, lapack_complex_double* a,
                           int64_t lda, double* s, lapack_complex_double* u,
                           int64_t ldu, lapack_complex_double* vt,
                           int64_t ldvt );

int64_t LAPACKE_sgesv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          float* a, int64_t lda, int64_t* ipiv, float* b,
                          int64_t ldb );
int64_t LAPACKE_dgesv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          double* a, int64_t lda, int64_t* ipiv,
                          double* b, int64_t ldb );
int64_t LAPACKE_cgesv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          lapack_complex_float* a, int64_t lda,
                          int64_t* ipiv, lapack_complex_float* b,
                          int64_t ldb );
int64_t LAPACKE_zgesv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          lapack_complex_double* a, int64_t lda,
                          int64_t* ipiv, lapack_complex_double* b,
                          int64_t ldb );
int64_t LAPACKE_dsgesv_64( int matrix_layout, int64_t n, int64_t nrhs,
                           double* a, int64_t lda, int64_t* ipiv,
                           double* b, int64_t ldb, double* x, int64_t ldx,
                           int64_t* iter );
int64_t LAPACKE_zcgesv_64( int matrix_layout, int64_t n, int64_t nrhs,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv, lapack_complex_double* b,
                           int64_t ldb, lapack_complex_double* x,
                           int64_t ldx, int64_t* iter );

int64_t LAPACKE_sgesvd_64( int matrix_layout, char jobu, char jobvt,
                           int64_t m, int64_t n, float* a, int64_t lda,
                           float* s, float* u, int64_t ldu, float* vt,
                           int64_t ldvt, float* superb );
int64_t LAPACKE_dgesvd_64( int matrix_layout, char jobu, char jobvt,
                           int64_t m, int64_t n, double* a,
                           int64_t lda, double* s, double* u, int64_t ldu,
                           double* vt, int64_t ldvt, double* superb );
int64_t LAPACKE_cgesvd_64( int matrix_layout, char jobu, char jobvt,
                           int64_t m, int64_t n, lapack_complex_float* a,
                           int64_t lda, float* s, lapack_complex_float* u,
                           int64_t ldu, lapack_complex_float* vt,
                           int64_t ldvt, float* superb );
int64_t LAPACKE_zgesvd_64( int matrix_layout, char jobu, char jobvt,
                           int64_t m, int64_t n, lapack_complex_double* a,
                           int64_t lda, double* s, lapack_complex_double* u,
                           int64_t ldu, lapack_complex_double* vt,
                           int64_t ldvt, double* superb );

int64_t LAPACKE_sgesvdx_64( int matrix_layout, char jobu, char jobvt, char range,
                           int64_t m, int64_t n, float* a,
                           int64_t lda, float vl, float vu,
                           int64_t il, int64_t iu, int64_t* ns,
                           float* s, float* u, int64_t ldu,
                           float* vt, int64_t ldvt,
                           int64_t* superb );
int64_t LAPACKE_dgesvdx_64( int matrix_layout, char jobu, char jobvt, char range,
                           int64_t m, int64_t n, double* a,
                           int64_t lda, double vl, double vu,
                           int64_t il, int64_t iu, int64_t* ns,
                           double* s, double* u, int64_t ldu,
                           double* vt, int64_t ldvt,
                           int64_t* superb );
int64_t LAPACKE_cgesvdx_64( int matrix_layout, char jobu, char jobvt, char range,
                           int64_t m, int64_t n, lapack_complex_float* a,
                           int64_t lda, float vl, float vu,
                           int64_t il, int64_t iu, int64_t* ns,
                           float* s, lapack_complex_float* u, int64_t ldu,
                           lapack_complex_float* vt, int64_t ldvt,
                           int64_t* superb );
int64_t LAPACKE_zgesvdx_64( int matrix_layout, char jobu, char jobvt, char range,
                           int64_t m, int64_t n, lapack_complex_double* a,
                           int64_t lda, double vl, double vu,
                           int64_t il, int64_t iu, int64_t* ns,
                           double* s, lapack_complex_double* u, int64_t ldu,
                           lapack_complex_double* vt, int64_t ldvt,
                           int64_t* superb );

int64_t LAPACKE_sgesvdq_64( int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,
                           int64_t m, int64_t n, float* a, int64_t lda,
                           float* s, float* u, int64_t ldu, float* v,
                           int64_t ldv, int64_t* numrank );
int64_t LAPACKE_dgesvdq_64( int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,
                           int64_t m, int64_t n, double* a,
                           int64_t lda, double* s, double* u, int64_t ldu,
                           double* v, int64_t ldv, int64_t* numrank);
int64_t LAPACKE_cgesvdq_64( int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,
                           int64_t m, int64_t n, lapack_complex_float* a,
                           int64_t lda, float* s, lapack_complex_float* u,
                           int64_t ldu, lapack_complex_float* v,
                           int64_t ldv, int64_t* numrank );
int64_t LAPACKE_zgesvdq_64( int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,
                           int64_t m, int64_t n, lapack_complex_double* a,
                           int64_t lda, double* s, lapack_complex_double* u,
                           int64_t ldu, lapack_complex_double* v,
                           int64_t ldv, int64_t* numrank );

int64_t LAPACKE_sgesvj_64( int matrix_layout, char joba, char jobu, char jobv,
                           int64_t m, int64_t n, float* a, int64_t lda,
                           float* sva, int64_t mv, float* v, int64_t ldv,
                           float* stat );
int64_t LAPACKE_dgesvj_64( int matrix_layout, char joba, char jobu, char jobv,
                           int64_t m, int64_t n, double* a,
                           int64_t lda, double* sva, int64_t mv,
                           double* v, int64_t ldv, double* stat );
int64_t LAPACKE_cgesvj_64( int matrix_layout, char joba, char jobu, char jobv,
                           int64_t m, int64_t n, lapack_complex_float* a,
                           int64_t lda, float* sva, int64_t mv,
                           lapack_complex_float* v, int64_t ldv, float* stat );
int64_t LAPACKE_zgesvj_64( int matrix_layout, char joba, char jobu, char jobv,
                           int64_t m, int64_t n, lapack_complex_double* a,
                           int64_t lda, double* sva, int64_t mv,
                           lapack_complex_double* v, int64_t ldv, double* stat );

int64_t LAPACKE_sgesvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t nrhs, float* a,
                           int64_t lda, float* af, int64_t ldaf,
                           int64_t* ipiv, char* equed, float* r, float* c,
                           float* b, int64_t ldb, float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr,
                           float* rpivot );
int64_t LAPACKE_dgesvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t nrhs, double* a,
                           int64_t lda, double* af, int64_t ldaf,
                           int64_t* ipiv, char* equed, double* r, double* c,
                           double* b, int64_t ldb, double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr,
                           double* rpivot );
int64_t LAPACKE_cgesvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t nrhs,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* af, int64_t ldaf,
                           int64_t* ipiv, char* equed, float* r, float* c,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr,
                           float* rpivot );
int64_t LAPACKE_zgesvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t nrhs,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* af, int64_t ldaf,
                           int64_t* ipiv, char* equed, double* r, double* c,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr,
                           double* rpivot );

int64_t LAPACKE_sgesvxx_64( int matrix_layout, char fact, char trans,
                            int64_t n, int64_t nrhs, float* a,
                            int64_t lda, float* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, float* r, float* c,
                            float* b, int64_t ldb, float* x, int64_t ldx,
                            float* rcond, float* rpvgrw, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_dgesvxx_64( int matrix_layout, char fact, char trans,
                            int64_t n, int64_t nrhs, double* a,
                            int64_t lda, double* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, double* r, double* c,
                            double* b, int64_t ldb, double* x,
                            int64_t ldx, double* rcond, double* rpvgrw,
                            double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );
int64_t LAPACKE_cgesvxx_64( int matrix_layout, char fact, char trans,
                            int64_t n, int64_t nrhs,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, float* r, float* c,
                            lapack_complex_float* b, int64_t ldb,
                            lapack_complex_float* x, int64_t ldx,
                            float* rcond, float* rpvgrw, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_zgesvxx_64( int matrix_layout, char fact, char trans,
                            int64_t n, int64_t nrhs,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, double* r, double* c,
                            lapack_complex_double* b, int64_t ldb,
                            lapack_complex_double* x, int64_t ldx,
                            double* rcond, double* rpvgrw, double* berr,
                            int64_t n_err_bnds, double* err_bnds_norm,
                            double* err_bnds_comp, int64_t nparams,
                            double* params );

int64_t LAPACKE_sgetf2_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dgetf2_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_cgetf2_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_zgetf2_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv );

int64_t LAPACKE_sgetrf_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dgetrf_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_cgetrf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_zgetrf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv );

int64_t LAPACKE_sgetrf2_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dgetrf2_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_cgetrf2_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_zgetrf2_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv );

int64_t LAPACKE_sgetri_64( int matrix_layout, int64_t n, float* a,
                           int64_t lda, const int64_t* ipiv );
int64_t LAPACKE_dgetri_64( int matrix_layout, int64_t n, double* a,
                           int64_t lda, const int64_t* ipiv );
int64_t LAPACKE_cgetri_64( int matrix_layout, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           const int64_t* ipiv );
int64_t LAPACKE_zgetri_64( int matrix_layout, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           const int64_t* ipiv );

int64_t LAPACKE_sgetrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dgetrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_cgetrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zgetrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sggbak_64( int matrix_layout, char job, char side, int64_t n,
                           int64_t ilo, int64_t ihi, const float* lscale,
                           const float* rscale, int64_t m, float* v,
                           int64_t ldv );
int64_t LAPACKE_dggbak_64( int matrix_layout, char job, char side, int64_t n,
                           int64_t ilo, int64_t ihi, const double* lscale,
                           const double* rscale, int64_t m, double* v,
                           int64_t ldv );
int64_t LAPACKE_cggbak_64( int matrix_layout, char job, char side, int64_t n,
                           int64_t ilo, int64_t ihi, const float* lscale,
                           const float* rscale, int64_t m,
                           lapack_complex_float* v, int64_t ldv );
int64_t LAPACKE_zggbak_64( int matrix_layout, char job, char side, int64_t n,
                           int64_t ilo, int64_t ihi, const double* lscale,
                           const double* rscale, int64_t m,
                           lapack_complex_double* v, int64_t ldv );

int64_t LAPACKE_sggbal_64( int matrix_layout, char job, int64_t n, float* a,
                           int64_t lda, float* b, int64_t ldb,
                           int64_t* ilo, int64_t* ihi, float* lscale,
                           float* rscale );
int64_t LAPACKE_dggbal_64( int matrix_layout, char job, int64_t n, double* a,
                           int64_t lda, double* b, int64_t ldb,
                           int64_t* ilo, int64_t* ihi, double* lscale,
                           double* rscale );
int64_t LAPACKE_cggbal_64( int matrix_layout, char job, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           int64_t* ilo, int64_t* ihi, float* lscale,
                           float* rscale );
int64_t LAPACKE_zggbal_64( int matrix_layout, char job, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           int64_t* ilo, int64_t* ihi, double* lscale,
                           double* rscale );

int64_t LAPACKE_sgges_64( int matrix_layout, char jobvsl, char jobvsr, char sort,
                          LAPACK_S_SELECT3 selctg, int64_t n, float* a,
                          int64_t lda, float* b, int64_t ldb,
                          int64_t* sdim, float* alphar, float* alphai,
                          float* beta, float* vsl, int64_t ldvsl, float* vsr,
                          int64_t ldvsr );
int64_t LAPACKE_dgges_64( int matrix_layout, char jobvsl, char jobvsr, char sort,
                          LAPACK_D_SELECT3 selctg, int64_t n, double* a,
                          int64_t lda, double* b, int64_t ldb,
                          int64_t* sdim, double* alphar, double* alphai,
                          double* beta, double* vsl, int64_t ldvsl,
                          double* vsr, int64_t ldvsr );
int64_t LAPACKE_cgges_64( int matrix_layout, char jobvsl, char jobvsr, char sort,
                          LAPACK_C_SELECT2 selctg, int64_t n,
                          lapack_complex_float* a, int64_t lda,
                          lapack_complex_float* b, int64_t ldb,
                          int64_t* sdim, lapack_complex_float* alpha,
                          lapack_complex_float* beta, lapack_complex_float* vsl,
                          int64_t ldvsl, lapack_complex_float* vsr,
                          int64_t ldvsr );
int64_t LAPACKE_zgges_64( int matrix_layout, char jobvsl, char jobvsr, char sort,
                          LAPACK_Z_SELECT2 selctg, int64_t n,
                          lapack_complex_double* a, int64_t lda,
                          lapack_complex_double* b, int64_t ldb,
                          int64_t* sdim, lapack_complex_double* alpha,
                          lapack_complex_double* beta,
                          lapack_complex_double* vsl, int64_t ldvsl,
                          lapack_complex_double* vsr, int64_t ldvsr );

int64_t LAPACKE_sgges3_64( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_S_SELECT3 selctg, int64_t n,
                           float* a, int64_t lda, float* b, int64_t ldb,
                           int64_t* sdim, float* alphar, float* alphai,
                           float* beta, float* vsl, int64_t ldvsl,
                           float* vsr, int64_t ldvsr );
int64_t LAPACKE_dgges3_64( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_D_SELECT3 selctg, int64_t n,
                           double* a, int64_t lda, double* b, int64_t ldb,
                           int64_t* sdim, double* alphar, double* alphai,
                           double* beta, double* vsl, int64_t ldvsl,
                           double* vsr, int64_t ldvsr );
int64_t LAPACKE_cgges3_64( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_C_SELECT2 selctg, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           int64_t* sdim, lapack_complex_float* alpha,
                           lapack_complex_float* beta,
                           lapack_complex_float* vsl, int64_t ldvsl,
                           lapack_complex_float* vsr, int64_t ldvsr );
int64_t LAPACKE_zgges3_64( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_Z_SELECT2 selctg, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           int64_t* sdim, lapack_complex_double* alpha,
                           lapack_complex_double* beta,
                           lapack_complex_double* vsl, int64_t ldvsl,
                           lapack_complex_double* vsr, int64_t ldvsr );

int64_t LAPACKE_sggesx_64( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_S_SELECT3 selctg, char sense,
                           int64_t n, float* a, int64_t lda, float* b,
                           int64_t ldb, int64_t* sdim, float* alphar,
                           float* alphai, float* beta, float* vsl,
                           int64_t ldvsl, float* vsr, int64_t ldvsr,
                           float* rconde, float* rcondv );
int64_t LAPACKE_dggesx_64( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_D_SELECT3 selctg, char sense,
                           int64_t n, double* a, int64_t lda, double* b,
                           int64_t ldb, int64_t* sdim, double* alphar,
                           double* alphai, double* beta, double* vsl,
                           int64_t ldvsl, double* vsr, int64_t ldvsr,
                           double* rconde, double* rcondv );
int64_t LAPACKE_cggesx_64( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_C_SELECT2 selctg, char sense,
                           int64_t n, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb, int64_t* sdim,
                           lapack_complex_float* alpha,
                           lapack_complex_float* beta,
                           lapack_complex_float* vsl, int64_t ldvsl,
                           lapack_complex_float* vsr, int64_t ldvsr,
                           float* rconde, float* rcondv );
int64_t LAPACKE_zggesx_64( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_Z_SELECT2 selctg, char sense,
                           int64_t n, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb, int64_t* sdim,
                           lapack_complex_double* alpha,
                           lapack_complex_double* beta,
                           lapack_complex_double* vsl, int64_t ldvsl,
                           lapack_complex_double* vsr, int64_t ldvsr,
                           double* rconde, double* rcondv );

int64_t LAPACKE_sggev_64( int matrix_layout, char jobvl, char jobvr,
                          int64_t n, float* a, int64_t lda, float* b,
                          int64_t ldb, float* alphar, float* alphai,
                          float* beta, float* vl, int64_t ldvl, float* vr,
                          int64_t ldvr );
int64_t LAPACKE_dggev_64( int matrix_layout, char jobvl, char jobvr,
                          int64_t n, double* a, int64_t lda, double* b,
                          int64_t ldb, double* alphar, double* alphai,
                          double* beta, double* vl, int64_t ldvl, double* vr,
                          int64_t ldvr );
int64_t LAPACKE_cggev_64( int matrix_layout, char jobvl, char jobvr,
                          int64_t n, lapack_complex_float* a, int64_t lda,
                          lapack_complex_float* b, int64_t ldb,
                          lapack_complex_float* alpha,
                          lapack_complex_float* beta, lapack_complex_float* vl,
                          int64_t ldvl, lapack_complex_float* vr,
                          int64_t ldvr );
int64_t LAPACKE_zggev_64( int matrix_layout, char jobvl, char jobvr,
                          int64_t n, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* b,
                          int64_t ldb, lapack_complex_double* alpha,
                          lapack_complex_double* beta,
                          lapack_complex_double* vl, int64_t ldvl,
                          lapack_complex_double* vr, int64_t ldvr );

int64_t LAPACKE_sggev3_64( int matrix_layout, char jobvl, char jobvr,
                           int64_t n, float* a, int64_t lda,
                           float* b, int64_t ldb,
                           float* alphar, float* alphai, float* beta,
                           float* vl, int64_t ldvl,
                           float* vr, int64_t ldvr );
int64_t LAPACKE_dggev3_64( int matrix_layout, char jobvl, char jobvr,
                           int64_t n, double* a, int64_t lda,
                           double* b, int64_t ldb,
                           double* alphar, double* alphai, double* beta,
                           double* vl, int64_t ldvl,
                           double* vr, int64_t ldvr );
int64_t LAPACKE_cggev3_64( int matrix_layout, char jobvl, char jobvr,
                           int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* alpha,
                           lapack_complex_float* beta,
                           lapack_complex_float* vl, int64_t ldvl,
                           lapack_complex_float* vr, int64_t ldvr );
int64_t LAPACKE_zggev3_64( int matrix_layout, char jobvl, char jobvr,
                           int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* alpha,
                           lapack_complex_double* beta,
                           lapack_complex_double* vl, int64_t ldvl,
                           lapack_complex_double* vr, int64_t ldvr );

int64_t LAPACKE_sggevx_64( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, int64_t n, float* a,
                           int64_t lda, float* b, int64_t ldb,
                           float* alphar, float* alphai, float* beta, float* vl,
                           int64_t ldvl, float* vr, int64_t ldvr,
                           int64_t* ilo, int64_t* ihi, float* lscale,
                           float* rscale, float* abnrm, float* bbnrm,
                           float* rconde, float* rcondv );
int64_t LAPACKE_dggevx_64( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, int64_t n, double* a,
                           int64_t lda, double* b, int64_t ldb,
                           double* alphar, double* alphai, double* beta,
                           double* vl, int64_t ldvl, double* vr,
                           int64_t ldvr, int64_t* ilo, int64_t* ihi,
                           double* lscale, double* rscale, double* abnrm,
                           double* bbnrm, double* rconde, double* rcondv );
int64_t LAPACKE_cggevx_64( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* alpha,
                           lapack_complex_float* beta, lapack_complex_float* vl,
                           int64_t ldvl, lapack_complex_float* vr,
                           int64_t ldvr, int64_t* ilo, int64_t* ihi,
                           float* lscale, float* rscale, float* abnrm,
                           float* bbnrm, float* rconde, float* rcondv );
int64_t LAPACKE_zggevx_64( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* alpha,
                           lapack_complex_double* beta,
                           lapack_complex_double* vl, int64_t ldvl,
                           lapack_complex_double* vr, int64_t ldvr,
                           int64_t* ilo, int64_t* ihi, double* lscale,
                           double* rscale, double* abnrm, double* bbnrm,
                           double* rconde, double* rcondv );

int64_t LAPACKE_sggglm_64( int matrix_layout, int64_t n, int64_t m,
                           int64_t p, float* a, int64_t lda, float* b,
                           int64_t ldb, float* d, float* x, float* y );
int64_t LAPACKE_dggglm_64( int matrix_layout, int64_t n, int64_t m,
                           int64_t p, double* a, int64_t lda, double* b,
                           int64_t ldb, double* d, double* x, double* y );
int64_t LAPACKE_cggglm_64( int matrix_layout, int64_t n, int64_t m,
                           int64_t p, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb, lapack_complex_float* d,
                           lapack_complex_float* x, lapack_complex_float* y );
int64_t LAPACKE_zggglm_64( int matrix_layout, int64_t n, int64_t m,
                           int64_t p, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb, lapack_complex_double* d,
                           lapack_complex_double* x, lapack_complex_double* y );

int64_t LAPACKE_sgghrd_64( int matrix_layout, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           float* a, int64_t lda, float* b, int64_t ldb,
                           float* q, int64_t ldq, float* z, int64_t ldz );
int64_t LAPACKE_dgghrd_64( int matrix_layout, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           double* a, int64_t lda, double* b, int64_t ldb,
                           double* q, int64_t ldq, double* z,
                           int64_t ldz );
int64_t LAPACKE_cgghrd_64( int matrix_layout, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* q, int64_t ldq,
                           lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zgghrd_64( int matrix_layout, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* q, int64_t ldq,
                           lapack_complex_double* z, int64_t ldz );

int64_t LAPACKE_sgghd3_64( int matrix_layout, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           float* a, int64_t lda, float* b, int64_t ldb,
                           float* q, int64_t ldq, float* z, int64_t ldz );
int64_t LAPACKE_dgghd3_64( int matrix_layout, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           double* a, int64_t lda, double* b, int64_t ldb,
                           double* q, int64_t ldq, double* z,
                           int64_t ldz );
int64_t LAPACKE_cgghd3_64( int matrix_layout, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* q, int64_t ldq,
                           lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zgghd3_64( int matrix_layout, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* q, int64_t ldq,
                           lapack_complex_double* z, int64_t ldz );

int64_t LAPACKE_sgglse_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t p, float* a, int64_t lda, float* b,
                           int64_t ldb, float* c, float* d, float* x );
int64_t LAPACKE_dgglse_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t p, double* a, int64_t lda, double* b,
                           int64_t ldb, double* c, double* d, double* x );
int64_t LAPACKE_cgglse_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t p, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb, lapack_complex_float* c,
                           lapack_complex_float* d, lapack_complex_float* x );
int64_t LAPACKE_zgglse_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t p, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb, lapack_complex_double* c,
                           lapack_complex_double* d, lapack_complex_double* x );

int64_t LAPACKE_sggqrf_64( int matrix_layout, int64_t n, int64_t m,
                           int64_t p, float* a, int64_t lda, float* taua,
                           float* b, int64_t ldb, float* taub );
int64_t LAPACKE_dggqrf_64( int matrix_layout, int64_t n, int64_t m,
                           int64_t p, double* a, int64_t lda,
                           double* taua, double* b, int64_t ldb,
                           double* taub );
int64_t LAPACKE_cggqrf_64( int matrix_layout, int64_t n, int64_t m,
                           int64_t p, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* taua,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* taub );
int64_t LAPACKE_zggqrf_64( int matrix_layout, int64_t n, int64_t m,
                           int64_t p, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* taua,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* taub );

int64_t LAPACKE_sggrqf_64( int matrix_layout, int64_t m, int64_t p,
                           int64_t n, float* a, int64_t lda, float* taua,
                           float* b, int64_t ldb, float* taub );
int64_t LAPACKE_dggrqf_64( int matrix_layout, int64_t m, int64_t p,
                           int64_t n, double* a, int64_t lda,
                           double* taua, double* b, int64_t ldb,
                           double* taub );
int64_t LAPACKE_cggrqf_64( int matrix_layout, int64_t m, int64_t p,
                           int64_t n, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* taua,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* taub );
int64_t LAPACKE_zggrqf_64( int matrix_layout, int64_t m, int64_t p,
                           int64_t n, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* taua,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* taub );

int64_t LAPACKE_sggsvd_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t n, int64_t p,
                           int64_t* k, int64_t* l, float* a,
                           int64_t lda, float* b, int64_t ldb,
                           float* alpha, float* beta, float* u, int64_t ldu,
                           float* v, int64_t ldv, float* q, int64_t ldq,
                           int64_t* iwork );
int64_t LAPACKE_dggsvd_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t n, int64_t p,
                           int64_t* k, int64_t* l, double* a,
                           int64_t lda, double* b, int64_t ldb,
                           double* alpha, double* beta, double* u,
                           int64_t ldu, double* v, int64_t ldv, double* q,
                           int64_t ldq, int64_t* iwork );
int64_t LAPACKE_cggsvd_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t n, int64_t p,
                           int64_t* k, int64_t* l,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           float* alpha, float* beta, lapack_complex_float* u,
                           int64_t ldu, lapack_complex_float* v,
                           int64_t ldv, lapack_complex_float* q,
                           int64_t ldq, int64_t* iwork );
int64_t LAPACKE_zggsvd_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t n, int64_t p,
                           int64_t* k, int64_t* l,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           double* alpha, double* beta,
                           lapack_complex_double* u, int64_t ldu,
                           lapack_complex_double* v, int64_t ldv,
                           lapack_complex_double* q, int64_t ldq,
                           int64_t* iwork );

int64_t LAPACKE_sggsvd3_64( int matrix_layout, char jobu, char jobv, char jobq,
                            int64_t m, int64_t n, int64_t p,
                            int64_t* k, int64_t* l, float* a,
                            int64_t lda, float* b, int64_t ldb,
                            float* alpha, float* beta, float* u, int64_t ldu,
                            float* v, int64_t ldv, float* q, int64_t ldq,
                            int64_t* iwork );
int64_t LAPACKE_dggsvd3_64( int matrix_layout, char jobu, char jobv, char jobq,
                            int64_t m, int64_t n, int64_t p,
                            int64_t* k, int64_t* l, double* a,
                            int64_t lda, double* b, int64_t ldb,
                            double* alpha, double* beta, double* u,
                            int64_t ldu, double* v, int64_t ldv, double* q,
                            int64_t ldq, int64_t* iwork );
int64_t LAPACKE_cggsvd3_64( int matrix_layout, char jobu, char jobv, char jobq,
                            int64_t m, int64_t n, int64_t p,
                            int64_t* k, int64_t* l,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* b, int64_t ldb,
                            float* alpha, float* beta, lapack_complex_float* u,
                            int64_t ldu, lapack_complex_float* v,
                            int64_t ldv, lapack_complex_float* q,
                            int64_t ldq, int64_t* iwork );
int64_t LAPACKE_zggsvd3_64( int matrix_layout, char jobu, char jobv, char jobq,
                            int64_t m, int64_t n, int64_t p,
                            int64_t* k, int64_t* l,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* b, int64_t ldb,
                            double* alpha, double* beta,
                            lapack_complex_double* u, int64_t ldu,
                            lapack_complex_double* v, int64_t ldv,
                            lapack_complex_double* q, int64_t ldq,
                            int64_t* iwork );

int64_t LAPACKE_sggsvp_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t p, int64_t n, float* a,
                           int64_t lda, float* b, int64_t ldb, float tola,
                           float tolb, int64_t* k, int64_t* l, float* u,
                           int64_t ldu, float* v, int64_t ldv, float* q,
                           int64_t ldq );
int64_t LAPACKE_dggsvp_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t p, int64_t n, double* a,
                           int64_t lda, double* b, int64_t ldb,
                           double tola, double tolb, int64_t* k,
                           int64_t* l, double* u, int64_t ldu, double* v,
                           int64_t ldv, double* q, int64_t ldq );
int64_t LAPACKE_cggsvp_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t p, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb, float tola,
                           float tolb, int64_t* k, int64_t* l,
                           lapack_complex_float* u, int64_t ldu,
                           lapack_complex_float* v, int64_t ldv,
                           lapack_complex_float* q, int64_t ldq );
int64_t LAPACKE_zggsvp_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t p, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           double tola, double tolb, int64_t* k,
                           int64_t* l, lapack_complex_double* u,
                           int64_t ldu, lapack_complex_double* v,
                           int64_t ldv, lapack_complex_double* q,
                           int64_t ldq );

int64_t LAPACKE_sggsvp3_64( int matrix_layout, char jobu, char jobv, char jobq,
                            int64_t m, int64_t p, int64_t n, float* a,
                            int64_t lda, float* b, int64_t ldb, float tola,
                            float tolb, int64_t* k, int64_t* l, float* u,
                            int64_t ldu, float* v, int64_t ldv, float* q,
                            int64_t ldq );
int64_t LAPACKE_dggsvp3_64( int matrix_layout, char jobu, char jobv, char jobq,
                            int64_t m, int64_t p, int64_t n, double* a,
                            int64_t lda, double* b, int64_t ldb,
                            double tola, double tolb, int64_t* k,
                            int64_t* l, double* u, int64_t ldu, double* v,
                            int64_t ldv, double* q, int64_t ldq );
int64_t LAPACKE_cggsvp3_64( int matrix_layout, char jobu, char jobv, char jobq,
                            int64_t m, int64_t p, int64_t n,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* b, int64_t ldb, float tola,
                            float tolb, int64_t* k, int64_t* l,
                            lapack_complex_float* u, int64_t ldu,
                            lapack_complex_float* v, int64_t ldv,
                            lapack_complex_float* q, int64_t ldq );
int64_t LAPACKE_zggsvp3_64( int matrix_layout, char jobu, char jobv, char jobq,
                            int64_t m, int64_t p, int64_t n,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* b, int64_t ldb,
                            double tola, double tolb, int64_t* k,
                            int64_t* l, lapack_complex_double* u,
                            int64_t ldu, lapack_complex_double* v,
                            int64_t ldv, lapack_complex_double* q,
                            int64_t ldq );

int64_t LAPACKE_sgtcon_64( char norm, int64_t n, const float* dl,
                           const float* d, const float* du, const float* du2,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_dgtcon_64( char norm, int64_t n, const double* dl,
                           const double* d, const double* du, const double* du2,
                           const int64_t* ipiv, double anorm,
                           double* rcond );
int64_t LAPACKE_cgtcon_64( char norm, int64_t n,
                           const lapack_complex_float* dl,
                           const lapack_complex_float* d,
                           const lapack_complex_float* du,
                           const lapack_complex_float* du2,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_zgtcon_64( char norm, int64_t n,
                           const lapack_complex_double* dl,
                           const lapack_complex_double* d,
                           const lapack_complex_double* du,
                           const lapack_complex_double* du2,
                           const int64_t* ipiv, double anorm,
                           double* rcond );

int64_t LAPACKE_sgtrfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const float* dl, const float* d,
                           const float* du, const float* dlf, const float* df,
                           const float* duf, const float* du2,
                           const int64_t* ipiv, const float* b,
                           int64_t ldb, float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_dgtrfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const double* dl, const double* d,
                           const double* du, const double* dlf,
                           const double* df, const double* duf,
                           const double* du2, const int64_t* ipiv,
                           const double* b, int64_t ldb, double* x,
                           int64_t ldx, double* ferr, double* berr );
int64_t LAPACKE_cgtrfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const lapack_complex_float* dl,
                           const lapack_complex_float* d,
                           const lapack_complex_float* du,
                           const lapack_complex_float* dlf,
                           const lapack_complex_float* df,
                           const lapack_complex_float* duf,
                           const lapack_complex_float* du2,
                           const int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zgtrfs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const lapack_complex_double* dl,
                           const lapack_complex_double* d,
                           const lapack_complex_double* du,
                           const lapack_complex_double* dlf,
                           const lapack_complex_double* df,
                           const lapack_complex_double* duf,
                           const lapack_complex_double* du2,
                           const int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_sgtsv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          float* dl, float* d, float* du, float* b,
                          int64_t ldb );
int64_t LAPACKE_dgtsv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          double* dl, double* d, double* du, double* b,
                          int64_t ldb );
int64_t LAPACKE_cgtsv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          lapack_complex_float* dl, lapack_complex_float* d,
                          lapack_complex_float* du, lapack_complex_float* b,
                          int64_t ldb );
int64_t LAPACKE_zgtsv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          lapack_complex_double* dl, lapack_complex_double* d,
                          lapack_complex_double* du, lapack_complex_double* b,
                          int64_t ldb );

int64_t LAPACKE_sgtsvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t nrhs, const float* dl,
                           const float* d, const float* du, float* dlf,
                           float* df, float* duf, float* du2, int64_t* ipiv,
                           const float* b, int64_t ldb, float* x,
                           int64_t ldx, float* rcond, float* ferr,
                           float* berr );
int64_t LAPACKE_dgtsvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t nrhs, const double* dl,
                           const double* d, const double* du, double* dlf,
                           double* df, double* duf, double* du2,
                           int64_t* ipiv, const double* b, int64_t ldb,
                           double* x, int64_t ldx, double* rcond,
                           double* ferr, double* berr );
int64_t LAPACKE_cgtsvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_float* dl,
                           const lapack_complex_float* d,
                           const lapack_complex_float* du,
                           lapack_complex_float* dlf, lapack_complex_float* df,
                           lapack_complex_float* duf, lapack_complex_float* du2,
                           int64_t* ipiv, const lapack_complex_float* b,
                           int64_t ldb, lapack_complex_float* x,
                           int64_t ldx, float* rcond, float* ferr,
                           float* berr );
int64_t LAPACKE_zgtsvx_64( int matrix_layout, char fact, char trans,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_double* dl,
                           const lapack_complex_double* d,
                           const lapack_complex_double* du,
                           lapack_complex_double* dlf,
                           lapack_complex_double* df,
                           lapack_complex_double* duf,
                           lapack_complex_double* du2, int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );

int64_t LAPACKE_sgttrf_64( int64_t n, float* dl, float* d, float* du,
                           float* du2, int64_t* ipiv );
int64_t LAPACKE_dgttrf_64( int64_t n, double* dl, double* d, double* du,
                           double* du2, int64_t* ipiv );
int64_t LAPACKE_cgttrf_64( int64_t n, lapack_complex_float* dl,
                           lapack_complex_float* d, lapack_complex_float* du,
                           lapack_complex_float* du2, int64_t* ipiv );
int64_t LAPACKE_zgttrf_64( int64_t n, lapack_complex_double* dl,
                           lapack_complex_double* d, lapack_complex_double* du,
                           lapack_complex_double* du2, int64_t* ipiv );

int64_t LAPACKE_sgttrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const float* dl, const float* d,
                           const float* du, const float* du2,
                           const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dgttrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const double* dl, const double* d,
                           const double* du, const double* du2,
                           const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_cgttrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const lapack_complex_float* dl,
                           const lapack_complex_float* d,
                           const lapack_complex_float* du,
                           const lapack_complex_float* du2,
                           const int64_t* ipiv, lapack_complex_float* b,
                           int64_t ldb );
int64_t LAPACKE_zgttrs_64( int matrix_layout, char trans, int64_t n,
                           int64_t nrhs, const lapack_complex_double* dl,
                           const lapack_complex_double* d,
                           const lapack_complex_double* du,
                           const lapack_complex_double* du2,
                           const int64_t* ipiv, lapack_complex_double* b,
                           int64_t ldb );

int64_t LAPACKE_chbev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t kd, lapack_complex_float* ab,
                          int64_t ldab, float* w, lapack_complex_float* z,
                          int64_t ldz );
int64_t LAPACKE_zhbev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t kd, lapack_complex_double* ab,
                          int64_t ldab, double* w, lapack_complex_double* z,
                          int64_t ldz );

int64_t LAPACKE_chbevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t kd, lapack_complex_float* ab,
                           int64_t ldab, float* w, lapack_complex_float* z,
                           int64_t ldz );
int64_t LAPACKE_zhbevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t kd, lapack_complex_double* ab,
                           int64_t ldab, double* w, lapack_complex_double* z,
                           int64_t ldz );

int64_t LAPACKE_chbevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t kd,
                           lapack_complex_float* ab, int64_t ldab,
                           lapack_complex_float* q, int64_t ldq, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, lapack_complex_float* z,
                           int64_t ldz, int64_t* ifail );
int64_t LAPACKE_zhbevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t kd,
                           lapack_complex_double* ab, int64_t ldab,
                           lapack_complex_double* q, int64_t ldq, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_chbgst_64( int matrix_layout, char vect, char uplo, int64_t n,
                           int64_t ka, int64_t kb,
                           lapack_complex_float* ab, int64_t ldab,
                           const lapack_complex_float* bb, int64_t ldbb,
                           lapack_complex_float* x, int64_t ldx );
int64_t LAPACKE_zhbgst_64( int matrix_layout, char vect, char uplo, int64_t n,
                           int64_t ka, int64_t kb,
                           lapack_complex_double* ab, int64_t ldab,
                           const lapack_complex_double* bb, int64_t ldbb,
                           lapack_complex_double* x, int64_t ldx );

int64_t LAPACKE_chbgv_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t ka, int64_t kb,
                          lapack_complex_float* ab, int64_t ldab,
                          lapack_complex_float* bb, int64_t ldbb, float* w,
                          lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zhbgv_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t ka, int64_t kb,
                          lapack_complex_double* ab, int64_t ldab,
                          lapack_complex_double* bb, int64_t ldbb, double* w,
                          lapack_complex_double* z, int64_t ldz );

int64_t LAPACKE_chbgvd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t ka, int64_t kb,
                           lapack_complex_float* ab, int64_t ldab,
                           lapack_complex_float* bb, int64_t ldbb, float* w,
                           lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zhbgvd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t ka, int64_t kb,
                           lapack_complex_double* ab, int64_t ldab,
                           lapack_complex_double* bb, int64_t ldbb,
                           double* w, lapack_complex_double* z,
                           int64_t ldz );

int64_t LAPACKE_chbgvx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t ka, int64_t kb,
                           lapack_complex_float* ab, int64_t ldab,
                           lapack_complex_float* bb, int64_t ldbb,
                           lapack_complex_float* q, int64_t ldq, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, lapack_complex_float* z,
                           int64_t ldz, int64_t* ifail );
int64_t LAPACKE_zhbgvx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t ka, int64_t kb,
                           lapack_complex_double* ab, int64_t ldab,
                           lapack_complex_double* bb, int64_t ldbb,
                           lapack_complex_double* q, int64_t ldq, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_chbtrd_64( int matrix_layout, char vect, char uplo, int64_t n,
                           int64_t kd, lapack_complex_float* ab,
                           int64_t ldab, float* d, float* e,
                           lapack_complex_float* q, int64_t ldq );
int64_t LAPACKE_zhbtrd_64( int matrix_layout, char vect, char uplo, int64_t n,
                           int64_t kd, lapack_complex_double* ab,
                           int64_t ldab, double* d, double* e,
                           lapack_complex_double* q, int64_t ldq );

int64_t LAPACKE_checon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_zhecon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           const int64_t* ipiv, double anorm,
                           double* rcond );

int64_t LAPACKE_cheequb_64( int matrix_layout, char uplo, int64_t n,
                            const lapack_complex_float* a, int64_t lda,
                            float* s, float* scond, float* amax );
int64_t LAPACKE_zheequb_64( int matrix_layout, char uplo, int64_t n,
                            const lapack_complex_double* a, int64_t lda,
                            double* s, double* scond, double* amax );

int64_t LAPACKE_cheev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          lapack_complex_float* a, int64_t lda, float* w );
int64_t LAPACKE_zheev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          lapack_complex_double* a, int64_t lda, double* w );

int64_t LAPACKE_cheevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda, float* w );
int64_t LAPACKE_zheevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           double* w );

int64_t LAPACKE_cheevr_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_float* a,
                           int64_t lda, float vl, float vu, int64_t il,
                           int64_t iu, float abstol, int64_t* m, float* w,
                           lapack_complex_float* z, int64_t ldz,
                           int64_t* isuppz );
int64_t LAPACKE_zheevr_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_double* a,
                           int64_t lda, double vl, double vu, int64_t il,
                           int64_t iu, double abstol, int64_t* m,
                           double* w, lapack_complex_double* z, int64_t ldz,
                           int64_t* isuppz );

int64_t LAPACKE_cheevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_float* a,
                           int64_t lda, float vl, float vu, int64_t il,
                           int64_t iu, float abstol, int64_t* m, float* w,
                           lapack_complex_float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_zheevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_double* a,
                           int64_t lda, double vl, double vu, int64_t il,
                           int64_t iu, double abstol, int64_t* m,
                           double* w, lapack_complex_double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_chegst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* b,
                           int64_t ldb );
int64_t LAPACKE_zhegst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* b,
                           int64_t ldb );

int64_t LAPACKE_chegv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* b,
                          int64_t ldb, float* w );
int64_t LAPACKE_zhegv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* b,
                          int64_t ldb, double* w );

int64_t LAPACKE_chegvd_64( int matrix_layout, int64_t itype, char jobz,
                           char uplo, int64_t n, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb, float* w );
int64_t LAPACKE_zhegvd_64( int matrix_layout, int64_t itype, char jobz,
                           char uplo, int64_t n, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb, double* w );

int64_t LAPACKE_chegvx_64( int matrix_layout, int64_t itype, char jobz,
                           char range, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, lapack_complex_float* z,
                           int64_t ldz, int64_t* ifail );
int64_t LAPACKE_zhegvx_64( int matrix_layout, int64_t itype, char jobz,
                           char range, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_cherfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* af,
                           int64_t ldaf, const int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zherfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* af,
                           int64_t ldaf, const int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_cherfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs,
                            const lapack_complex_float* a, int64_t lda,
                            const lapack_complex_float* af, int64_t ldaf,
                            const int64_t* ipiv, const float* s,
                            const lapack_complex_float* b, int64_t ldb,
                            lapack_complex_float* x, int64_t ldx,
                            float* rcond, float* berr, int64_t n_err_bnds,
                            float* err_bnds_norm, float* err_bnds_comp,
                            int64_t nparams, float* params );
int64_t LAPACKE_zherfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs,
                            const lapack_complex_double* a, int64_t lda,
                            const lapack_complex_double* af, int64_t ldaf,
                            const int64_t* ipiv, const double* s,
                            const lapack_complex_double* b, int64_t ldb,
                            lapack_complex_double* x, int64_t ldx,
                            double* rcond, double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );

int64_t LAPACKE_chesv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, int64_t* ipiv,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zhesv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, int64_t* ipiv,
                          lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_chesvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* af,
                           int64_t ldaf, int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr );
int64_t LAPACKE_zhesvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* af,
                           int64_t ldaf, int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );

int64_t LAPACKE_chesvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, float* s,
                            lapack_complex_float* b, int64_t ldb,
                            lapack_complex_float* x, int64_t ldx,
                            float* rcond, float* rpvgrw, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_zhesvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, double* s,
                            lapack_complex_double* b, int64_t ldb,
                            lapack_complex_double* x, int64_t ldx,
                            double* rcond, double* rpvgrw, double* berr,
                            int64_t n_err_bnds, double* err_bnds_norm,
                            double* err_bnds_comp, int64_t nparams,
                            double* params );

int64_t LAPACKE_chetrd_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda, float* d,
                           float* e, lapack_complex_float* tau );
int64_t LAPACKE_zhetrd_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda, double* d,
                           double* e, lapack_complex_double* tau );

int64_t LAPACKE_chetrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_zhetrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv );

int64_t LAPACKE_chetri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           const int64_t* ipiv );
int64_t LAPACKE_zhetri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           const int64_t* ipiv );

int64_t LAPACKE_chetrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zhetrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_chfrk_64( int matrix_layout, char transr, char uplo, char trans,
                          int64_t n, int64_t k, float alpha,
                          const lapack_complex_float* a, int64_t lda,
                          float beta, lapack_complex_float* c );
int64_t LAPACKE_zhfrk_64( int matrix_layout, char transr, char uplo, char trans,
                          int64_t n, int64_t k, double alpha,
                          const lapack_complex_double* a, int64_t lda,
                          double beta, lapack_complex_double* c );

int64_t LAPACKE_shgeqz_64( int matrix_layout, char job, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           float* h, int64_t ldh, float* t, int64_t ldt,
                           float* alphar, float* alphai, float* beta, float* q,
                           int64_t ldq, float* z, int64_t ldz );
int64_t LAPACKE_dhgeqz_64( int matrix_layout, char job, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           double* h, int64_t ldh, double* t, int64_t ldt,
                           double* alphar, double* alphai, double* beta,
                           double* q, int64_t ldq, double* z,
                           int64_t ldz );
int64_t LAPACKE_chgeqz_64( int matrix_layout, char job, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           lapack_complex_float* h, int64_t ldh,
                           lapack_complex_float* t, int64_t ldt,
                           lapack_complex_float* alpha,
                           lapack_complex_float* beta, lapack_complex_float* q,
                           int64_t ldq, lapack_complex_float* z,
                           int64_t ldz );
int64_t LAPACKE_zhgeqz_64( int matrix_layout, char job, char compq, char compz,
                           int64_t n, int64_t ilo, int64_t ihi,
                           lapack_complex_double* h, int64_t ldh,
                           lapack_complex_double* t, int64_t ldt,
                           lapack_complex_double* alpha,
                           lapack_complex_double* beta,
                           lapack_complex_double* q, int64_t ldq,
                           lapack_complex_double* z, int64_t ldz );

int64_t LAPACKE_chpcon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* ap,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_zhpcon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* ap,
                           const int64_t* ipiv, double anorm,
                           double* rcond );

int64_t LAPACKE_chpev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          lapack_complex_float* ap, float* w,
                          lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zhpev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          lapack_complex_double* ap, double* w,
                          lapack_complex_double* z, int64_t ldz );

int64_t LAPACKE_chpevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           lapack_complex_float* ap, float* w,
                           lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zhpevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           lapack_complex_double* ap, double* w,
                           lapack_complex_double* z, int64_t ldz );

int64_t LAPACKE_chpevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_float* ap, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, lapack_complex_float* z,
                           int64_t ldz, int64_t* ifail );
int64_t LAPACKE_zhpevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_double* ap, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_chpgst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, lapack_complex_float* ap,
                           const lapack_complex_float* bp );
int64_t LAPACKE_zhpgst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, lapack_complex_double* ap,
                           const lapack_complex_double* bp );

int64_t LAPACKE_chpgv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, lapack_complex_float* ap,
                          lapack_complex_float* bp, float* w,
                          lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zhpgv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, lapack_complex_double* ap,
                          lapack_complex_double* bp, double* w,
                          lapack_complex_double* z, int64_t ldz );

int64_t LAPACKE_chpgvd_64( int matrix_layout, int64_t itype, char jobz,
                           char uplo, int64_t n, lapack_complex_float* ap,
                           lapack_complex_float* bp, float* w,
                           lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zhpgvd_64( int matrix_layout, int64_t itype, char jobz,
                           char uplo, int64_t n, lapack_complex_double* ap,
                           lapack_complex_double* bp, double* w,
                           lapack_complex_double* z, int64_t ldz );

int64_t LAPACKE_chpgvx_64( int matrix_layout, int64_t itype, char jobz,
                           char range, char uplo, int64_t n,
                           lapack_complex_float* ap, lapack_complex_float* bp,
                           float vl, float vu, int64_t il, int64_t iu,
                           float abstol, int64_t* m, float* w,
                           lapack_complex_float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_zhpgvx_64( int matrix_layout, int64_t itype, char jobz,
                           char range, char uplo, int64_t n,
                           lapack_complex_double* ap, lapack_complex_double* bp,
                           double vl, double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_chprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* ap,
                           const lapack_complex_float* afp,
                           const int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zhprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* ap,
                           const lapack_complex_double* afp,
                           const int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_chpsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* ap,
                          int64_t* ipiv, lapack_complex_float* b,
                          int64_t ldb );
int64_t LAPACKE_zhpsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* ap,
                          int64_t* ipiv, lapack_complex_double* b,
                          int64_t ldb );

int64_t LAPACKE_chpsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* ap,
                           lapack_complex_float* afp, int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr );
int64_t LAPACKE_zhpsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* ap,
                           lapack_complex_double* afp, int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );

int64_t LAPACKE_chptrd_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* ap, float* d, float* e,
                           lapack_complex_float* tau );
int64_t LAPACKE_zhptrd_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* ap, double* d, double* e,
                           lapack_complex_double* tau );

int64_t LAPACKE_chptrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* ap, int64_t* ipiv );
int64_t LAPACKE_zhptrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* ap, int64_t* ipiv );

int64_t LAPACKE_chptri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* ap, const int64_t* ipiv );
int64_t LAPACKE_zhptri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* ap, const int64_t* ipiv );

int64_t LAPACKE_chptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* ap,
                           const int64_t* ipiv, lapack_complex_float* b,
                           int64_t ldb );
int64_t LAPACKE_zhptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* ap,
                           const int64_t* ipiv, lapack_complex_double* b,
                           int64_t ldb );

int64_t LAPACKE_shsein_64( int matrix_layout, char job, char eigsrc, char initv,
                           lapack_logical* select, int64_t n, const float* h,
                           int64_t ldh, float* wr, const float* wi,
                           float* vl, int64_t ldvl, float* vr,
                           int64_t ldvr, int64_t mm, int64_t* m,
                           int64_t* ifaill, int64_t* ifailr );
int64_t LAPACKE_dhsein_64( int matrix_layout, char job, char eigsrc, char initv,
                           lapack_logical* select, int64_t n,
                           const double* h, int64_t ldh, double* wr,
                           const double* wi, double* vl, int64_t ldvl,
                           double* vr, int64_t ldvr, int64_t mm,
                           int64_t* m, int64_t* ifaill,
                           int64_t* ifailr );
int64_t LAPACKE_chsein_64( int matrix_layout, char job, char eigsrc, char initv,
                           const lapack_logical* select, int64_t n,
                           const lapack_complex_float* h, int64_t ldh,
                           lapack_complex_float* w, lapack_complex_float* vl,
                           int64_t ldvl, lapack_complex_float* vr,
                           int64_t ldvr, int64_t mm, int64_t* m,
                           int64_t* ifaill, int64_t* ifailr );
int64_t LAPACKE_zhsein_64( int matrix_layout, char job, char eigsrc, char initv,
                           const lapack_logical* select, int64_t n,
                           const lapack_complex_double* h, int64_t ldh,
                           lapack_complex_double* w, lapack_complex_double* vl,
                           int64_t ldvl, lapack_complex_double* vr,
                           int64_t ldvr, int64_t mm, int64_t* m,
                           int64_t* ifaill, int64_t* ifailr );

int64_t LAPACKE_shseqr_64( int matrix_layout, char job, char compz, int64_t n,
                           int64_t ilo, int64_t ihi, float* h,
                           int64_t ldh, float* wr, float* wi, float* z,
                           int64_t ldz );
int64_t LAPACKE_dhseqr_64( int matrix_layout, char job, char compz, int64_t n,
                           int64_t ilo, int64_t ihi, double* h,
                           int64_t ldh, double* wr, double* wi, double* z,
                           int64_t ldz );
int64_t LAPACKE_chseqr_64( int matrix_layout, char job, char compz, int64_t n,
                           int64_t ilo, int64_t ihi,
                           lapack_complex_float* h, int64_t ldh,
                           lapack_complex_float* w, lapack_complex_float* z,
                           int64_t ldz );
int64_t LAPACKE_zhseqr_64( int matrix_layout, char job, char compz, int64_t n,
                           int64_t ilo, int64_t ihi,
                           lapack_complex_double* h, int64_t ldh,
                           lapack_complex_double* w, lapack_complex_double* z,
                           int64_t ldz );

int64_t LAPACKE_clacgv_64( int64_t n, lapack_complex_float* x,
                           int64_t incx );
int64_t LAPACKE_zlacgv_64( int64_t n, lapack_complex_double* x,
                           int64_t incx );

int64_t LAPACKE_slacn2_64( int64_t n, float* v, float* x, int64_t* isgn,
                           float* est, int64_t* kase, int64_t* isave );
int64_t LAPACKE_dlacn2_64( int64_t n, double* v, double* x, int64_t* isgn,
                           double* est, int64_t* kase, int64_t* isave );
int64_t LAPACKE_clacn2_64( int64_t n, lapack_complex_float* v,
                           lapack_complex_float* x,
                           float* est, int64_t* kase, int64_t* isave );
int64_t LAPACKE_zlacn2_64( int64_t n, lapack_complex_double* v,
                           lapack_complex_double* x,
                           double* est, int64_t* kase, int64_t* isave );

int64_t LAPACKE_slacpy_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, const float* a, int64_t lda, float* b,
                           int64_t ldb );
int64_t LAPACKE_dlacpy_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, const double* a, int64_t lda, double* b,
                           int64_t ldb );
int64_t LAPACKE_clacpy_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, const lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb );
int64_t LAPACKE_zlacpy_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, const lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb );

int64_t LAPACKE_clacp2_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, const float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zlacp2_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, const double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_zlag2c_64( int matrix_layout, int64_t m, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           lapack_complex_float* sa, int64_t ldsa );

int64_t LAPACKE_slag2d_64( int matrix_layout, int64_t m, int64_t n,
                           const float* sa, int64_t ldsa, double* a,
                           int64_t lda );

int64_t LAPACKE_dlag2s_64( int matrix_layout, int64_t m, int64_t n,
                           const double* a, int64_t lda, float* sa,
                           int64_t ldsa );

int64_t LAPACKE_clag2z_64( int matrix_layout, int64_t m, int64_t n,
                           const lapack_complex_float* sa, int64_t ldsa,
                           lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_slagge_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku, const float* d,
                           float* a, int64_t lda, int64_t* iseed );
int64_t LAPACKE_dlagge_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku, const double* d,
                           double* a, int64_t lda, int64_t* iseed );
int64_t LAPACKE_clagge_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku, const float* d,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* iseed );
int64_t LAPACKE_zlagge_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t kl, int64_t ku, const double* d,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* iseed );

float LAPACKE_slamch_64( char cmach );
double LAPACKE_dlamch_64( char cmach );

float LAPACKE_slangb_64( int matrix_layout, char norm, int64_t n,
                      int64_t kl, int64_t ku, const float* ab,
                      int64_t ldab );
double LAPACKE_dlangb_64( int matrix_layout, char norm, int64_t n,
                       int64_t kl, int64_t ku, const double* ab,
                       int64_t ldab );
float LAPACKE_clangb_64( int matrix_layout, char norm, int64_t n,
                      int64_t kl, int64_t ku,
                      const lapack_complex_float* ab, int64_t ldab );
double LAPACKE_zlangb_64( int matrix_layout, char norm, int64_t n,
                       int64_t kl, int64_t ku,
                       const lapack_complex_double* ab, int64_t ldab );

float LAPACKE_slange_64( int matrix_layout, char norm, int64_t m,
                           int64_t n, const float* a, int64_t lda );
double LAPACKE_dlange_64( int matrix_layout, char norm, int64_t m,
                           int64_t n, const double* a, int64_t lda );
float LAPACKE_clange_64( int matrix_layout, char norm, int64_t m,
                           int64_t n, const lapack_complex_float* a,
                           int64_t lda );
double LAPACKE_zlange_64( int matrix_layout, char norm, int64_t m,
                           int64_t n, const lapack_complex_double* a,
                           int64_t lda );

float LAPACKE_clanhe_64( int matrix_layout, char norm, char uplo, int64_t n,
                           const lapack_complex_float* a, int64_t lda );
double LAPACKE_zlanhe_64( int matrix_layout, char norm, char uplo, int64_t n,
                           const lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_clacrm_64( int matrix_layout, int64_t m, int64_t n,
                          const lapack_complex_float* a,
                          int64_t lda, const float* b,
                          int64_t ldb, lapack_complex_float* c,
                          int64_t ldc );
int64_t LAPACKE_zlacrm_64( int matrix_layout, int64_t m, int64_t n,
                           const lapack_complex_double* a,
                           int64_t lda, const double* b,
                           int64_t ldb, lapack_complex_double* c,
                           int64_t ldc );

int64_t LAPACKE_clarcm_64( int matrix_layout, int64_t m, int64_t n,
                          const float* a, int64_t lda,
                          const lapack_complex_float* b,
                          int64_t ldb, lapack_complex_float* c,
                          int64_t ldc );
int64_t LAPACKE_zlarcm_64( int matrix_layout, int64_t m, int64_t n,
                           const double* a, int64_t lda,
                           const lapack_complex_double* b,
                           int64_t ldb, lapack_complex_double* c,
                           int64_t ldc );

float LAPACKE_slansy_64( int matrix_layout, char norm, char uplo, int64_t n,
                           const float* a, int64_t lda );
double LAPACKE_dlansy_64( int matrix_layout, char norm, char uplo, int64_t n,
                           const double* a, int64_t lda );
float LAPACKE_clansy_64( int matrix_layout, char norm, char uplo, int64_t n,
                           const lapack_complex_float* a, int64_t lda );
double LAPACKE_zlansy_64( int matrix_layout, char norm, char uplo, int64_t n,
                           const lapack_complex_double* a, int64_t lda );

float LAPACKE_slanky_64( int matrix_layout, char norm, char uplo, int64_t n,
                           const float* a, int64_t lda );
double LAPACKE_dlanky_64( int matrix_layout, char norm, char uplo, int64_t n,
                           const double* a, int64_t lda );

float LAPACKE_slantr_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t m, int64_t n, const float* a,
                           int64_t lda );
double LAPACKE_dlantr_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t m, int64_t n, const double* a,
                           int64_t lda );
float LAPACKE_clantr_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t m, int64_t n, const lapack_complex_float* a,
                           int64_t lda );
double LAPACKE_zlantr_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t m, int64_t n, const lapack_complex_double* a,
                           int64_t lda );


int64_t LAPACKE_slarfb_64( int matrix_layout, char side, char trans, char direct,
                           char storev, int64_t m, int64_t n,
                           int64_t k, const float* v, int64_t ldv,
                           const float* t, int64_t ldt, float* c,
                           int64_t ldc );
int64_t LAPACKE_dlarfb_64( int matrix_layout, char side, char trans, char direct,
                           char storev, int64_t m, int64_t n,
                           int64_t k, const double* v, int64_t ldv,
                           const double* t, int64_t ldt, double* c,
                           int64_t ldc );
int64_t LAPACKE_clarfb_64( int matrix_layout, char side, char trans, char direct,
                           char storev, int64_t m, int64_t n,
                           int64_t k, const lapack_complex_float* v,
                           int64_t ldv, const lapack_complex_float* t,
                           int64_t ldt, lapack_complex_float* c,
                           int64_t ldc );
int64_t LAPACKE_zlarfb_64( int matrix_layout, char side, char trans, char direct,
                           char storev, int64_t m, int64_t n,
                           int64_t k, const lapack_complex_double* v,
                           int64_t ldv, const lapack_complex_double* t,
                           int64_t ldt, lapack_complex_double* c,
                           int64_t ldc );

int64_t LAPACKE_slarfg_64( int64_t n, float* alpha, float* x,
                           int64_t incx, float* tau );
int64_t LAPACKE_dlarfg_64( int64_t n, double* alpha, double* x,
                           int64_t incx, double* tau );
int64_t LAPACKE_clarfg_64( int64_t n, lapack_complex_float* alpha,
                           lapack_complex_float* x, int64_t incx,
                           lapack_complex_float* tau );
int64_t LAPACKE_zlarfg_64( int64_t n, lapack_complex_double* alpha,
                           lapack_complex_double* x, int64_t incx,
                           lapack_complex_double* tau );

int64_t LAPACKE_slarft_64( int matrix_layout, char direct, char storev,
                           int64_t n, int64_t k, const float* v,
                           int64_t ldv, const float* tau, float* t,
                           int64_t ldt );
int64_t LAPACKE_dlarft_64( int matrix_layout, char direct, char storev,
                           int64_t n, int64_t k, const double* v,
                           int64_t ldv, const double* tau, double* t,
                           int64_t ldt );
int64_t LAPACKE_clarft_64( int matrix_layout, char direct, char storev,
                           int64_t n, int64_t k,
                           const lapack_complex_float* v, int64_t ldv,
                           const lapack_complex_float* tau,
                           lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_zlarft_64( int matrix_layout, char direct, char storev,
                           int64_t n, int64_t k,
                           const lapack_complex_double* v, int64_t ldv,
                           const lapack_complex_double* tau,
                           lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_slarfx_64( int matrix_layout, char side, int64_t m,
                           int64_t n, const float* v, float tau, float* c,
                           int64_t ldc, float* work );
int64_t LAPACKE_dlarfx_64( int matrix_layout, char side, int64_t m,
                           int64_t n, const double* v, double tau, double* c,
                           int64_t ldc, double* work );
int64_t LAPACKE_clarfx_64( int matrix_layout, char side, int64_t m,
                           int64_t n, const lapack_complex_float* v,
                           lapack_complex_float tau, lapack_complex_float* c,
                           int64_t ldc, lapack_complex_float* work );
int64_t LAPACKE_zlarfx_64( int matrix_layout, char side, int64_t m,
                           int64_t n, const lapack_complex_double* v,
                           lapack_complex_double tau, lapack_complex_double* c,
                           int64_t ldc, lapack_complex_double* work );

int64_t LAPACKE_slarnv_64( int64_t idist, int64_t* iseed, int64_t n,
                           float* x );
int64_t LAPACKE_dlarnv_64( int64_t idist, int64_t* iseed, int64_t n,
                           double* x );
int64_t LAPACKE_clarnv_64( int64_t idist, int64_t* iseed, int64_t n,
                           lapack_complex_float* x );
int64_t LAPACKE_zlarnv_64( int64_t idist, int64_t* iseed, int64_t n,
                           lapack_complex_double* x );

int64_t LAPACKE_slascl_64( int matrix_layout, char type, int64_t kl,
                           int64_t ku, float cfrom, float cto,
                           int64_t m, int64_t n, float* a,
                           int64_t lda );
int64_t LAPACKE_dlascl_64( int matrix_layout, char type, int64_t kl,
                           int64_t ku, double cfrom, double cto,
                           int64_t m, int64_t n, double* a,
                           int64_t lda );
int64_t LAPACKE_clascl_64( int matrix_layout, char type, int64_t kl,
                           int64_t ku, float cfrom, float cto,
                           int64_t m, int64_t n, lapack_complex_float* a,
                           int64_t lda );
int64_t LAPACKE_zlascl_64( int matrix_layout, char type, int64_t kl,
                           int64_t ku, double cfrom, double cto,
                           int64_t m, int64_t n, lapack_complex_double* a,
                           int64_t lda );

int64_t LAPACKE_slaset_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, float alpha, float beta, float* a,
                           int64_t lda );
int64_t LAPACKE_dlaset_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, double alpha, double beta, double* a,
                           int64_t lda );
int64_t LAPACKE_claset_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, lapack_complex_float alpha,
                           lapack_complex_float beta, lapack_complex_float* a,
                           int64_t lda );
int64_t LAPACKE_zlaset_64( int matrix_layout, char uplo, int64_t m,
                           int64_t n, lapack_complex_double alpha,
                           lapack_complex_double beta, lapack_complex_double* a,
                           int64_t lda );

int64_t LAPACKE_slasrt_64( char id, int64_t n, float* d );
int64_t LAPACKE_dlasrt_64( char id, int64_t n, double* d );

int64_t LAPACKE_slassq_64( int64_t n,                 float* x, int64_t incx,  float* scale,  float* sumsq );
int64_t LAPACKE_dlassq_64( int64_t n,                double* x, int64_t incx, double* scale, double* sumsq );
int64_t LAPACKE_classq_64( int64_t n,  lapack_complex_float* x, int64_t incx,  float* scale,  float* sumsq );
int64_t LAPACKE_zlassq_64( int64_t n, lapack_complex_double* x, int64_t incx, double* scale, double* sumsq );

int64_t LAPACKE_slaswp_64( int matrix_layout, int64_t n, float* a,
                           int64_t lda, int64_t k1, int64_t k2,
                           const int64_t* ipiv, int64_t incx );
int64_t LAPACKE_dlaswp_64( int matrix_layout, int64_t n, double* a,
                           int64_t lda, int64_t k1, int64_t k2,
                           const int64_t* ipiv, int64_t incx );
int64_t LAPACKE_claswp_64( int matrix_layout, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t k1, int64_t k2, const int64_t* ipiv,
                           int64_t incx );
int64_t LAPACKE_zlaswp_64( int matrix_layout, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t k1, int64_t k2, const int64_t* ipiv,
                           int64_t incx );

int64_t LAPACKE_slatms_64( int matrix_layout, int64_t m, int64_t n,
                           char dist, int64_t* iseed, char sym, float* d,
                           int64_t mode, float cond, float dmax,
                           int64_t kl, int64_t ku, char pack, float* a,
                           int64_t lda );
int64_t LAPACKE_dlatms_64( int matrix_layout, int64_t m, int64_t n,
                           char dist, int64_t* iseed, char sym, double* d,
                           int64_t mode, double cond, double dmax,
                           int64_t kl, int64_t ku, char pack, double* a,
                           int64_t lda );
int64_t LAPACKE_clatms_64( int matrix_layout, int64_t m, int64_t n,
                           char dist, int64_t* iseed, char sym, float* d,
                           int64_t mode, float cond, float dmax,
                           int64_t kl, int64_t ku, char pack,
                           lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zlatms_64( int matrix_layout, int64_t m, int64_t n,
                           char dist, int64_t* iseed, char sym, double* d,
                           int64_t mode, double cond, double dmax,
                           int64_t kl, int64_t ku, char pack,
                           lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_slauum_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda );
int64_t LAPACKE_dlauum_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda );
int64_t LAPACKE_clauum_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zlauum_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_sopgtr_64( int matrix_layout, char uplo, int64_t n,
                           const float* ap, const float* tau, float* q,
                           int64_t ldq );
int64_t LAPACKE_dopgtr_64( int matrix_layout, char uplo, int64_t n,
                           const double* ap, const double* tau, double* q,
                           int64_t ldq );

int64_t LAPACKE_sopmtr_64( int matrix_layout, char side, char uplo, char trans,
                           int64_t m, int64_t n, const float* ap,
                           const float* tau, float* c, int64_t ldc );
int64_t LAPACKE_dopmtr_64( int matrix_layout, char side, char uplo, char trans,
                           int64_t m, int64_t n, const double* ap,
                           const double* tau, double* c, int64_t ldc );

int64_t LAPACKE_sorgbr_64( int matrix_layout, char vect, int64_t m,
                           int64_t n, int64_t k, float* a, int64_t lda,
                           const float* tau );
int64_t LAPACKE_dorgbr_64( int matrix_layout, char vect, int64_t m,
                           int64_t n, int64_t k, double* a,
                           int64_t lda, const double* tau );

int64_t LAPACKE_sorghr_64( int matrix_layout, int64_t n, int64_t ilo,
                           int64_t ihi, float* a, int64_t lda,
                           const float* tau );
int64_t LAPACKE_dorghr_64( int matrix_layout, int64_t n, int64_t ilo,
                           int64_t ihi, double* a, int64_t lda,
                           const double* tau );

int64_t LAPACKE_sorglq_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, float* a, int64_t lda,
                           const float* tau );
int64_t LAPACKE_dorglq_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, double* a, int64_t lda,
                           const double* tau );

int64_t LAPACKE_sorgql_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, float* a, int64_t lda,
                           const float* tau );
int64_t LAPACKE_dorgql_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, double* a, int64_t lda,
                           const double* tau );

int64_t LAPACKE_sorgqr_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, float* a, int64_t lda,
                           const float* tau );
int64_t LAPACKE_dorgqr_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, double* a, int64_t lda,
                           const double* tau );

int64_t LAPACKE_sorgrq_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, float* a, int64_t lda,
                           const float* tau );
int64_t LAPACKE_dorgrq_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, double* a, int64_t lda,
                           const double* tau );

int64_t LAPACKE_sorgtr_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, const float* tau );
int64_t LAPACKE_dorgtr_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, const double* tau );

int64_t LAPACKE_sorgtsqr_row_64( int matrix_layout, int64_t m, int64_t n,
                                 int64_t mb, int64_t nb,
                                 float* a, int64_t lda,
                                 const float* t, int64_t ldt );
int64_t LAPACKE_dorgtsqr_row_64( int matrix_layout, int64_t m, int64_t n,
                                 int64_t mb, int64_t nb,
                                 double* a, int64_t lda,
                                 const double* t, int64_t ldt );

int64_t LAPACKE_sormbr_64( int matrix_layout, char vect, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const float* a, int64_t lda, const float* tau,
                           float* c, int64_t ldc );
int64_t LAPACKE_dormbr_64( int matrix_layout, char vect, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const double* a, int64_t lda, const double* tau,
                           double* c, int64_t ldc );

int64_t LAPACKE_sormhr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t ilo,
                           int64_t ihi, const float* a, int64_t lda,
                           const float* tau, float* c, int64_t ldc );
int64_t LAPACKE_dormhr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t ilo,
                           int64_t ihi, const double* a, int64_t lda,
                           const double* tau, double* c, int64_t ldc );

int64_t LAPACKE_sormlq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const float* a, int64_t lda, const float* tau,
                           float* c, int64_t ldc );
int64_t LAPACKE_dormlq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const double* a, int64_t lda, const double* tau,
                           double* c, int64_t ldc );

int64_t LAPACKE_sormql_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const float* a, int64_t lda, const float* tau,
                           float* c, int64_t ldc );
int64_t LAPACKE_dormql_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const double* a, int64_t lda, const double* tau,
                           double* c, int64_t ldc );

int64_t LAPACKE_sormqr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const float* a, int64_t lda, const float* tau,
                           float* c, int64_t ldc );
int64_t LAPACKE_dormqr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const double* a, int64_t lda, const double* tau,
                           double* c, int64_t ldc );

int64_t LAPACKE_sormrq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const float* a, int64_t lda, const float* tau,
                           float* c, int64_t ldc );
int64_t LAPACKE_dormrq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const double* a, int64_t lda, const double* tau,
                           double* c, int64_t ldc );

int64_t LAPACKE_sormrz_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           int64_t l, const float* a, int64_t lda,
                           const float* tau, float* c, int64_t ldc );
int64_t LAPACKE_dormrz_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           int64_t l, const double* a, int64_t lda,
                           const double* tau, double* c, int64_t ldc );

int64_t LAPACKE_sormtr_64( int matrix_layout, char side, char uplo, char trans,
                           int64_t m, int64_t n, const float* a,
                           int64_t lda, const float* tau, float* c,
                           int64_t ldc );
int64_t LAPACKE_dormtr_64( int matrix_layout, char side, char uplo, char trans,
                           int64_t m, int64_t n, const double* a,
                           int64_t lda, const double* tau, double* c,
                           int64_t ldc );

int64_t LAPACKE_spbcon_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, const float* ab, int64_t ldab,
                           float anorm, float* rcond );
int64_t LAPACKE_dpbcon_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, const double* ab, int64_t ldab,
                           double anorm, double* rcond );
int64_t LAPACKE_cpbcon_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, const lapack_complex_float* ab,
                           int64_t ldab, float anorm, float* rcond );
int64_t LAPACKE_zpbcon_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, const lapack_complex_double* ab,
                           int64_t ldab, double anorm, double* rcond );

int64_t LAPACKE_spbequ_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, const float* ab, int64_t ldab,
                           float* s, float* scond, float* amax );
int64_t LAPACKE_dpbequ_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, const double* ab, int64_t ldab,
                           double* s, double* scond, double* amax );
int64_t LAPACKE_cpbequ_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, const lapack_complex_float* ab,
                           int64_t ldab, float* s, float* scond,
                           float* amax );
int64_t LAPACKE_zpbequ_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, const lapack_complex_double* ab,
                           int64_t ldab, double* s, double* scond,
                           double* amax );

int64_t LAPACKE_spbrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs, const float* ab,
                           int64_t ldab, const float* afb, int64_t ldafb,
                           const float* b, int64_t ldb, float* x,
                           int64_t ldx, float* ferr, float* berr );
int64_t LAPACKE_dpbrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs, const double* ab,
                           int64_t ldab, const double* afb, int64_t ldafb,
                           const double* b, int64_t ldb, double* x,
                           int64_t ldx, double* ferr, double* berr );
int64_t LAPACKE_cpbrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs,
                           const lapack_complex_float* ab, int64_t ldab,
                           const lapack_complex_float* afb, int64_t ldafb,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zpbrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs,
                           const lapack_complex_double* ab, int64_t ldab,
                           const lapack_complex_double* afb, int64_t ldafb,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_spbstf_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kb, float* bb, int64_t ldbb );
int64_t LAPACKE_dpbstf_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kb, double* bb, int64_t ldbb );
int64_t LAPACKE_cpbstf_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kb, lapack_complex_float* bb,
                           int64_t ldbb );
int64_t LAPACKE_zpbstf_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kb, lapack_complex_double* bb,
                           int64_t ldbb );

int64_t LAPACKE_spbsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t kd, int64_t nrhs, float* ab,
                          int64_t ldab, float* b, int64_t ldb );
int64_t LAPACKE_dpbsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t kd, int64_t nrhs, double* ab,
                          int64_t ldab, double* b, int64_t ldb );
int64_t LAPACKE_cpbsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t kd, int64_t nrhs,
                          lapack_complex_float* ab, int64_t ldab,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpbsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t kd, int64_t nrhs,
                          lapack_complex_double* ab, int64_t ldab,
                          lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_spbsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs, float* ab,
                           int64_t ldab, float* afb, int64_t ldafb,
                           char* equed, float* s, float* b, int64_t ldb,
                           float* x, int64_t ldx, float* rcond, float* ferr,
                           float* berr );
int64_t LAPACKE_dpbsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs, double* ab,
                           int64_t ldab, double* afb, int64_t ldafb,
                           char* equed, double* s, double* b, int64_t ldb,
                           double* x, int64_t ldx, double* rcond,
                           double* ferr, double* berr );
int64_t LAPACKE_cpbsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs,
                           lapack_complex_float* ab, int64_t ldab,
                           lapack_complex_float* afb, int64_t ldafb,
                           char* equed, float* s, lapack_complex_float* b,
                           int64_t ldb, lapack_complex_float* x,
                           int64_t ldx, float* rcond, float* ferr,
                           float* berr );
int64_t LAPACKE_zpbsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs,
                           lapack_complex_double* ab, int64_t ldab,
                           lapack_complex_double* afb, int64_t ldafb,
                           char* equed, double* s, lapack_complex_double* b,
                           int64_t ldb, lapack_complex_double* x,
                           int64_t ldx, double* rcond, double* ferr,
                           double* berr );

int64_t LAPACKE_spbtrf_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, float* ab, int64_t ldab );
int64_t LAPACKE_dpbtrf_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, double* ab, int64_t ldab );
int64_t LAPACKE_cpbtrf_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, lapack_complex_float* ab,
                           int64_t ldab );
int64_t LAPACKE_zpbtrf_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, lapack_complex_double* ab,
                           int64_t ldab );

int64_t LAPACKE_spbtrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs, const float* ab,
                           int64_t ldab, float* b, int64_t ldb );
int64_t LAPACKE_dpbtrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs, const double* ab,
                           int64_t ldab, double* b, int64_t ldb );
int64_t LAPACKE_cpbtrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs,
                           const lapack_complex_float* ab, int64_t ldab,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpbtrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t kd, int64_t nrhs,
                           const lapack_complex_double* ab, int64_t ldab,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_spftrf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, float* a );
int64_t LAPACKE_dpftrf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, double* a );
int64_t LAPACKE_cpftrf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, lapack_complex_float* a );
int64_t LAPACKE_zpftrf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, lapack_complex_double* a );

int64_t LAPACKE_spftri_64( int matrix_layout, char transr, char uplo,
                           int64_t n, float* a );
int64_t LAPACKE_dpftri_64( int matrix_layout, char transr, char uplo,
                           int64_t n, double* a );
int64_t LAPACKE_cpftri_64( int matrix_layout, char transr, char uplo,
                           int64_t n, lapack_complex_float* a );
int64_t LAPACKE_zpftri_64( int matrix_layout, char transr, char uplo,
                           int64_t n, lapack_complex_double* a );

int64_t LAPACKE_spftrs_64( int matrix_layout, char transr, char uplo,
                           int64_t n, int64_t nrhs, const float* a,
                           float* b, int64_t ldb );
int64_t LAPACKE_dpftrs_64( int matrix_layout, char transr, char uplo,
                           int64_t n, int64_t nrhs, const double* a,
                           double* b, int64_t ldb );
int64_t LAPACKE_cpftrs_64( int matrix_layout, char transr, char uplo,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_float* a,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpftrs_64( int matrix_layout, char transr, char uplo,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_double* a,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_spocon_64( int matrix_layout, char uplo, int64_t n,
                           const float* a, int64_t lda, float anorm,
                           float* rcond );
int64_t LAPACKE_dpocon_64( int matrix_layout, char uplo, int64_t n,
                           const double* a, int64_t lda, double anorm,
                           double* rcond );
int64_t LAPACKE_cpocon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           float anorm, float* rcond );
int64_t LAPACKE_zpocon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           double anorm, double* rcond );

int64_t LAPACKE_spoequ_64( int matrix_layout, int64_t n, const float* a,
                           int64_t lda, float* s, float* scond,
                           float* amax );
int64_t LAPACKE_dpoequ_64( int matrix_layout, int64_t n, const double* a,
                           int64_t lda, double* s, double* scond,
                           double* amax );
int64_t LAPACKE_cpoequ_64( int matrix_layout, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           float* s, float* scond, float* amax );
int64_t LAPACKE_zpoequ_64( int matrix_layout, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           double* s, double* scond, double* amax );

int64_t LAPACKE_spoequb_64( int matrix_layout, int64_t n, const float* a,
                            int64_t lda, float* s, float* scond,
                            float* amax );
int64_t LAPACKE_dpoequb_64( int matrix_layout, int64_t n, const double* a,
                            int64_t lda, double* s, double* scond,
                            double* amax );
int64_t LAPACKE_cpoequb_64( int matrix_layout, int64_t n,
                            const lapack_complex_float* a, int64_t lda,
                            float* s, float* scond, float* amax );
int64_t LAPACKE_zpoequb_64( int matrix_layout, int64_t n,
                            const lapack_complex_double* a, int64_t lda,
                            double* s, double* scond, double* amax );

int64_t LAPACKE_sporfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const float* af, int64_t ldaf, const float* b,
                           int64_t ldb, float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_dporfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           const double* af, int64_t ldaf, const double* b,
                           int64_t ldb, double* x, int64_t ldx,
                           double* ferr, double* berr );
int64_t LAPACKE_cporfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* af,
                           int64_t ldaf, const lapack_complex_float* b,
                           int64_t ldb, lapack_complex_float* x,
                           int64_t ldx, float* ferr, float* berr );
int64_t LAPACKE_zporfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* af,
                           int64_t ldaf, const lapack_complex_double* b,
                           int64_t ldb, lapack_complex_double* x,
                           int64_t ldx, double* ferr, double* berr );

int64_t LAPACKE_sporfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs, const float* a,
                            int64_t lda, const float* af, int64_t ldaf,
                            const float* s, const float* b, int64_t ldb,
                            float* x, int64_t ldx, float* rcond, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_dporfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs, const double* a,
                            int64_t lda, const double* af, int64_t ldaf,
                            const double* s, const double* b, int64_t ldb,
                            double* x, int64_t ldx, double* rcond,
                            double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );
int64_t LAPACKE_cporfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs,
                            const lapack_complex_float* a, int64_t lda,
                            const lapack_complex_float* af, int64_t ldaf,
                            const float* s, const lapack_complex_float* b,
                            int64_t ldb, lapack_complex_float* x,
                            int64_t ldx, float* rcond, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_zporfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs,
                            const lapack_complex_double* a, int64_t lda,
                            const lapack_complex_double* af, int64_t ldaf,
                            const double* s, const lapack_complex_double* b,
                            int64_t ldb, lapack_complex_double* x,
                            int64_t ldx, double* rcond, double* berr,
                            int64_t n_err_bnds, double* err_bnds_norm,
                            double* err_bnds_comp, int64_t nparams,
                            double* params );

int64_t LAPACKE_sposv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, float* a, int64_t lda, float* b,
                          int64_t ldb );
int64_t LAPACKE_dposv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, double* a, int64_t lda, double* b,
                          int64_t ldb );
int64_t LAPACKE_cposv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* b,
                          int64_t ldb );
int64_t LAPACKE_zposv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* b,
                          int64_t ldb );
int64_t LAPACKE_dsposv_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, double* a, int64_t lda,
                           double* b, int64_t ldb, double* x, int64_t ldx,
                           int64_t* iter );
int64_t LAPACKE_zcposv_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb, lapack_complex_double* x,
                           int64_t ldx, int64_t* iter );

int64_t LAPACKE_sposvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, float* a, int64_t lda, float* af,
                           int64_t ldaf, char* equed, float* s, float* b,
                           int64_t ldb, float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr );
int64_t LAPACKE_dposvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, double* a, int64_t lda,
                           double* af, int64_t ldaf, char* equed, double* s,
                           double* b, int64_t ldb, double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );
int64_t LAPACKE_cposvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* af,
                           int64_t ldaf, char* equed, float* s,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr );
int64_t LAPACKE_zposvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* af,
                           int64_t ldaf, char* equed, double* s,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );

int64_t LAPACKE_sposvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs, float* a,
                            int64_t lda, float* af, int64_t ldaf,
                            char* equed, float* s, float* b, int64_t ldb,
                            float* x, int64_t ldx, float* rcond,
                            float* rpvgrw, float* berr, int64_t n_err_bnds,
                            float* err_bnds_norm, float* err_bnds_comp,
                            int64_t nparams, float* params );
int64_t LAPACKE_dposvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs, double* a,
                            int64_t lda, double* af, int64_t ldaf,
                            char* equed, double* s, double* b, int64_t ldb,
                            double* x, int64_t ldx, double* rcond,
                            double* rpvgrw, double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );
int64_t LAPACKE_cposvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* af, int64_t ldaf,
                            char* equed, float* s, lapack_complex_float* b,
                            int64_t ldb, lapack_complex_float* x,
                            int64_t ldx, float* rcond, float* rpvgrw,
                            float* berr, int64_t n_err_bnds,
                            float* err_bnds_norm, float* err_bnds_comp,
                            int64_t nparams, float* params );
int64_t LAPACKE_zposvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* af, int64_t ldaf,
                            char* equed, double* s, lapack_complex_double* b,
                            int64_t ldb, lapack_complex_double* x,
                            int64_t ldx, double* rcond, double* rpvgrw,
                            double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );

int64_t LAPACKE_spotrf2_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda );
int64_t LAPACKE_dpotrf2_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda );
int64_t LAPACKE_cpotrf2_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zpotrf2_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_spotrf_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda );
int64_t LAPACKE_dpotrf_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda );
int64_t LAPACKE_cpotrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zpotrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_spotri_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda );
int64_t LAPACKE_dpotri_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda );
int64_t LAPACKE_cpotri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zpotri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_spotrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           float* b, int64_t ldb );
int64_t LAPACKE_dpotrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           double* b, int64_t ldb );
int64_t LAPACKE_cpotrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb );
int64_t LAPACKE_zpotrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb );

int64_t LAPACKE_sppcon_64( int matrix_layout, char uplo, int64_t n,
                           const float* ap, float anorm, float* rcond );
int64_t LAPACKE_dppcon_64( int matrix_layout, char uplo, int64_t n,
                           const double* ap, double anorm, double* rcond );
int64_t LAPACKE_cppcon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* ap, float anorm,
                           float* rcond );
int64_t LAPACKE_zppcon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* ap, double anorm,
                           double* rcond );

int64_t LAPACKE_sppequ_64( int matrix_layout, char uplo, int64_t n,
                           const float* ap, float* s, float* scond,
                           float* amax );
int64_t LAPACKE_dppequ_64( int matrix_layout, char uplo, int64_t n,
                           const double* ap, double* s, double* scond,
                           double* amax );
int64_t LAPACKE_cppequ_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* ap, float* s,
                           float* scond, float* amax );
int64_t LAPACKE_zppequ_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* ap, double* s,
                           double* scond, double* amax );

int64_t LAPACKE_spprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* ap, const float* afp,
                           const float* b, int64_t ldb, float* x,
                           int64_t ldx, float* ferr, float* berr );
int64_t LAPACKE_dpprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* ap, const double* afp,
                           const double* b, int64_t ldb, double* x,
                           int64_t ldx, double* ferr, double* berr );
int64_t LAPACKE_cpprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* ap,
                           const lapack_complex_float* afp,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zpprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* ap,
                           const lapack_complex_double* afp,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_sppsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, float* ap, float* b,
                          int64_t ldb );
int64_t LAPACKE_dppsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, double* ap, double* b,
                          int64_t ldb );
int64_t LAPACKE_cppsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* ap,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zppsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* ap,
                          lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sppsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, float* ap, float* afp, char* equed,
                           float* s, float* b, int64_t ldb, float* x,
                           int64_t ldx, float* rcond, float* ferr,
                           float* berr );
int64_t LAPACKE_dppsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, double* ap, double* afp,
                           char* equed, double* s, double* b, int64_t ldb,
                           double* x, int64_t ldx, double* rcond,
                           double* ferr, double* berr );
int64_t LAPACKE_cppsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, lapack_complex_float* ap,
                           lapack_complex_float* afp, char* equed, float* s,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr );
int64_t LAPACKE_zppsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, lapack_complex_double* ap,
                           lapack_complex_double* afp, char* equed, double* s,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );

int64_t LAPACKE_spptrf_64( int matrix_layout, char uplo, int64_t n,
                           float* ap );
int64_t LAPACKE_dpptrf_64( int matrix_layout, char uplo, int64_t n,
                           double* ap );
int64_t LAPACKE_cpptrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* ap );
int64_t LAPACKE_zpptrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* ap );

int64_t LAPACKE_spptri_64( int matrix_layout, char uplo, int64_t n,
                           float* ap );
int64_t LAPACKE_dpptri_64( int matrix_layout, char uplo, int64_t n,
                           double* ap );
int64_t LAPACKE_cpptri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* ap );
int64_t LAPACKE_zpptri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* ap );

int64_t LAPACKE_spptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* ap, float* b,
                           int64_t ldb );
int64_t LAPACKE_dpptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* ap, double* b,
                           int64_t ldb );
int64_t LAPACKE_cpptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* ap,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* ap,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_spstrf_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, int64_t* piv, int64_t* rank,
                           float tol );
int64_t LAPACKE_dpstrf_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, int64_t* piv, int64_t* rank,
                           double tol );
int64_t LAPACKE_cpstrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* piv, int64_t* rank, float tol );
int64_t LAPACKE_zpstrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* piv, int64_t* rank, double tol );

int64_t LAPACKE_sptcon_64( int64_t n, const float* d, const float* e,
                           float anorm, float* rcond );
int64_t LAPACKE_dptcon_64( int64_t n, const double* d, const double* e,
                           double anorm, double* rcond );
int64_t LAPACKE_cptcon_64( int64_t n, const float* d,
                           const lapack_complex_float* e, float anorm,
                           float* rcond );
int64_t LAPACKE_zptcon_64( int64_t n, const double* d,
                           const lapack_complex_double* e, double anorm,
                           double* rcond );

int64_t LAPACKE_spteqr_64( int matrix_layout, char compz, int64_t n, float* d,
                           float* e, float* z, int64_t ldz );
int64_t LAPACKE_dpteqr_64( int matrix_layout, char compz, int64_t n,
                           double* d, double* e, double* z, int64_t ldz );
int64_t LAPACKE_cpteqr_64( int matrix_layout, char compz, int64_t n, float* d,
                           float* e, lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zpteqr_64( int matrix_layout, char compz, int64_t n,
                           double* d, double* e, lapack_complex_double* z,
                           int64_t ldz );

int64_t LAPACKE_sptrfs_64( int matrix_layout, int64_t n, int64_t nrhs,
                           const float* d, const float* e, const float* df,
                           const float* ef, const float* b, int64_t ldb,
                           float* x, int64_t ldx, float* ferr, float* berr );
int64_t LAPACKE_dptrfs_64( int matrix_layout, int64_t n, int64_t nrhs,
                           const double* d, const double* e, const double* df,
                           const double* ef, const double* b, int64_t ldb,
                           double* x, int64_t ldx, double* ferr,
                           double* berr );
int64_t LAPACKE_cptrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* d,
                           const lapack_complex_float* e, const float* df,
                           const lapack_complex_float* ef,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zptrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* d,
                           const lapack_complex_double* e, const double* df,
                           const lapack_complex_double* ef,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_sptsv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          float* d, float* e, float* b, int64_t ldb );
int64_t LAPACKE_dptsv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          double* d, double* e, double* b, int64_t ldb );
int64_t LAPACKE_cptsv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          float* d, lapack_complex_float* e,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zptsv_64( int matrix_layout, int64_t n, int64_t nrhs,
                          double* d, lapack_complex_double* e,
                          lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sptsvx_64( int matrix_layout, char fact, int64_t n,
                           int64_t nrhs, const float* d, const float* e,
                           float* df, float* ef, const float* b, int64_t ldb,
                           float* x, int64_t ldx, float* rcond, float* ferr,
                           float* berr );
int64_t LAPACKE_dptsvx_64( int matrix_layout, char fact, int64_t n,
                           int64_t nrhs, const double* d, const double* e,
                           double* df, double* ef, const double* b,
                           int64_t ldb, double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );
int64_t LAPACKE_cptsvx_64( int matrix_layout, char fact, int64_t n,
                           int64_t nrhs, const float* d,
                           const lapack_complex_float* e, float* df,
                           lapack_complex_float* ef,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr );
int64_t LAPACKE_zptsvx_64( int matrix_layout, char fact, int64_t n,
                           int64_t nrhs, const double* d,
                           const lapack_complex_double* e, double* df,
                           lapack_complex_double* ef,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );

int64_t LAPACKE_spttrf_64( int64_t n, float* d, float* e );
int64_t LAPACKE_dpttrf_64( int64_t n, double* d, double* e );
int64_t LAPACKE_cpttrf_64( int64_t n, float* d, lapack_complex_float* e );
int64_t LAPACKE_zpttrf_64( int64_t n, double* d, lapack_complex_double* e );

int64_t LAPACKE_spttrs_64( int matrix_layout, int64_t n, int64_t nrhs,
                           const float* d, const float* e, float* b,
                           int64_t ldb );
int64_t LAPACKE_dpttrs_64( int matrix_layout, int64_t n, int64_t nrhs,
                           const double* d, const double* e, double* b,
                           int64_t ldb );
int64_t LAPACKE_cpttrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* d,
                           const lapack_complex_float* e,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpttrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* d,
                           const lapack_complex_double* e,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_ssbev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t kd, float* ab, int64_t ldab, float* w,
                          float* z, int64_t ldz );
int64_t LAPACKE_dsbev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t kd, double* ab, int64_t ldab, double* w,
                          double* z, int64_t ldz );

int64_t LAPACKE_ssbevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t kd, float* ab, int64_t ldab, float* w,
                           float* z, int64_t ldz );
int64_t LAPACKE_dsbevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t kd, double* ab, int64_t ldab,
                           double* w, double* z, int64_t ldz );

int64_t LAPACKE_ssbevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t kd, float* ab,
                           int64_t ldab, float* q, int64_t ldq, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_dsbevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t kd, double* ab,
                           int64_t ldab, double* q, int64_t ldq,
                           double vl, double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* ifail );

int64_t LAPACKE_ssbgst_64( int matrix_layout, char vect, char uplo, int64_t n,
                           int64_t ka, int64_t kb, float* ab,
                           int64_t ldab, const float* bb, int64_t ldbb,
                           float* x, int64_t ldx );
int64_t LAPACKE_dsbgst_64( int matrix_layout, char vect, char uplo, int64_t n,
                           int64_t ka, int64_t kb, double* ab,
                           int64_t ldab, const double* bb, int64_t ldbb,
                           double* x, int64_t ldx );

int64_t LAPACKE_ssbgv_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t ka, int64_t kb, float* ab,
                          int64_t ldab, float* bb, int64_t ldbb, float* w,
                          float* z, int64_t ldz );
int64_t LAPACKE_dsbgv_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t ka, int64_t kb, double* ab,
                          int64_t ldab, double* bb, int64_t ldbb,
                          double* w, double* z, int64_t ldz );

int64_t LAPACKE_ssbgvd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t ka, int64_t kb, float* ab,
                           int64_t ldab, float* bb, int64_t ldbb,
                           float* w, float* z, int64_t ldz );
int64_t LAPACKE_dsbgvd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t ka, int64_t kb, double* ab,
                           int64_t ldab, double* bb, int64_t ldbb,
                           double* w, double* z, int64_t ldz );

int64_t LAPACKE_ssbgvx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t ka, int64_t kb,
                           float* ab, int64_t ldab, float* bb,
                           int64_t ldbb, float* q, int64_t ldq, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_dsbgvx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t ka, int64_t kb,
                           double* ab, int64_t ldab, double* bb,
                           int64_t ldbb, double* q, int64_t ldq,
                           double vl, double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* ifail );

int64_t LAPACKE_ssbtrd_64( int matrix_layout, char vect, char uplo, int64_t n,
                           int64_t kd, float* ab, int64_t ldab, float* d,
                           float* e, float* q, int64_t ldq );
int64_t LAPACKE_dsbtrd_64( int matrix_layout, char vect, char uplo, int64_t n,
                           int64_t kd, double* ab, int64_t ldab,
                           double* d, double* e, double* q, int64_t ldq );

int64_t LAPACKE_ssfrk_64( int matrix_layout, char transr, char uplo, char trans,
                          int64_t n, int64_t k, float alpha,
                          const float* a, int64_t lda, float beta,
                          float* c );
int64_t LAPACKE_dsfrk_64( int matrix_layout, char transr, char uplo, char trans,
                          int64_t n, int64_t k, double alpha,
                          const double* a, int64_t lda, double beta,
                          double* c );

int64_t LAPACKE_sspcon_64( int matrix_layout, char uplo, int64_t n,
                           const float* ap, const int64_t* ipiv, float anorm,
                           float* rcond );
int64_t LAPACKE_dspcon_64( int matrix_layout, char uplo, int64_t n,
                           const double* ap, const int64_t* ipiv,
                           double anorm, double* rcond );
int64_t LAPACKE_cspcon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* ap,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_zspcon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* ap,
                           const int64_t* ipiv, double anorm,
                           double* rcond );

int64_t LAPACKE_sspev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          float* ap, float* w, float* z, int64_t ldz );
int64_t LAPACKE_dspev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          double* ap, double* w, double* z, int64_t ldz );

int64_t LAPACKE_sspevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           float* ap, float* w, float* z, int64_t ldz );
int64_t LAPACKE_dspevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           double* ap, double* w, double* z, int64_t ldz );

int64_t LAPACKE_sspevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, float* ap, float vl, float vu,
                           int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_dspevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, double* ap, double vl, double vu,
                           int64_t il, int64_t iu, double abstol,
                           int64_t* m, double* w, double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_sspgst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, float* ap, const float* bp );
int64_t LAPACKE_dspgst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, double* ap, const double* bp );

int64_t LAPACKE_sspgv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, float* ap, float* bp,
                          float* w, float* z, int64_t ldz );
int64_t LAPACKE_dspgv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, double* ap, double* bp,
                          double* w, double* z, int64_t ldz );

int64_t LAPACKE_sspgvd_64( int matrix_layout, int64_t itype, char jobz,
                           char uplo, int64_t n, float* ap, float* bp,
                           float* w, float* z, int64_t ldz );
int64_t LAPACKE_dspgvd_64( int matrix_layout, int64_t itype, char jobz,
                           char uplo, int64_t n, double* ap, double* bp,
                           double* w, double* z, int64_t ldz );

int64_t LAPACKE_sspgvx_64( int matrix_layout, int64_t itype, char jobz,
                           char range, char uplo, int64_t n, float* ap,
                           float* bp, float vl, float vu, int64_t il,
                           int64_t iu, float abstol, int64_t* m, float* w,
                           float* z, int64_t ldz, int64_t* ifail );
int64_t LAPACKE_dspgvx_64( int matrix_layout, int64_t itype, char jobz,
                           char range, char uplo, int64_t n, double* ap,
                           double* bp, double vl, double vu, int64_t il,
                           int64_t iu, double abstol, int64_t* m,
                           double* w, double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_ssprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* ap, const float* afp,
                           const int64_t* ipiv, const float* b,
                           int64_t ldb, float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_dsprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* ap, const double* afp,
                           const int64_t* ipiv, const double* b,
                           int64_t ldb, double* x, int64_t ldx,
                           double* ferr, double* berr );
int64_t LAPACKE_csprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* ap,
                           const lapack_complex_float* afp,
                           const int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zsprfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* ap,
                           const lapack_complex_double* afp,
                           const int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_sspsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, float* ap, int64_t* ipiv,
                          float* b, int64_t ldb );
int64_t LAPACKE_dspsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, double* ap, int64_t* ipiv,
                          double* b, int64_t ldb );
int64_t LAPACKE_cspsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* ap,
                          int64_t* ipiv, lapack_complex_float* b,
                          int64_t ldb );
int64_t LAPACKE_zspsv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* ap,
                          int64_t* ipiv, lapack_complex_double* b,
                          int64_t ldb );

int64_t LAPACKE_sspsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const float* ap, float* afp,
                           int64_t* ipiv, const float* b, int64_t ldb,
                           float* x, int64_t ldx, float* rcond, float* ferr,
                           float* berr );
int64_t LAPACKE_dspsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const double* ap, double* afp,
                           int64_t* ipiv, const double* b, int64_t ldb,
                           double* x, int64_t ldx, double* rcond,
                           double* ferr, double* berr );
int64_t LAPACKE_cspsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* ap,
                           lapack_complex_float* afp, int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr );
int64_t LAPACKE_zspsvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* ap,
                           lapack_complex_double* afp, int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );

int64_t LAPACKE_ssptrd_64( int matrix_layout, char uplo, int64_t n, float* ap,
                           float* d, float* e, float* tau );
int64_t LAPACKE_dsptrd_64( int matrix_layout, char uplo, int64_t n,
                           double* ap, double* d, double* e, double* tau );

int64_t LAPACKE_ssptrf_64( int matrix_layout, char uplo, int64_t n, float* ap,
                           int64_t* ipiv );
int64_t LAPACKE_dsptrf_64( int matrix_layout, char uplo, int64_t n,
                           double* ap, int64_t* ipiv );
int64_t LAPACKE_csptrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* ap, int64_t* ipiv );
int64_t LAPACKE_zsptrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* ap, int64_t* ipiv );

int64_t LAPACKE_ssptri_64( int matrix_layout, char uplo, int64_t n, float* ap,
                           const int64_t* ipiv );
int64_t LAPACKE_dsptri_64( int matrix_layout, char uplo, int64_t n,
                           double* ap, const int64_t* ipiv );
int64_t LAPACKE_csptri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* ap, const int64_t* ipiv );
int64_t LAPACKE_zsptri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* ap, const int64_t* ipiv );

int64_t LAPACKE_ssptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* ap,
                           const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dsptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* ap,
                           const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_csptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* ap,
                           const int64_t* ipiv, lapack_complex_float* b,
                           int64_t ldb );
int64_t LAPACKE_zsptrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* ap,
                           const int64_t* ipiv, lapack_complex_double* b,
                           int64_t ldb );

int64_t LAPACKE_sstebz_64( char range, char order, int64_t n, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           const float* d, const float* e, int64_t* m,
                           int64_t* nsplit, float* w, int64_t* iblock,
                           int64_t* isplit );
int64_t LAPACKE_dstebz_64( char range, char order, int64_t n, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, const double* d, const double* e,
                           int64_t* m, int64_t* nsplit, double* w,
                           int64_t* iblock, int64_t* isplit );

int64_t LAPACKE_sstedc_64( int matrix_layout, char compz, int64_t n, float* d,
                           float* e, float* z, int64_t ldz );
int64_t LAPACKE_dstedc_64( int matrix_layout, char compz, int64_t n,
                           double* d, double* e, double* z, int64_t ldz );
int64_t LAPACKE_cstedc_64( int matrix_layout, char compz, int64_t n, float* d,
                           float* e, lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zstedc_64( int matrix_layout, char compz, int64_t n,
                           double* d, double* e, lapack_complex_double* z,
                           int64_t ldz );

int64_t LAPACKE_sstegr_64( int matrix_layout, char jobz, char range,
                           int64_t n, float* d, float* e, float vl, float vu,
                           int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* isuppz );
int64_t LAPACKE_dstegr_64( int matrix_layout, char jobz, char range,
                           int64_t n, double* d, double* e, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* isuppz );
int64_t LAPACKE_cstegr_64( int matrix_layout, char jobz, char range,
                           int64_t n, float* d, float* e, float vl, float vu,
                           int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, lapack_complex_float* z,
                           int64_t ldz, int64_t* isuppz );
int64_t LAPACKE_zstegr_64( int matrix_layout, char jobz, char range,
                           int64_t n, double* d, double* e, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t* isuppz );

int64_t LAPACKE_sstein_64( int matrix_layout, int64_t n, const float* d,
                           const float* e, int64_t m, const float* w,
                           const int64_t* iblock, const int64_t* isplit,
                           float* z, int64_t ldz, int64_t* ifailv );
int64_t LAPACKE_dstein_64( int matrix_layout, int64_t n, const double* d,
                           const double* e, int64_t m, const double* w,
                           const int64_t* iblock, const int64_t* isplit,
                           double* z, int64_t ldz, int64_t* ifailv );
int64_t LAPACKE_cstein_64( int matrix_layout, int64_t n, const float* d,
                           const float* e, int64_t m, const float* w,
                           const int64_t* iblock, const int64_t* isplit,
                           lapack_complex_float* z, int64_t ldz,
                           int64_t* ifailv );
int64_t LAPACKE_zstein_64( int matrix_layout, int64_t n, const double* d,
                           const double* e, int64_t m, const double* w,
                           const int64_t* iblock, const int64_t* isplit,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t* ifailv );

int64_t LAPACKE_sstemr_64( int matrix_layout, char jobz, char range,
                           int64_t n, float* d, float* e, float vl, float vu,
                           int64_t il, int64_t iu, int64_t* m,
                           float* w, float* z, int64_t ldz, int64_t nzc,
                           int64_t* isuppz, lapack_logical* tryrac );
int64_t LAPACKE_dstemr_64( int matrix_layout, char jobz, char range,
                           int64_t n, double* d, double* e, double vl,
                           double vu, int64_t il, int64_t iu,
                           int64_t* m, double* w, double* z, int64_t ldz,
                           int64_t nzc, int64_t* isuppz,
                           lapack_logical* tryrac );
int64_t LAPACKE_cstemr_64( int matrix_layout, char jobz, char range,
                           int64_t n, float* d, float* e, float vl, float vu,
                           int64_t il, int64_t iu, int64_t* m,
                           float* w, lapack_complex_float* z, int64_t ldz,
                           int64_t nzc, int64_t* isuppz,
                           lapack_logical* tryrac );
int64_t LAPACKE_zstemr_64( int matrix_layout, char jobz, char range,
                           int64_t n, double* d, double* e, double vl,
                           double vu, int64_t il, int64_t iu,
                           int64_t* m, double* w, lapack_complex_double* z,
                           int64_t ldz, int64_t nzc, int64_t* isuppz,
                           lapack_logical* tryrac );

int64_t LAPACKE_ssteqr_64( int matrix_layout, char compz, int64_t n, float* d,
                           float* e, float* z, int64_t ldz );
int64_t LAPACKE_dsteqr_64( int matrix_layout, char compz, int64_t n,
                           double* d, double* e, double* z, int64_t ldz );
int64_t LAPACKE_csteqr_64( int matrix_layout, char compz, int64_t n, float* d,
                           float* e, lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zsteqr_64( int matrix_layout, char compz, int64_t n,
                           double* d, double* e, lapack_complex_double* z,
                           int64_t ldz );

int64_t LAPACKE_skteqr_64( int matrix_layout, char compz, int64_t n,
                           float* e, float* z, int64_t ldz );
int64_t LAPACKE_dkteqr_64( int matrix_layout, char compz, int64_t n,
                           double* e, double* z, int64_t ldz );

int64_t LAPACKE_ssterf_64( int64_t n, float* d, float* e );
int64_t LAPACKE_dsterf_64( int64_t n, double* d, double* e );

int64_t LAPACKE_sstev_64( int matrix_layout, char jobz, int64_t n, float* d,
                          float* e, float* z, int64_t ldz );
int64_t LAPACKE_dstev_64( int matrix_layout, char jobz, int64_t n, double* d,
                          double* e, double* z, int64_t ldz );

int64_t LAPACKE_sktev_64( int matrix_layout, char jobz, int64_t n, float* d,
                          float* e, float* z, int64_t ldz );
int64_t LAPACKE_dktev_64( int matrix_layout, char jobz, int64_t n, double* d,
                          double* e, double* z, int64_t ldz );

int64_t LAPACKE_sstevd_64( int matrix_layout, char jobz, int64_t n, float* d,
                           float* e, float* z, int64_t ldz );
int64_t LAPACKE_dstevd_64( int matrix_layout, char jobz, int64_t n, double* d,
                           double* e, double* z, int64_t ldz );

int64_t LAPACKE_sstevr_64( int matrix_layout, char jobz, char range,
                           int64_t n, float* d, float* e, float vl, float vu,
                           int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* isuppz );
int64_t LAPACKE_dstevr_64( int matrix_layout, char jobz, char range,
                           int64_t n, double* d, double* e, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* isuppz );

int64_t LAPACKE_sstevx_64( int matrix_layout, char jobz, char range,
                           int64_t n, float* d, float* e, float vl, float vu,
                           int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_dstevx_64( int matrix_layout, char jobz, char range,
                           int64_t n, double* d, double* e, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* ifail );

int64_t LAPACKE_ssycon_64( int matrix_layout, char uplo, int64_t n,
                           const float* a, int64_t lda,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_dsycon_64( int matrix_layout, char uplo, int64_t n,
                           const double* a, int64_t lda,
                           const int64_t* ipiv, double anorm,
                           double* rcond );
int64_t LAPACKE_csycon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_zsycon_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           const int64_t* ipiv, double anorm,
                           double* rcond );

int64_t LAPACKE_skycon_64( int matrix_layout, char uplo, int64_t n,
                           const float* a, int64_t lda,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_dkycon_64( int matrix_layout, char uplo, int64_t n,
                           const double* a, int64_t lda,
                           const int64_t* ipiv, double anorm,
                           double* rcond );

int64_t LAPACKE_ssyequb_64( int matrix_layout, char uplo, int64_t n,
                            const float* a, int64_t lda, float* s,
                            float* scond, float* amax );
int64_t LAPACKE_dsyequb_64( int matrix_layout, char uplo, int64_t n,
                            const double* a, int64_t lda, double* s,
                            double* scond, double* amax );
int64_t LAPACKE_csyequb_64( int matrix_layout, char uplo, int64_t n,
                            const lapack_complex_float* a, int64_t lda,
                            float* s, float* scond, float* amax );
int64_t LAPACKE_zsyequb_64( int matrix_layout, char uplo, int64_t n,
                            const lapack_complex_double* a, int64_t lda,
                            double* s, double* scond, double* amax );

int64_t LAPACKE_ssyev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          float* a, int64_t lda, float* w );
int64_t LAPACKE_dsyev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          double* a, int64_t lda, double* w );

int64_t LAPACKE_skyev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          float* a, int64_t lda, float* w );
int64_t LAPACKE_dkyev_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          double* a, int64_t lda, double* w );

int64_t LAPACKE_ssyevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           float* a, int64_t lda, float* w );
int64_t LAPACKE_dsyevd_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           double* a, int64_t lda, double* w );

int64_t LAPACKE_ssyevr_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, float* a, int64_t lda, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* isuppz );
int64_t LAPACKE_dsyevr_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, double* a, int64_t lda, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* isuppz );

int64_t LAPACKE_ssyevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, float* a, int64_t lda, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_dsyevx_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, double* a, int64_t lda, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* ifail );

int64_t LAPACKE_ssygst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, float* a, int64_t lda,
                           const float* b, int64_t ldb );
int64_t LAPACKE_dsygst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, double* a, int64_t lda,
                           const double* b, int64_t ldb );

int64_t LAPACKE_skygst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, float* a, int64_t lda,
                           const float* b, int64_t ldb );
int64_t LAPACKE_dkygst_64( int matrix_layout, int64_t itype, char uplo,
                           int64_t n, double* a, int64_t lda,
                           const double* b, int64_t ldb );

int64_t LAPACKE_ssygv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, float* a, int64_t lda,
                          float* b, int64_t ldb, float* w );
int64_t LAPACKE_dsygv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, double* a, int64_t lda,
                          double* b, int64_t ldb, double* w );

int64_t LAPACKE_skygv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, float* a, int64_t lda,
                          float* b, int64_t ldb, float* w );
int64_t LAPACKE_dkygv_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, double* a, int64_t lda,
                          double* b, int64_t ldb, double* w );

int64_t LAPACKE_ssygvd_64( int matrix_layout, int64_t itype, char jobz,
                           char uplo, int64_t n, float* a, int64_t lda,
                           float* b, int64_t ldb, float* w );
int64_t LAPACKE_dsygvd_64( int matrix_layout, int64_t itype, char jobz,
                           char uplo, int64_t n, double* a, int64_t lda,
                           double* b, int64_t ldb, double* w );

int64_t LAPACKE_ssygvx_64( int matrix_layout, int64_t itype, char jobz,
                           char range, char uplo, int64_t n, float* a,
                           int64_t lda, float* b, int64_t ldb, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_dsygvx_64( int matrix_layout, int64_t itype, char jobz,
                           char range, char uplo, int64_t n, double* a,
                           int64_t lda, double* b, int64_t ldb, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* ifail );

int64_t LAPACKE_ssyrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const float* af, int64_t ldaf,
                           const int64_t* ipiv, const float* b,
                           int64_t ldb, float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_dsyrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           const double* af, int64_t ldaf,
                           const int64_t* ipiv, const double* b,
                           int64_t ldb, double* x, int64_t ldx,
                           double* ferr, double* berr );
int64_t LAPACKE_csyrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* af,
                           int64_t ldaf, const int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_zsyrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* af,
                           int64_t ldaf, const int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_skyrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const float* af, int64_t ldaf,
                           const int64_t* ipiv, const float* b,
                           int64_t ldb, float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_dkyrfs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           const double* af, int64_t ldaf,
                           const int64_t* ipiv, const double* b,
                           int64_t ldb, double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_ssyrfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs, const float* a,
                            int64_t lda, const float* af, int64_t ldaf,
                            const int64_t* ipiv, const float* s,
                            const float* b, int64_t ldb, float* x,
                            int64_t ldx, float* rcond, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_dsyrfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs, const double* a,
                            int64_t lda, const double* af, int64_t ldaf,
                            const int64_t* ipiv, const double* s,
                            const double* b, int64_t ldb, double* x,
                            int64_t ldx, double* rcond, double* berr,
                            int64_t n_err_bnds, double* err_bnds_norm,
                            double* err_bnds_comp, int64_t nparams,
                            double* params );
int64_t LAPACKE_csyrfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs,
                            const lapack_complex_float* a, int64_t lda,
                            const lapack_complex_float* af, int64_t ldaf,
                            const int64_t* ipiv, const float* s,
                            const lapack_complex_float* b, int64_t ldb,
                            lapack_complex_float* x, int64_t ldx,
                            float* rcond, float* berr, int64_t n_err_bnds,
                            float* err_bnds_norm, float* err_bnds_comp,
                            int64_t nparams, float* params );
int64_t LAPACKE_zsyrfsx_64( int matrix_layout, char uplo, char equed,
                            int64_t n, int64_t nrhs,
                            const lapack_complex_double* a, int64_t lda,
                            const lapack_complex_double* af, int64_t ldaf,
                            const int64_t* ipiv, const double* s,
                            const lapack_complex_double* b, int64_t ldb,
                            lapack_complex_double* x, int64_t ldx,
                            double* rcond, double* berr, int64_t n_err_bnds,
                            double* err_bnds_norm, double* err_bnds_comp,
                            int64_t nparams, double* params );

int64_t LAPACKE_ssysv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, float* a, int64_t lda,
                          int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dsysv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, double* a, int64_t lda,
                          int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_csysv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, int64_t* ipiv,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zsysv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, int64_t* ipiv,
                          lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_skysv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, float* a, int64_t lda,
                          int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dkysv_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, double* a, int64_t lda,
                          int64_t* ipiv, double* b, int64_t ldb );

int64_t LAPACKE_ssysvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           float* af, int64_t ldaf, int64_t* ipiv,
                           const float* b, int64_t ldb, float* x,
                           int64_t ldx, float* rcond, float* ferr,
                           float* berr );
int64_t LAPACKE_dsysvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           double* af, int64_t ldaf, int64_t* ipiv,
                           const double* b, int64_t ldb, double* x,
                           int64_t ldx, double* rcond, double* ferr,
                           double* berr );
int64_t LAPACKE_csysvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* af,
                           int64_t ldaf, int64_t* ipiv,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* x, int64_t ldx,
                           float* rcond, float* ferr, float* berr );
int64_t LAPACKE_zsysvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* af,
                           int64_t ldaf, int64_t* ipiv,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* x, int64_t ldx,
                           double* rcond, double* ferr, double* berr );

int64_t LAPACKE_skysvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           float* af, int64_t ldaf, int64_t* ipiv,
                           const float* b, int64_t ldb, float* x,
                           int64_t ldx, float* rcond, float* ferr,
                           float* berr );
int64_t LAPACKE_dkysvx_64( int matrix_layout, char fact, char uplo, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           double* af, int64_t ldaf, int64_t* ipiv,
                           const double* b, int64_t ldb, double* x,
                           int64_t ldx, double* rcond, double* ferr,
                           double* berr );

int64_t LAPACKE_ssysvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs, float* a,
                            int64_t lda, float* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, float* s, float* b,
                            int64_t ldb, float* x, int64_t ldx,
                            float* rcond, float* rpvgrw, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_dsysvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs, double* a,
                            int64_t lda, double* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, double* s, double* b,
                            int64_t ldb, double* x, int64_t ldx,
                            double* rcond, double* rpvgrw, double* berr,
                            int64_t n_err_bnds, double* err_bnds_norm,
                            double* err_bnds_comp, int64_t nparams,
                            double* params );
int64_t LAPACKE_csysvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, float* s,
                            lapack_complex_float* b, int64_t ldb,
                            lapack_complex_float* x, int64_t ldx,
                            float* rcond, float* rpvgrw, float* berr,
                            int64_t n_err_bnds, float* err_bnds_norm,
                            float* err_bnds_comp, int64_t nparams,
                            float* params );
int64_t LAPACKE_zsysvxx_64( int matrix_layout, char fact, char uplo,
                            int64_t n, int64_t nrhs,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* af, int64_t ldaf,
                            int64_t* ipiv, char* equed, double* s,
                            lapack_complex_double* b, int64_t ldb,
                            lapack_complex_double* x, int64_t ldx,
                            double* rcond, double* rpvgrw, double* berr,
                            int64_t n_err_bnds, double* err_bnds_norm,
                            double* err_bnds_comp, int64_t nparams,
                            double* params );

int64_t LAPACKE_ssytrd_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, float* d, float* e, float* tau );
int64_t LAPACKE_dsytrd_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, double* d, double* e, double* tau );

int64_t LAPACKE_skytrd_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, float* e, float* tau );
int64_t LAPACKE_dkytrd_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, double* e, double* tau );

int64_t LAPACKE_ssytrf_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dsytrf_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, int64_t* ipiv );
int64_t LAPACKE_csytrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_zsytrf_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv );

int64_t LAPACKE_skytrf_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dkytrf_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, int64_t* ipiv );

int64_t LAPACKE_ssytri_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, const int64_t* ipiv );
int64_t LAPACKE_dsytri_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, const int64_t* ipiv );
int64_t LAPACKE_csytri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           const int64_t* ipiv );
int64_t LAPACKE_zsytri_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           const int64_t* ipiv );

int64_t LAPACKE_skytri_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, const int64_t* ipiv );
int64_t LAPACKE_dkytri_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, const int64_t* ipiv );

int64_t LAPACKE_ssytrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dsytrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_csytrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zsytrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_skytrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dkytrs_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           const int64_t* ipiv, double* b, int64_t ldb );

int64_t LAPACKE_stbcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, int64_t kd, const float* ab,
                           int64_t ldab, float* rcond );
int64_t LAPACKE_dtbcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, int64_t kd, const double* ab,
                           int64_t ldab, double* rcond );
int64_t LAPACKE_ctbcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, int64_t kd,
                           const lapack_complex_float* ab, int64_t ldab,
                           float* rcond );
int64_t LAPACKE_ztbcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, int64_t kd,
                           const lapack_complex_double* ab, int64_t ldab,
                           double* rcond );

int64_t LAPACKE_stbrfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t kd, int64_t nrhs,
                           const float* ab, int64_t ldab, const float* b,
                           int64_t ldb, const float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_dtbrfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t kd, int64_t nrhs,
                           const double* ab, int64_t ldab, const double* b,
                           int64_t ldb, const double* x, int64_t ldx,
                           double* ferr, double* berr );
int64_t LAPACKE_ctbrfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t kd, int64_t nrhs,
                           const lapack_complex_float* ab, int64_t ldab,
                           const lapack_complex_float* b, int64_t ldb,
                           const lapack_complex_float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_ztbrfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t kd, int64_t nrhs,
                           const lapack_complex_double* ab, int64_t ldab,
                           const lapack_complex_double* b, int64_t ldb,
                           const lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_stbtrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t kd, int64_t nrhs,
                           const float* ab, int64_t ldab, float* b,
                           int64_t ldb );
int64_t LAPACKE_dtbtrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t kd, int64_t nrhs,
                           const double* ab, int64_t ldab, double* b,
                           int64_t ldb );
int64_t LAPACKE_ctbtrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t kd, int64_t nrhs,
                           const lapack_complex_float* ab, int64_t ldab,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_ztbtrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t kd, int64_t nrhs,
                           const lapack_complex_double* ab, int64_t ldab,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_stfsm_64( int matrix_layout, char transr, char side, char uplo,
                          char trans, char diag, int64_t m, int64_t n,
                          float alpha, const float* a, float* b,
                          int64_t ldb );
int64_t LAPACKE_dtfsm_64( int matrix_layout, char transr, char side, char uplo,
                          char trans, char diag, int64_t m, int64_t n,
                          double alpha, const double* a, double* b,
                          int64_t ldb );
int64_t LAPACKE_ctfsm_64( int matrix_layout, char transr, char side, char uplo,
                          char trans, char diag, int64_t m, int64_t n,
                          lapack_complex_float alpha,
                          const lapack_complex_float* a,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_ztfsm_64( int matrix_layout, char transr, char side, char uplo,
                          char trans, char diag, int64_t m, int64_t n,
                          lapack_complex_double alpha,
                          const lapack_complex_double* a,
                          lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_stftri_64( int matrix_layout, char transr, char uplo, char diag,
                           int64_t n, float* a );
int64_t LAPACKE_dtftri_64( int matrix_layout, char transr, char uplo, char diag,
                           int64_t n, double* a );
int64_t LAPACKE_ctftri_64( int matrix_layout, char transr, char uplo, char diag,
                           int64_t n, lapack_complex_float* a );
int64_t LAPACKE_ztftri_64( int matrix_layout, char transr, char uplo, char diag,
                           int64_t n, lapack_complex_double* a );

int64_t LAPACKE_stfttp_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const float* arf, float* ap );
int64_t LAPACKE_dtfttp_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const double* arf, double* ap );
int64_t LAPACKE_ctfttp_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const lapack_complex_float* arf,
                           lapack_complex_float* ap );
int64_t LAPACKE_ztfttp_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const lapack_complex_double* arf,
                           lapack_complex_double* ap );

int64_t LAPACKE_stfttr_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const float* arf, float* a,
                           int64_t lda );
int64_t LAPACKE_dtfttr_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const double* arf, double* a,
                           int64_t lda );
int64_t LAPACKE_ctfttr_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const lapack_complex_float* arf,
                           lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_ztfttr_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const lapack_complex_double* arf,
                           lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_stgevc_64( int matrix_layout, char side, char howmny,
                           const lapack_logical* select, int64_t n,
                           const float* s, int64_t lds, const float* p,
                           int64_t ldp, float* vl, int64_t ldvl,
                           float* vr, int64_t ldvr, int64_t mm,
                           int64_t* m );
int64_t LAPACKE_dtgevc_64( int matrix_layout, char side, char howmny,
                           const lapack_logical* select, int64_t n,
                           const double* s, int64_t lds, const double* p,
                           int64_t ldp, double* vl, int64_t ldvl,
                           double* vr, int64_t ldvr, int64_t mm,
                           int64_t* m );
int64_t LAPACKE_ctgevc_64( int matrix_layout, char side, char howmny,
                           const lapack_logical* select, int64_t n,
                           const lapack_complex_float* s, int64_t lds,
                           const lapack_complex_float* p, int64_t ldp,
                           lapack_complex_float* vl, int64_t ldvl,
                           lapack_complex_float* vr, int64_t ldvr,
                           int64_t mm, int64_t* m );
int64_t LAPACKE_ztgevc_64( int matrix_layout, char side, char howmny,
                           const lapack_logical* select, int64_t n,
                           const lapack_complex_double* s, int64_t lds,
                           const lapack_complex_double* p, int64_t ldp,
                           lapack_complex_double* vl, int64_t ldvl,
                           lapack_complex_double* vr, int64_t ldvr,
                           int64_t mm, int64_t* m );

int64_t LAPACKE_stgexc_64( int matrix_layout, lapack_logical wantq,
                           lapack_logical wantz, int64_t n, float* a,
                           int64_t lda, float* b, int64_t ldb, float* q,
                           int64_t ldq, float* z, int64_t ldz,
                           int64_t* ifst, int64_t* ilst );
int64_t LAPACKE_dtgexc_64( int matrix_layout, lapack_logical wantq,
                           lapack_logical wantz, int64_t n, double* a,
                           int64_t lda, double* b, int64_t ldb, double* q,
                           int64_t ldq, double* z, int64_t ldz,
                           int64_t* ifst, int64_t* ilst );
int64_t LAPACKE_ctgexc_64( int matrix_layout, lapack_logical wantq,
                           lapack_logical wantz, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* q, int64_t ldq,
                           lapack_complex_float* z, int64_t ldz,
                           int64_t ifst, int64_t ilst );
int64_t LAPACKE_ztgexc_64( int matrix_layout, lapack_logical wantq,
                           lapack_logical wantz, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* q, int64_t ldq,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t ifst, int64_t ilst );

int64_t LAPACKE_stgsen_64( int matrix_layout, int64_t ijob,
                           lapack_logical wantq, lapack_logical wantz,
                           const lapack_logical* select, int64_t n, float* a,
                           int64_t lda, float* b, int64_t ldb,
                           float* alphar, float* alphai, float* beta, float* q,
                           int64_t ldq, float* z, int64_t ldz,
                           int64_t* m, float* pl, float* pr, float* dif );
int64_t LAPACKE_dtgsen_64( int matrix_layout, int64_t ijob,
                           lapack_logical wantq, lapack_logical wantz,
                           const lapack_logical* select, int64_t n,
                           double* a, int64_t lda, double* b, int64_t ldb,
                           double* alphar, double* alphai, double* beta,
                           double* q, int64_t ldq, double* z, int64_t ldz,
                           int64_t* m, double* pl, double* pr, double* dif );
int64_t LAPACKE_ctgsen_64( int matrix_layout, int64_t ijob,
                           lapack_logical wantq, lapack_logical wantz,
                           const lapack_logical* select, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* alpha,
                           lapack_complex_float* beta, lapack_complex_float* q,
                           int64_t ldq, lapack_complex_float* z,
                           int64_t ldz, int64_t* m, float* pl, float* pr,
                           float* dif );
int64_t LAPACKE_ztgsen_64( int matrix_layout, int64_t ijob,
                           lapack_logical wantq, lapack_logical wantz,
                           const lapack_logical* select, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* alpha,
                           lapack_complex_double* beta,
                           lapack_complex_double* q, int64_t ldq,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t* m, double* pl, double* pr, double* dif );

int64_t LAPACKE_stgsja_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t p, int64_t n,
                           int64_t k, int64_t l, float* a, int64_t lda,
                           float* b, int64_t ldb, float tola, float tolb,
                           float* alpha, float* beta, float* u, int64_t ldu,
                           float* v, int64_t ldv, float* q, int64_t ldq,
                           int64_t* ncycle );
int64_t LAPACKE_dtgsja_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t p, int64_t n,
                           int64_t k, int64_t l, double* a,
                           int64_t lda, double* b, int64_t ldb,
                           double tola, double tolb, double* alpha,
                           double* beta, double* u, int64_t ldu, double* v,
                           int64_t ldv, double* q, int64_t ldq,
                           int64_t* ncycle );
int64_t LAPACKE_ctgsja_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t p, int64_t n,
                           int64_t k, int64_t l, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* b,
                           int64_t ldb, float tola, float tolb, float* alpha,
                           float* beta, lapack_complex_float* u, int64_t ldu,
                           lapack_complex_float* v, int64_t ldv,
                           lapack_complex_float* q, int64_t ldq,
                           int64_t* ncycle );
int64_t LAPACKE_ztgsja_64( int matrix_layout, char jobu, char jobv, char jobq,
                           int64_t m, int64_t p, int64_t n,
                           int64_t k, int64_t l, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* b,
                           int64_t ldb, double tola, double tolb,
                           double* alpha, double* beta,
                           lapack_complex_double* u, int64_t ldu,
                           lapack_complex_double* v, int64_t ldv,
                           lapack_complex_double* q, int64_t ldq,
                           int64_t* ncycle );

int64_t LAPACKE_stgsna_64( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, int64_t n,
                           const float* a, int64_t lda, const float* b,
                           int64_t ldb, const float* vl, int64_t ldvl,
                           const float* vr, int64_t ldvr, float* s,
                           float* dif, int64_t mm, int64_t* m );
int64_t LAPACKE_dtgsna_64( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, int64_t n,
                           const double* a, int64_t lda, const double* b,
                           int64_t ldb, const double* vl, int64_t ldvl,
                           const double* vr, int64_t ldvr, double* s,
                           double* dif, int64_t mm, int64_t* m );
int64_t LAPACKE_ctgsna_64( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* b, int64_t ldb,
                           const lapack_complex_float* vl, int64_t ldvl,
                           const lapack_complex_float* vr, int64_t ldvr,
                           float* s, float* dif, int64_t mm, int64_t* m );
int64_t LAPACKE_ztgsna_64( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* b, int64_t ldb,
                           const lapack_complex_double* vl, int64_t ldvl,
                           const lapack_complex_double* vr, int64_t ldvr,
                           double* s, double* dif, int64_t mm,
                           int64_t* m );

int64_t LAPACKE_stgsyl_64( int matrix_layout, char trans, int64_t ijob,
                           int64_t m, int64_t n, const float* a,
                           int64_t lda, const float* b, int64_t ldb,
                           float* c, int64_t ldc, const float* d,
                           int64_t ldd, const float* e, int64_t lde,
                           float* f, int64_t ldf, float* scale, float* dif );
int64_t LAPACKE_dtgsyl_64( int matrix_layout, char trans, int64_t ijob,
                           int64_t m, int64_t n, const double* a,
                           int64_t lda, const double* b, int64_t ldb,
                           double* c, int64_t ldc, const double* d,
                           int64_t ldd, const double* e, int64_t lde,
                           double* f, int64_t ldf, double* scale,
                           double* dif );
int64_t LAPACKE_ctgsyl_64( int matrix_layout, char trans, int64_t ijob,
                           int64_t m, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* c, int64_t ldc,
                           const lapack_complex_float* d, int64_t ldd,
                           const lapack_complex_float* e, int64_t lde,
                           lapack_complex_float* f, int64_t ldf,
                           float* scale, float* dif );
int64_t LAPACKE_ztgsyl_64( int matrix_layout, char trans, int64_t ijob,
                           int64_t m, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* c, int64_t ldc,
                           const lapack_complex_double* d, int64_t ldd,
                           const lapack_complex_double* e, int64_t lde,
                           lapack_complex_double* f, int64_t ldf,
                           double* scale, double* dif );

int64_t LAPACKE_stpcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, const float* ap, float* rcond );
int64_t LAPACKE_dtpcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, const double* ap, double* rcond );
int64_t LAPACKE_ctpcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, const lapack_complex_float* ap,
                           float* rcond );
int64_t LAPACKE_ztpcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, const lapack_complex_double* ap,
                           double* rcond );

int64_t LAPACKE_stprfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs, const float* ap,
                           const float* b, int64_t ldb, const float* x,
                           int64_t ldx, float* ferr, float* berr );
int64_t LAPACKE_dtprfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs, const double* ap,
                           const double* b, int64_t ldb, const double* x,
                           int64_t ldx, double* ferr, double* berr );
int64_t LAPACKE_ctprfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_float* ap,
                           const lapack_complex_float* b, int64_t ldb,
                           const lapack_complex_float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_ztprfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_double* ap,
                           const lapack_complex_double* b, int64_t ldb,
                           const lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_stptri_64( int matrix_layout, char uplo, char diag, int64_t n,
                           float* ap );
int64_t LAPACKE_dtptri_64( int matrix_layout, char uplo, char diag, int64_t n,
                           double* ap );
int64_t LAPACKE_ctptri_64( int matrix_layout, char uplo, char diag, int64_t n,
                           lapack_complex_float* ap );
int64_t LAPACKE_ztptri_64( int matrix_layout, char uplo, char diag, int64_t n,
                           lapack_complex_double* ap );

int64_t LAPACKE_stptrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs, const float* ap,
                           float* b, int64_t ldb );
int64_t LAPACKE_dtptrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs, const double* ap,
                           double* b, int64_t ldb );
int64_t LAPACKE_ctptrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_float* ap,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_ztptrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_double* ap,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_stpttf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const float* ap, float* arf );
int64_t LAPACKE_dtpttf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const double* ap, double* arf );
int64_t LAPACKE_ctpttf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const lapack_complex_float* ap,
                           lapack_complex_float* arf );
int64_t LAPACKE_ztpttf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const lapack_complex_double* ap,
                           lapack_complex_double* arf );

int64_t LAPACKE_stpttr_64( int matrix_layout, char uplo, int64_t n,
                           const float* ap, float* a, int64_t lda );
int64_t LAPACKE_dtpttr_64( int matrix_layout, char uplo, int64_t n,
                           const double* ap, double* a, int64_t lda );
int64_t LAPACKE_ctpttr_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* ap,
                           lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_ztpttr_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* ap,
                           lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_strcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, const float* a, int64_t lda,
                           float* rcond );
int64_t LAPACKE_dtrcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, const double* a, int64_t lda,
                           double* rcond );
int64_t LAPACKE_ctrcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, const lapack_complex_float* a,
                           int64_t lda, float* rcond );
int64_t LAPACKE_ztrcon_64( int matrix_layout, char norm, char uplo, char diag,
                           int64_t n, const lapack_complex_double* a,
                           int64_t lda, double* rcond );

int64_t LAPACKE_strevc_64( int matrix_layout, char side, char howmny,
                           lapack_logical* select, int64_t n, const float* t,
                           int64_t ldt, float* vl, int64_t ldvl,
                           float* vr, int64_t ldvr, int64_t mm,
                           int64_t* m );
int64_t LAPACKE_dtrevc_64( int matrix_layout, char side, char howmny,
                           lapack_logical* select, int64_t n,
                           const double* t, int64_t ldt, double* vl,
                           int64_t ldvl, double* vr, int64_t ldvr,
                           int64_t mm, int64_t* m );
int64_t LAPACKE_ctrevc_64( int matrix_layout, char side, char howmny,
                           const lapack_logical* select, int64_t n,
                           lapack_complex_float* t, int64_t ldt,
                           lapack_complex_float* vl, int64_t ldvl,
                           lapack_complex_float* vr, int64_t ldvr,
                           int64_t mm, int64_t* m );
int64_t LAPACKE_ztrevc_64( int matrix_layout, char side, char howmny,
                           const lapack_logical* select, int64_t n,
                           lapack_complex_double* t, int64_t ldt,
                           lapack_complex_double* vl, int64_t ldvl,
                           lapack_complex_double* vr, int64_t ldvr,
                           int64_t mm, int64_t* m );

int64_t LAPACKE_strexc_64( int matrix_layout, char compq, int64_t n, float* t,
                           int64_t ldt, float* q, int64_t ldq,
                           int64_t* ifst, int64_t* ilst );
int64_t LAPACKE_dtrexc_64( int matrix_layout, char compq, int64_t n,
                           double* t, int64_t ldt, double* q, int64_t ldq,
                           int64_t* ifst, int64_t* ilst );
int64_t LAPACKE_ctrexc_64( int matrix_layout, char compq, int64_t n,
                           lapack_complex_float* t, int64_t ldt,
                           lapack_complex_float* q, int64_t ldq,
                           int64_t ifst, int64_t ilst );
int64_t LAPACKE_ztrexc_64( int matrix_layout, char compq, int64_t n,
                           lapack_complex_double* t, int64_t ldt,
                           lapack_complex_double* q, int64_t ldq,
                           int64_t ifst, int64_t ilst );

int64_t LAPACKE_strrfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs, const float* a,
                           int64_t lda, const float* b, int64_t ldb,
                           const float* x, int64_t ldx, float* ferr,
                           float* berr );
int64_t LAPACKE_dtrrfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs, const double* a,
                           int64_t lda, const double* b, int64_t ldb,
                           const double* x, int64_t ldx, double* ferr,
                           double* berr );
int64_t LAPACKE_ctrrfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* b, int64_t ldb,
                           const lapack_complex_float* x, int64_t ldx,
                           float* ferr, float* berr );
int64_t LAPACKE_ztrrfs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* b, int64_t ldb,
                           const lapack_complex_double* x, int64_t ldx,
                           double* ferr, double* berr );

int64_t LAPACKE_strsen_64( int matrix_layout, char job, char compq,
                           const lapack_logical* select, int64_t n, float* t,
                           int64_t ldt, float* q, int64_t ldq, float* wr,
                           float* wi, int64_t* m, float* s, float* sep );
int64_t LAPACKE_dtrsen_64( int matrix_layout, char job, char compq,
                           const lapack_logical* select, int64_t n,
                           double* t, int64_t ldt, double* q, int64_t ldq,
                           double* wr, double* wi, int64_t* m, double* s,
                           double* sep );
int64_t LAPACKE_ctrsen_64( int matrix_layout, char job, char compq,
                           const lapack_logical* select, int64_t n,
                           lapack_complex_float* t, int64_t ldt,
                           lapack_complex_float* q, int64_t ldq,
                           lapack_complex_float* w, int64_t* m, float* s,
                           float* sep );
int64_t LAPACKE_ztrsen_64( int matrix_layout, char job, char compq,
                           const lapack_logical* select, int64_t n,
                           lapack_complex_double* t, int64_t ldt,
                           lapack_complex_double* q, int64_t ldq,
                           lapack_complex_double* w, int64_t* m, double* s,
                           double* sep );

int64_t LAPACKE_strsna_64( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, int64_t n,
                           const float* t, int64_t ldt, const float* vl,
                           int64_t ldvl, const float* vr, int64_t ldvr,
                           float* s, float* sep, int64_t mm, int64_t* m );
int64_t LAPACKE_dtrsna_64( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, int64_t n,
                           const double* t, int64_t ldt, const double* vl,
                           int64_t ldvl, const double* vr, int64_t ldvr,
                           double* s, double* sep, int64_t mm,
                           int64_t* m );
int64_t LAPACKE_ctrsna_64( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, int64_t n,
                           const lapack_complex_float* t, int64_t ldt,
                           const lapack_complex_float* vl, int64_t ldvl,
                           const lapack_complex_float* vr, int64_t ldvr,
                           float* s, float* sep, int64_t mm, int64_t* m );
int64_t LAPACKE_ztrsna_64( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, int64_t n,
                           const lapack_complex_double* t, int64_t ldt,
                           const lapack_complex_double* vl, int64_t ldvl,
                           const lapack_complex_double* vr, int64_t ldvr,
                           double* s, double* sep, int64_t mm,
                           int64_t* m );

int64_t LAPACKE_strsyl_64( int matrix_layout, char trana, char tranb,
                           int64_t isgn, int64_t m, int64_t n,
                           const float* a, int64_t lda, const float* b,
                           int64_t ldb, float* c, int64_t ldc,
                           float* scale );
int64_t LAPACKE_dtrsyl_64( int matrix_layout, char trana, char tranb,
                           int64_t isgn, int64_t m, int64_t n,
                           const double* a, int64_t lda, const double* b,
                           int64_t ldb, double* c, int64_t ldc,
                           double* scale );
int64_t LAPACKE_ctrsyl_64( int matrix_layout, char trana, char tranb,
                           int64_t isgn, int64_t m, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* c, int64_t ldc,
                           float* scale );
int64_t LAPACKE_ztrsyl_64( int matrix_layout, char trana, char tranb,
                           int64_t isgn, int64_t m, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* c, int64_t ldc,
                           double* scale );

int64_t LAPACKE_strsyl3_64( int matrix_layout, char trana, char tranb,
                            int64_t isgn, int64_t m, int64_t n,
                            const float* a, int64_t lda, const float* b,
                            int64_t ldb, float* c, int64_t ldc,
                            float* scale );
int64_t LAPACKE_dtrsyl3_64( int matrix_layout, char trana, char tranb,
                            int64_t isgn, int64_t m, int64_t n,
                            const double* a, int64_t lda, const double* b,
                            int64_t ldb, double* c, int64_t ldc,
                            double* scale );
int64_t LAPACKE_ztrsyl3_64( int matrix_layout, char trana, char tranb,
                            int64_t isgn, int64_t m, int64_t n,
                            const lapack_complex_double* a, int64_t lda,
                            const lapack_complex_double* b, int64_t ldb,
                            lapack_complex_double* c, int64_t ldc,
                            double* scale );

int64_t LAPACKE_strtri_64( int matrix_layout, char uplo, char diag, int64_t n,
                           float* a, int64_t lda );
int64_t LAPACKE_dtrtri_64( int matrix_layout, char uplo, char diag, int64_t n,
                           double* a, int64_t lda );
int64_t LAPACKE_ctrtri_64( int matrix_layout, char uplo, char diag, int64_t n,
                           lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_ztrtri_64( int matrix_layout, char uplo, char diag, int64_t n,
                           lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_strtrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs, const float* a,
                           int64_t lda, float* b, int64_t ldb );
int64_t LAPACKE_dtrtrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs, const double* a,
                           int64_t lda, double* b, int64_t ldb );
int64_t LAPACKE_ctrtrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_ztrtrs_64( int matrix_layout, char uplo, char trans, char diag,
                           int64_t n, int64_t nrhs,
                           const lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_strttf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const float* a, int64_t lda,
                           float* arf );
int64_t LAPACKE_dtrttf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const double* a, int64_t lda,
                           double* arf );
int64_t LAPACKE_ctrttf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* arf );
int64_t LAPACKE_ztrttf_64( int matrix_layout, char transr, char uplo,
                           int64_t n, const lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* arf );

int64_t LAPACKE_strttp_64( int matrix_layout, char uplo, int64_t n,
                           const float* a, int64_t lda, float* ap );
int64_t LAPACKE_dtrttp_64( int matrix_layout, char uplo, int64_t n,
                           const double* a, int64_t lda, double* ap );
int64_t LAPACKE_ctrttp_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* ap );
int64_t LAPACKE_ztrttp_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* ap );

int64_t LAPACKE_stzrzf_64( int matrix_layout, int64_t m, int64_t n,
                           float* a, int64_t lda, float* tau );
int64_t LAPACKE_dtzrzf_64( int matrix_layout, int64_t m, int64_t n,
                           double* a, int64_t lda, double* tau );
int64_t LAPACKE_ctzrzf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* tau );
int64_t LAPACKE_ztzrzf_64( int matrix_layout, int64_t m, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* tau );

int64_t LAPACKE_cungbr_64( int matrix_layout, char vect, int64_t m,
                           int64_t n, int64_t k, lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* tau );
int64_t LAPACKE_zungbr_64( int matrix_layout, char vect, int64_t m,
                           int64_t n, int64_t k, lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* tau );

int64_t LAPACKE_cunghr_64( int matrix_layout, int64_t n, int64_t ilo,
                           int64_t ihi, lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* tau );
int64_t LAPACKE_zunghr_64( int matrix_layout, int64_t n, int64_t ilo,
                           int64_t ihi, lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* tau );

int64_t LAPACKE_cunglq_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* tau );
int64_t LAPACKE_zunglq_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* tau );

int64_t LAPACKE_cungql_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* tau );
int64_t LAPACKE_zungql_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* tau );

int64_t LAPACKE_cungqr_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* tau );
int64_t LAPACKE_zungqr_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* tau );

int64_t LAPACKE_cungrq_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* tau );
int64_t LAPACKE_zungrq_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t k, lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* tau );

int64_t LAPACKE_cungtr_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* tau );
int64_t LAPACKE_zungtr_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* tau );

int64_t LAPACKE_cungtsqr_row_64( int matrix_layout, int64_t m, int64_t n,
                                 int64_t mb, int64_t nb,
                                 lapack_complex_float* a, int64_t lda,
                                 const lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_zungtsqr_row_64( int matrix_layout, int64_t m, int64_t n,
                                 int64_t mb, int64_t nb,
                                 lapack_complex_double* a, int64_t lda,
                                 const lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_cunmbr_64( int matrix_layout, char vect, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* tau,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zunmbr_64( int matrix_layout, char vect, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_cunmhr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t ilo,
                           int64_t ihi, const lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* tau,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zunmhr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t ilo,
                           int64_t ihi, const lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* tau,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_cunmlq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* tau,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zunmlq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_cunmql_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* tau,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zunmql_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_cunmqr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* tau,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zunmqr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_cunmrq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* tau,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zunmrq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_cunmrz_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           int64_t l, const lapack_complex_float* a,
                           int64_t lda, const lapack_complex_float* tau,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zunmrz_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           int64_t l, const lapack_complex_double* a,
                           int64_t lda, const lapack_complex_double* tau,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_cunmtr_64( int matrix_layout, char side, char uplo, char trans,
                           int64_t m, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* tau,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zunmtr_64( int matrix_layout, char side, char uplo, char trans,
                           int64_t m, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_cupgtr_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* ap,
                           const lapack_complex_float* tau,
                           lapack_complex_float* q, int64_t ldq );
int64_t LAPACKE_zupgtr_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* ap,
                           const lapack_complex_double* tau,
                           lapack_complex_double* q, int64_t ldq );

int64_t LAPACKE_cupmtr_64( int matrix_layout, char side, char uplo, char trans,
                           int64_t m, int64_t n,
                           const lapack_complex_float* ap,
                           const lapack_complex_float* tau,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zupmtr_64( int matrix_layout, char side, char uplo, char trans,
                           int64_t m, int64_t n,
                           const lapack_complex_double* ap,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_sbdsdc_work_64( int matrix_layout, char uplo, char compq,
                                int64_t n, float* d, float* e, float* u,
                                int64_t ldu, float* vt, int64_t ldvt,
                                float* q, int64_t* iq, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dbdsdc_work_64( int matrix_layout, char uplo, char compq,
                                int64_t n, double* d, double* e, double* u,
                                int64_t ldu, double* vt, int64_t ldvt,
                                double* q, int64_t* iq, double* work,
                                int64_t* iwork );

int64_t LAPACKE_sbdsvdx_work_64( int matrix_layout, char uplo, char jobz, char range,
                                 int64_t n, float* d, float* e,
                                 float vl, float vu,
                                 int64_t il, int64_t iu, int64_t* ns,
                                 float* s, float* z, int64_t ldz,
                                 float* work, int64_t* iwork );
int64_t LAPACKE_dbdsvdx_work_64( int matrix_layout, char uplo, char jobz, char range,
                                 int64_t n, double* d, double* e,
                                 double vl, double vu,
                                 int64_t il, int64_t iu, int64_t* ns,
                                 double* s, double* z, int64_t ldz,
                                 double* work, int64_t* iwork );

int64_t LAPACKE_sbdsqr_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t ncvt, int64_t nru, int64_t ncc,
                                float* d, float* e, float* vt, int64_t ldvt,
                                float* u, int64_t ldu, float* c,
                                int64_t ldc, float* work );
int64_t LAPACKE_dbdsqr_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t ncvt, int64_t nru, int64_t ncc,
                                double* d, double* e, double* vt,
                                int64_t ldvt, double* u, int64_t ldu,
                                double* c, int64_t ldc, double* work );
int64_t LAPACKE_cbdsqr_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t ncvt, int64_t nru, int64_t ncc,
                                float* d, float* e, lapack_complex_float* vt,
                                int64_t ldvt, lapack_complex_float* u,
                                int64_t ldu, lapack_complex_float* c,
                                int64_t ldc, float* work );
int64_t LAPACKE_zbdsqr_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t ncvt, int64_t nru, int64_t ncc,
                                double* d, double* e, lapack_complex_double* vt,
                                int64_t ldvt, lapack_complex_double* u,
                                int64_t ldu, lapack_complex_double* c,
                                int64_t ldc, double* work );

int64_t LAPACKE_sdisna_work_64( char job, int64_t m, int64_t n,
                                const float* d, float* sep );
int64_t LAPACKE_ddisna_work_64( char job, int64_t m, int64_t n,
                                const double* d, double* sep );

int64_t LAPACKE_sgbbrd_work_64( int matrix_layout, char vect, int64_t m,
                                int64_t n, int64_t ncc, int64_t kl,
                                int64_t ku, float* ab, int64_t ldab,
                                float* d, float* e, float* q, int64_t ldq,
                                float* pt, int64_t ldpt, float* c,
                                int64_t ldc, float* work );
int64_t LAPACKE_dgbbrd_work_64( int matrix_layout, char vect, int64_t m,
                                int64_t n, int64_t ncc, int64_t kl,
                                int64_t ku, double* ab, int64_t ldab,
                                double* d, double* e, double* q, int64_t ldq,
                                double* pt, int64_t ldpt, double* c,
                                int64_t ldc, double* work );
int64_t LAPACKE_cgbbrd_work_64( int matrix_layout, char vect, int64_t m,
                                int64_t n, int64_t ncc, int64_t kl,
                                int64_t ku, lapack_complex_float* ab,
                                int64_t ldab, float* d, float* e,
                                lapack_complex_float* q, int64_t ldq,
                                lapack_complex_float* pt, int64_t ldpt,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zgbbrd_work_64( int matrix_layout, char vect, int64_t m,
                                int64_t n, int64_t ncc, int64_t kl,
                                int64_t ku, lapack_complex_double* ab,
                                int64_t ldab, double* d, double* e,
                                lapack_complex_double* q, int64_t ldq,
                                lapack_complex_double* pt, int64_t ldpt,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sgbcon_work_64( int matrix_layout, char norm, int64_t n,
                                int64_t kl, int64_t ku, const float* ab,
                                int64_t ldab, const int64_t* ipiv,
                                float anorm, float* rcond, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dgbcon_work_64( int matrix_layout, char norm, int64_t n,
                                int64_t kl, int64_t ku, const double* ab,
                                int64_t ldab, const int64_t* ipiv,
                                double anorm, double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_cgbcon_work_64( int matrix_layout, char norm, int64_t n,
                                int64_t kl, int64_t ku,
                                const lapack_complex_float* ab, int64_t ldab,
                                const int64_t* ipiv, float anorm,
                                float* rcond, lapack_complex_float* work,
                                float* rwork );
int64_t LAPACKE_zgbcon_work_64( int matrix_layout, char norm, int64_t n,
                                int64_t kl, int64_t ku,
                                const lapack_complex_double* ab,
                                int64_t ldab, const int64_t* ipiv,
                                double anorm, double* rcond,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sgbequ_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku, const float* ab,
                                int64_t ldab, float* r, float* c,
                                float* rowcnd, float* colcnd, float* amax );
int64_t LAPACKE_dgbequ_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku, const double* ab,
                                int64_t ldab, double* r, double* c,
                                double* rowcnd, double* colcnd, double* amax );
int64_t LAPACKE_cgbequ_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku,
                                const lapack_complex_float* ab, int64_t ldab,
                                float* r, float* c, float* rowcnd,
                                float* colcnd, float* amax );
int64_t LAPACKE_zgbequ_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku,
                                const lapack_complex_double* ab,
                                int64_t ldab, double* r, double* c,
                                double* rowcnd, double* colcnd, double* amax );

int64_t LAPACKE_sgbequb_work_64( int matrix_layout, int64_t m, int64_t n,
                                 int64_t kl, int64_t ku, const float* ab,
                                 int64_t ldab, float* r, float* c,
                                 float* rowcnd, float* colcnd, float* amax );
int64_t LAPACKE_dgbequb_work_64( int matrix_layout, int64_t m, int64_t n,
                                 int64_t kl, int64_t ku, const double* ab,
                                 int64_t ldab, double* r, double* c,
                                 double* rowcnd, double* colcnd, double* amax );
int64_t LAPACKE_cgbequb_work_64( int matrix_layout, int64_t m, int64_t n,
                                 int64_t kl, int64_t ku,
                                 const lapack_complex_float* ab,
                                 int64_t ldab, float* r, float* c,
                                 float* rowcnd, float* colcnd, float* amax );
int64_t LAPACKE_zgbequb_work_64( int matrix_layout, int64_t m, int64_t n,
                                 int64_t kl, int64_t ku,
                                 const lapack_complex_double* ab,
                                 int64_t ldab, double* r, double* c,
                                 double* rowcnd, double* colcnd, double* amax );

int64_t LAPACKE_sgbrfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t kl, int64_t ku, int64_t nrhs,
                                const float* ab, int64_t ldab,
                                const float* afb, int64_t ldafb,
                                const int64_t* ipiv, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* ferr, float* berr, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dgbrfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t kl, int64_t ku, int64_t nrhs,
                                const double* ab, int64_t ldab,
                                const double* afb, int64_t ldafb,
                                const int64_t* ipiv, const double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* ferr, double* berr, double* work,
                                int64_t* iwork );
int64_t LAPACKE_cgbrfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t kl, int64_t ku, int64_t nrhs,
                                const lapack_complex_float* ab, int64_t ldab,
                                const lapack_complex_float* afb,
                                int64_t ldafb, const int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zgbrfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t kl, int64_t ku, int64_t nrhs,
                                const lapack_complex_double* ab,
                                int64_t ldab,
                                const lapack_complex_double* afb,
                                int64_t ldafb, const int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sgbrfsx_work_64( int matrix_layout, char trans, char equed,
                                 int64_t n, int64_t kl, int64_t ku,
                                 int64_t nrhs, const float* ab,
                                 int64_t ldab, const float* afb,
                                 int64_t ldafb, const int64_t* ipiv,
                                 const float* r, const float* c, const float* b,
                                 int64_t ldb, float* x, int64_t ldx,
                                 float* rcond, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, float* work,
                                 int64_t* iwork );
int64_t LAPACKE_dgbrfsx_work_64( int matrix_layout, char trans, char equed,
                                 int64_t n, int64_t kl, int64_t ku,
                                 int64_t nrhs, const double* ab,
                                 int64_t ldab, const double* afb,
                                 int64_t ldafb, const int64_t* ipiv,
                                 const double* r, const double* c,
                                 const double* b, int64_t ldb, double* x,
                                 int64_t ldx, double* rcond, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, double* work,
                                 int64_t* iwork );
int64_t LAPACKE_cgbrfsx_work_64( int matrix_layout, char trans, char equed,
                                 int64_t n, int64_t kl, int64_t ku,
                                 int64_t nrhs,
                                 const lapack_complex_float* ab,
                                 int64_t ldab,
                                 const lapack_complex_float* afb,
                                 int64_t ldafb, const int64_t* ipiv,
                                 const float* r, const float* c,
                                 const lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* x, int64_t ldx,
                                 float* rcond, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, lapack_complex_float* work,
                                 float* rwork );
int64_t LAPACKE_zgbrfsx_work_64( int matrix_layout, char trans, char equed,
                                 int64_t n, int64_t kl, int64_t ku,
                                 int64_t nrhs,
                                 const lapack_complex_double* ab,
                                 int64_t ldab,
                                 const lapack_complex_double* afb,
                                 int64_t ldafb, const int64_t* ipiv,
                                 const double* r, const double* c,
                                 const lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* x, int64_t ldx,
                                 double* rcond, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, lapack_complex_double* work,
                                 double* rwork );

int64_t LAPACKE_sgbsv_work_64( int matrix_layout, int64_t n, int64_t kl,
                               int64_t ku, int64_t nrhs, float* ab,
                               int64_t ldab, int64_t* ipiv, float* b,
                               int64_t ldb );
int64_t LAPACKE_dgbsv_work_64( int matrix_layout, int64_t n, int64_t kl,
                               int64_t ku, int64_t nrhs, double* ab,
                               int64_t ldab, int64_t* ipiv, double* b,
                               int64_t ldb );
int64_t LAPACKE_cgbsv_work_64( int matrix_layout, int64_t n, int64_t kl,
                               int64_t ku, int64_t nrhs,
                               lapack_complex_float* ab, int64_t ldab,
                               int64_t* ipiv, lapack_complex_float* b,
                               int64_t ldb );
int64_t LAPACKE_zgbsv_work_64( int matrix_layout, int64_t n, int64_t kl,
                               int64_t ku, int64_t nrhs,
                               lapack_complex_double* ab, int64_t ldab,
                               int64_t* ipiv, lapack_complex_double* b,
                               int64_t ldb );

int64_t LAPACKE_sgbsvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t kl, int64_t ku,
                                int64_t nrhs, float* ab, int64_t ldab,
                                float* afb, int64_t ldafb, int64_t* ipiv,
                                char* equed, float* r, float* c, float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dgbsvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t kl, int64_t ku,
                                int64_t nrhs, double* ab, int64_t ldab,
                                double* afb, int64_t ldafb, int64_t* ipiv,
                                char* equed, double* r, double* c, double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cgbsvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t kl, int64_t ku,
                                int64_t nrhs, lapack_complex_float* ab,
                                int64_t ldab, lapack_complex_float* afb,
                                int64_t ldafb, int64_t* ipiv, char* equed,
                                float* r, float* c, lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* x,
                                int64_t ldx, float* rcond, float* ferr,
                                float* berr, lapack_complex_float* work,
                                float* rwork );
int64_t LAPACKE_zgbsvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t kl, int64_t ku,
                                int64_t nrhs, lapack_complex_double* ab,
                                int64_t ldab, lapack_complex_double* afb,
                                int64_t ldafb, int64_t* ipiv, char* equed,
                                double* r, double* c, lapack_complex_double* b,
                                int64_t ldb, lapack_complex_double* x,
                                int64_t ldx, double* rcond, double* ferr,
                                double* berr, lapack_complex_double* work,
                                double* rwork );

int64_t LAPACKE_sgbsvxx_work_64( int matrix_layout, char fact, char trans,
                                 int64_t n, int64_t kl, int64_t ku,
                                 int64_t nrhs, float* ab, int64_t ldab,
                                 float* afb, int64_t ldafb, int64_t* ipiv,
                                 char* equed, float* r, float* c, float* b,
                                 int64_t ldb, float* x, int64_t ldx,
                                 float* rcond, float* rpvgrw, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, float* work,
                                 int64_t* iwork );
int64_t LAPACKE_dgbsvxx_work_64( int matrix_layout, char fact, char trans,
                                 int64_t n, int64_t kl, int64_t ku,
                                 int64_t nrhs, double* ab, int64_t ldab,
                                 double* afb, int64_t ldafb,
                                 int64_t* ipiv, char* equed, double* r,
                                 double* c, double* b, int64_t ldb,
                                 double* x, int64_t ldx, double* rcond,
                                 double* rpvgrw, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, double* work,
                                 int64_t* iwork );
int64_t LAPACKE_cgbsvxx_work_64( int matrix_layout, char fact, char trans,
                                 int64_t n, int64_t kl, int64_t ku,
                                 int64_t nrhs, lapack_complex_float* ab,
                                 int64_t ldab, lapack_complex_float* afb,
                                 int64_t ldafb, int64_t* ipiv,
                                 char* equed, float* r, float* c,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* x, int64_t ldx,
                                 float* rcond, float* rpvgrw, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, lapack_complex_float* work,
                                 float* rwork );
int64_t LAPACKE_zgbsvxx_work_64( int matrix_layout, char fact, char trans,
                                 int64_t n, int64_t kl, int64_t ku,
                                 int64_t nrhs, lapack_complex_double* ab,
                                 int64_t ldab, lapack_complex_double* afb,
                                 int64_t ldafb, int64_t* ipiv,
                                 char* equed, double* r, double* c,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* x, int64_t ldx,
                                 double* rcond, double* rpvgrw, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, lapack_complex_double* work,
                                 double* rwork );

int64_t LAPACKE_sgbtrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku, float* ab,
                                int64_t ldab, int64_t* ipiv );
int64_t LAPACKE_dgbtrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku, double* ab,
                                int64_t ldab, int64_t* ipiv );
int64_t LAPACKE_cgbtrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku,
                                lapack_complex_float* ab, int64_t ldab,
                                int64_t* ipiv );
int64_t LAPACKE_zgbtrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku,
                                lapack_complex_double* ab, int64_t ldab,
                                int64_t* ipiv );

int64_t LAPACKE_sgbtrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t kl, int64_t ku, int64_t nrhs,
                                const float* ab, int64_t ldab,
                                const int64_t* ipiv, float* b,
                                int64_t ldb );
int64_t LAPACKE_dgbtrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t kl, int64_t ku, int64_t nrhs,
                                const double* ab, int64_t ldab,
                                const int64_t* ipiv, double* b,
                                int64_t ldb );
int64_t LAPACKE_cgbtrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t kl, int64_t ku, int64_t nrhs,
                                const lapack_complex_float* ab, int64_t ldab,
                                const int64_t* ipiv, lapack_complex_float* b,
                                int64_t ldb );
int64_t LAPACKE_zgbtrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t kl, int64_t ku, int64_t nrhs,
                                const lapack_complex_double* ab,
                                int64_t ldab, const int64_t* ipiv,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sgebak_work_64( int matrix_layout, char job, char side,
                                int64_t n, int64_t ilo, int64_t ihi,
                                const float* scale, int64_t m, float* v,
                                int64_t ldv );
int64_t LAPACKE_dgebak_work_64( int matrix_layout, char job, char side,
                                int64_t n, int64_t ilo, int64_t ihi,
                                const double* scale, int64_t m, double* v,
                                int64_t ldv );
int64_t LAPACKE_cgebak_work_64( int matrix_layout, char job, char side,
                                int64_t n, int64_t ilo, int64_t ihi,
                                const float* scale, int64_t m,
                                lapack_complex_float* v, int64_t ldv );
int64_t LAPACKE_zgebak_work_64( int matrix_layout, char job, char side,
                                int64_t n, int64_t ilo, int64_t ihi,
                                const double* scale, int64_t m,
                                lapack_complex_double* v, int64_t ldv );

int64_t LAPACKE_sgebal_work_64( int matrix_layout, char job, int64_t n,
                                float* a, int64_t lda, int64_t* ilo,
                                int64_t* ihi, float* scale );
int64_t LAPACKE_dgebal_work_64( int matrix_layout, char job, int64_t n,
                                double* a, int64_t lda, int64_t* ilo,
                                int64_t* ihi, double* scale );
int64_t LAPACKE_cgebal_work_64( int matrix_layout, char job, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ilo, int64_t* ihi,
                                float* scale );
int64_t LAPACKE_zgebal_work_64( int matrix_layout, char job, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ilo, int64_t* ihi,
                                double* scale );

int64_t LAPACKE_sgebrd_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, float* d, float* e,
                                float* tauq, float* taup, float* work,
                                int64_t lwork );
int64_t LAPACKE_dgebrd_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, double* d, double* e,
                                double* tauq, double* taup, double* work,
                                int64_t lwork );
int64_t LAPACKE_cgebrd_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                float* d, float* e, lapack_complex_float* tauq,
                                lapack_complex_float* taup,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgebrd_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                double* d, double* e,
                                lapack_complex_double* tauq,
                                lapack_complex_double* taup,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgecon_work_64( int matrix_layout, char norm, int64_t n,
                                const float* a, int64_t lda, float anorm,
                                float* rcond, float* work, int64_t* iwork );
int64_t LAPACKE_dgecon_work_64( int matrix_layout, char norm, int64_t n,
                                const double* a, int64_t lda, double anorm,
                                double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_cgecon_work_64( int matrix_layout, char norm, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                float anorm, float* rcond,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zgecon_work_64( int matrix_layout, char norm, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                double anorm, double* rcond,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sgeequ_work_64( int matrix_layout, int64_t m, int64_t n,
                                const float* a, int64_t lda, float* r,
                                float* c, float* rowcnd, float* colcnd,
                                float* amax );
int64_t LAPACKE_dgeequ_work_64( int matrix_layout, int64_t m, int64_t n,
                                const double* a, int64_t lda, double* r,
                                double* c, double* rowcnd, double* colcnd,
                                double* amax );
int64_t LAPACKE_cgeequ_work_64( int matrix_layout, int64_t m, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                float* r, float* c, float* rowcnd,
                                float* colcnd, float* amax );
int64_t LAPACKE_zgeequ_work_64( int matrix_layout, int64_t m, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                double* r, double* c, double* rowcnd,
                                double* colcnd, double* amax );

int64_t LAPACKE_sgeequb_work_64( int matrix_layout, int64_t m, int64_t n,
                                 const float* a, int64_t lda, float* r,
                                 float* c, float* rowcnd, float* colcnd,
                                 float* amax );
int64_t LAPACKE_dgeequb_work_64( int matrix_layout, int64_t m, int64_t n,
                                 const double* a, int64_t lda, double* r,
                                 double* c, double* rowcnd, double* colcnd,
                                 double* amax );
int64_t LAPACKE_cgeequb_work_64( int matrix_layout, int64_t m, int64_t n,
                                 const lapack_complex_float* a, int64_t lda,
                                 float* r, float* c, float* rowcnd,
                                 float* colcnd, float* amax );
int64_t LAPACKE_zgeequb_work_64( int matrix_layout, int64_t m, int64_t n,
                                 const lapack_complex_double* a, int64_t lda,
                                 double* r, double* c, double* rowcnd,
                                 double* colcnd, double* amax );

int64_t LAPACKE_sgees_work_64( int matrix_layout, char jobvs, char sort,
                               LAPACK_S_SELECT2 select, int64_t n, float* a,
                               int64_t lda, int64_t* sdim, float* wr,
                               float* wi, float* vs, int64_t ldvs,
                               float* work, int64_t lwork,
                               lapack_logical* bwork );
int64_t LAPACKE_dgees_work_64( int matrix_layout, char jobvs, char sort,
                               LAPACK_D_SELECT2 select, int64_t n, double* a,
                               int64_t lda, int64_t* sdim, double* wr,
                               double* wi, double* vs, int64_t ldvs,
                               double* work, int64_t lwork,
                               lapack_logical* bwork );
int64_t LAPACKE_cgees_work_64( int matrix_layout, char jobvs, char sort,
                               LAPACK_C_SELECT1 select, int64_t n,
                               lapack_complex_float* a, int64_t lda,
                               int64_t* sdim, lapack_complex_float* w,
                               lapack_complex_float* vs, int64_t ldvs,
                               lapack_complex_float* work, int64_t lwork,
                               float* rwork, lapack_logical* bwork );
int64_t LAPACKE_zgees_work_64( int matrix_layout, char jobvs, char sort,
                               LAPACK_Z_SELECT1 select, int64_t n,
                               lapack_complex_double* a, int64_t lda,
                               int64_t* sdim, lapack_complex_double* w,
                               lapack_complex_double* vs, int64_t ldvs,
                               lapack_complex_double* work, int64_t lwork,
                               double* rwork, lapack_logical* bwork );

int64_t LAPACKE_sgeesx_work_64( int matrix_layout, char jobvs, char sort,
                                LAPACK_S_SELECT2 select, char sense,
                                int64_t n, float* a, int64_t lda,
                                int64_t* sdim, float* wr, float* wi,
                                float* vs, int64_t ldvs, float* rconde,
                                float* rcondv, float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork,
                                lapack_logical* bwork );
int64_t LAPACKE_dgeesx_work_64( int matrix_layout, char jobvs, char sort,
                                LAPACK_D_SELECT2 select, char sense,
                                int64_t n, double* a, int64_t lda,
                                int64_t* sdim, double* wr, double* wi,
                                double* vs, int64_t ldvs, double* rconde,
                                double* rcondv, double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork,
                                lapack_logical* bwork );
int64_t LAPACKE_cgeesx_work_64( int matrix_layout, char jobvs, char sort,
                                LAPACK_C_SELECT1 select, char sense,
                                int64_t n, lapack_complex_float* a,
                                int64_t lda, int64_t* sdim,
                                lapack_complex_float* w,
                                lapack_complex_float* vs, int64_t ldvs,
                                float* rconde, float* rcondv,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, lapack_logical* bwork );
int64_t LAPACKE_zgeesx_work_64( int matrix_layout, char jobvs, char sort,
                                LAPACK_Z_SELECT1 select, char sense,
                                int64_t n, lapack_complex_double* a,
                                int64_t lda, int64_t* sdim,
                                lapack_complex_double* w,
                                lapack_complex_double* vs, int64_t ldvs,
                                double* rconde, double* rcondv,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, lapack_logical* bwork );

int64_t LAPACKE_sgeev_work_64( int matrix_layout, char jobvl, char jobvr,
                               int64_t n, float* a, int64_t lda,
                               float* wr, float* wi, float* vl, int64_t ldvl,
                               float* vr, int64_t ldvr, float* work,
                               int64_t lwork );
int64_t LAPACKE_dgeev_work_64( int matrix_layout, char jobvl, char jobvr,
                               int64_t n, double* a, int64_t lda,
                               double* wr, double* wi, double* vl,
                               int64_t ldvl, double* vr, int64_t ldvr,
                               double* work, int64_t lwork );
int64_t LAPACKE_cgeev_work_64( int matrix_layout, char jobvl, char jobvr,
                               int64_t n, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* w,
                               lapack_complex_float* vl, int64_t ldvl,
                               lapack_complex_float* vr, int64_t ldvr,
                               lapack_complex_float* work, int64_t lwork,
                               float* rwork );
int64_t LAPACKE_zgeev_work_64( int matrix_layout, char jobvl, char jobvr,
                               int64_t n, lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* w,
                               lapack_complex_double* vl, int64_t ldvl,
                               lapack_complex_double* vr, int64_t ldvr,
                               lapack_complex_double* work, int64_t lwork,
                               double* rwork );

int64_t LAPACKE_sgeevx_work_64( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, int64_t n, float* a,
                                int64_t lda, float* wr, float* wi, float* vl,
                                int64_t ldvl, float* vr, int64_t ldvr,
                                int64_t* ilo, int64_t* ihi, float* scale,
                                float* abnrm, float* rconde, float* rcondv,
                                float* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_dgeevx_work_64( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, int64_t n, double* a,
                                int64_t lda, double* wr, double* wi,
                                double* vl, int64_t ldvl, double* vr,
                                int64_t ldvr, int64_t* ilo,
                                int64_t* ihi, double* scale, double* abnrm,
                                double* rconde, double* rcondv, double* work,
                                int64_t lwork, int64_t* iwork );
int64_t LAPACKE_cgeevx_work_64( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* w,
                                lapack_complex_float* vl, int64_t ldvl,
                                lapack_complex_float* vr, int64_t ldvr,
                                int64_t* ilo, int64_t* ihi, float* scale,
                                float* abnrm, float* rconde, float* rcondv,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork );
int64_t LAPACKE_zgeevx_work_64( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* w,
                                lapack_complex_double* vl, int64_t ldvl,
                                lapack_complex_double* vr, int64_t ldvr,
                                int64_t* ilo, int64_t* ihi, double* scale,
                                double* abnrm, double* rconde, double* rcondv,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork );

int64_t LAPACKE_sgehrd_work_64( int matrix_layout, int64_t n, int64_t ilo,
                                int64_t ihi, float* a, int64_t lda,
                                float* tau, float* work, int64_t lwork );
int64_t LAPACKE_dgehrd_work_64( int matrix_layout, int64_t n, int64_t ilo,
                                int64_t ihi, double* a, int64_t lda,
                                double* tau, double* work, int64_t lwork );
int64_t LAPACKE_cgehrd_work_64( int matrix_layout, int64_t n, int64_t ilo,
                                int64_t ihi, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgehrd_work_64( int matrix_layout, int64_t n, int64_t ilo,
                                int64_t ihi, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgejsv_work_64( int matrix_layout, char joba, char jobu,
                                char jobv, char jobr, char jobt, char jobp,
                                int64_t m, int64_t n, float* a,
                                int64_t lda, float* sva, float* u,
                                int64_t ldu, float* v, int64_t ldv,
                                float* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_dgejsv_work_64( int matrix_layout, char joba, char jobu,
                                char jobv, char jobr, char jobt, char jobp,
                                int64_t m, int64_t n, double* a,
                                int64_t lda, double* sva, double* u,
                                int64_t ldu, double* v, int64_t ldv,
                                double* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_cgejsv_work_64( int matrix_layout, char joba, char jobu,
                                char jobv, char jobr, char jobt, char jobp,
                                int64_t m, int64_t n, lapack_complex_float* a,
                                int64_t lda, float* sva, lapack_complex_float* u,
                                int64_t ldu, lapack_complex_float* v, int64_t ldv,
                                lapack_complex_float* cwork, int64_t lwork,
                                float* work, int64_t lrwork,
                                int64_t* iwork );
int64_t LAPACKE_zgejsv_work_64( int matrix_layout, char joba, char jobu,
                                char jobv, char jobr, char jobt, char jobp,
                                int64_t m, int64_t n, lapack_complex_double* a,
                                int64_t lda, double* sva, lapack_complex_double* u,
                                int64_t ldu, lapack_complex_double* v, int64_t ldv,
                                lapack_complex_double* cwork, int64_t lwork,
                                double* work, int64_t lrwork,
                                int64_t* iwork );

int64_t LAPACKE_sgelq2_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, float* tau,
                                float* work );
int64_t LAPACKE_dgelq2_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, double* tau,
                                double* work );
int64_t LAPACKE_cgelq2_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* tau,
                                lapack_complex_float* work );
int64_t LAPACKE_zgelq2_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* tau,
                                lapack_complex_double* work );

int64_t LAPACKE_sgelqf_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, float* tau,
                                float* work, int64_t lwork );
int64_t LAPACKE_dgelqf_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, double* tau,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgelqf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgelqf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgels_work_64( int matrix_layout, char trans, int64_t m,
                               int64_t n, int64_t nrhs, float* a,
                               int64_t lda, float* b, int64_t ldb,
                               float* work, int64_t lwork );
int64_t LAPACKE_dgels_work_64( int matrix_layout, char trans, int64_t m,
                               int64_t n, int64_t nrhs, double* a,
                               int64_t lda, double* b, int64_t ldb,
                               double* work, int64_t lwork );
int64_t LAPACKE_cgels_work_64( int matrix_layout, char trans, int64_t m,
                               int64_t n, int64_t nrhs,
                               lapack_complex_float* a, int64_t lda,
                               lapack_complex_float* b, int64_t ldb,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgels_work_64( int matrix_layout, char trans, int64_t m,
                               int64_t n, int64_t nrhs,
                               lapack_complex_double* a, int64_t lda,
                               lapack_complex_double* b, int64_t ldb,
                               lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgelsd_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, float* a, int64_t lda,
                                float* b, int64_t ldb, float* s, float rcond,
                                int64_t* rank, float* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_dgelsd_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, double* a, int64_t lda,
                                double* b, int64_t ldb, double* s,
                                double rcond, int64_t* rank, double* work,
                                int64_t lwork, int64_t* iwork );
int64_t LAPACKE_cgelsd_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* b,
                                int64_t ldb, float* s, float rcond,
                                int64_t* rank, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t* iwork );
int64_t LAPACKE_zgelsd_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb, double* s, double rcond,
                                int64_t* rank, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t* iwork );

int64_t LAPACKE_sgelss_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, float* a, int64_t lda,
                                float* b, int64_t ldb, float* s, float rcond,
                                int64_t* rank, float* work,
                                int64_t lwork );
int64_t LAPACKE_dgelss_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, double* a, int64_t lda,
                                double* b, int64_t ldb, double* s,
                                double rcond, int64_t* rank, double* work,
                                int64_t lwork );
int64_t LAPACKE_cgelss_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* b,
                                int64_t ldb, float* s, float rcond,
                                int64_t* rank, lapack_complex_float* work,
                                int64_t lwork, float* rwork );
int64_t LAPACKE_zgelss_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb, double* s, double rcond,
                                int64_t* rank, lapack_complex_double* work,
                                int64_t lwork, double* rwork );

int64_t LAPACKE_sgelsy_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, float* a, int64_t lda,
                                float* b, int64_t ldb, int64_t* jpvt,
                                float rcond, int64_t* rank, float* work,
                                int64_t lwork );
int64_t LAPACKE_dgelsy_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, double* a, int64_t lda,
                                double* b, int64_t ldb, int64_t* jpvt,
                                double rcond, int64_t* rank, double* work,
                                int64_t lwork );
int64_t LAPACKE_cgelsy_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* b,
                                int64_t ldb, int64_t* jpvt, float rcond,
                                int64_t* rank, lapack_complex_float* work,
                                int64_t lwork, float* rwork );
int64_t LAPACKE_zgelsy_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nrhs, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb, int64_t* jpvt, double rcond,
                                int64_t* rank, lapack_complex_double* work,
                                int64_t lwork, double* rwork );

int64_t LAPACKE_sgeqlf_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, float* tau,
                                float* work, int64_t lwork );
int64_t LAPACKE_dgeqlf_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, double* tau,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgeqlf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgeqlf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgeqp3_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, int64_t* jpvt,
                                float* tau, float* work, int64_t lwork );
int64_t LAPACKE_dgeqp3_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, int64_t* jpvt,
                                double* tau, double* work, int64_t lwork );
int64_t LAPACKE_cgeqp3_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* jpvt, lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork );
int64_t LAPACKE_zgeqp3_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* jpvt, lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork );

int64_t LAPACKE_sgeqpf_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, int64_t* jpvt,
                                float* tau, float* work );
int64_t LAPACKE_dgeqpf_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, int64_t* jpvt,
                                double* tau, double* work );
int64_t LAPACKE_cgeqpf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* jpvt, lapack_complex_float* tau,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zgeqpf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* jpvt, lapack_complex_double* tau,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sgeqr2_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, float* tau,
                                float* work );
int64_t LAPACKE_dgeqr2_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, double* tau,
                                double* work );
int64_t LAPACKE_cgeqr2_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* tau,
                                lapack_complex_float* work );
int64_t LAPACKE_zgeqr2_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* tau,
                                lapack_complex_double* work );

int64_t LAPACKE_sgeqrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, float* tau,
                                float* work, int64_t lwork );
int64_t LAPACKE_dgeqrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, double* tau,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgeqrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgeqrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgeqrfp_work_64( int matrix_layout, int64_t m, int64_t n,
                                 float* a, int64_t lda, float* tau,
                                 float* work, int64_t lwork );
int64_t LAPACKE_dgeqrfp_work_64( int matrix_layout, int64_t m, int64_t n,
                                 double* a, int64_t lda, double* tau,
                                 double* work, int64_t lwork );
int64_t LAPACKE_cgeqrfp_work_64( int matrix_layout, int64_t m, int64_t n,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* tau,
                                 lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgeqrfp_work_64( int matrix_layout, int64_t m, int64_t n,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* tau,
                                 lapack_complex_double* work,
                                 int64_t lwork );

int64_t LAPACKE_sgerfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const float* af, int64_t ldaf,
                                const int64_t* ipiv, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* ferr, float* berr, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dgerfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const double* a,
                                int64_t lda, const double* af,
                                int64_t ldaf, const int64_t* ipiv,
                                const double* b, int64_t ldb, double* x,
                                int64_t ldx, double* ferr, double* berr,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cgerfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* af,
                                int64_t ldaf, const int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zgerfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, const lapack_complex_double* af,
                                int64_t ldaf, const int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sgerfsx_work_64( int matrix_layout, char trans, char equed,
                                 int64_t n, int64_t nrhs, const float* a,
                                 int64_t lda, const float* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const float* r, const float* c, const float* b,
                                 int64_t ldb, float* x, int64_t ldx,
                                 float* rcond, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, float* work,
                                 int64_t* iwork );
int64_t LAPACKE_dgerfsx_work_64( int matrix_layout, char trans, char equed,
                                 int64_t n, int64_t nrhs, const double* a,
                                 int64_t lda, const double* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const double* r, const double* c,
                                 const double* b, int64_t ldb, double* x,
                                 int64_t ldx, double* rcond, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, double* work,
                                 int64_t* iwork );
int64_t LAPACKE_cgerfsx_work_64( int matrix_layout, char trans, char equed,
                                 int64_t n, int64_t nrhs,
                                 const lapack_complex_float* a, int64_t lda,
                                 const lapack_complex_float* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const float* r, const float* c,
                                 const lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* x, int64_t ldx,
                                 float* rcond, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, lapack_complex_float* work,
                                 float* rwork );
int64_t LAPACKE_zgerfsx_work_64( int matrix_layout, char trans, char equed,
                                 int64_t n, int64_t nrhs,
                                 const lapack_complex_double* a, int64_t lda,
                                 const lapack_complex_double* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const double* r, const double* c,
                                 const lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* x, int64_t ldx,
                                 double* rcond, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, lapack_complex_double* work,
                                 double* rwork );

int64_t LAPACKE_sgerqf_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, float* tau,
                                float* work, int64_t lwork );
int64_t LAPACKE_dgerqf_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, double* tau,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgerqf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgerqf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgesdd_work_64( int matrix_layout, char jobz, int64_t m,
                                int64_t n, float* a, int64_t lda,
                                float* s, float* u, int64_t ldu, float* vt,
                                int64_t ldvt, float* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_dgesdd_work_64( int matrix_layout, char jobz, int64_t m,
                                int64_t n, double* a, int64_t lda,
                                double* s, double* u, int64_t ldu,
                                double* vt, int64_t ldvt, double* work,
                                int64_t lwork, int64_t* iwork );
int64_t LAPACKE_cgesdd_work_64( int matrix_layout, char jobz, int64_t m,
                                int64_t n, lapack_complex_float* a,
                                int64_t lda, float* s,
                                lapack_complex_float* u, int64_t ldu,
                                lapack_complex_float* vt, int64_t ldvt,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, int64_t* iwork );
int64_t LAPACKE_zgesdd_work_64( int matrix_layout, char jobz, int64_t m,
                                int64_t n, lapack_complex_double* a,
                                int64_t lda, double* s,
                                lapack_complex_double* u, int64_t ldu,
                                lapack_complex_double* vt, int64_t ldvt,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, int64_t* iwork );

int64_t LAPACKE_sgedmd_work_64( int matrix_layout, char jobs, char jobz,
                                char jobr, char jobf, int64_t whtsvd,
                                int64_t m, int64_t n, float* x,
                                int64_t ldx, float* y, int64_t ldy,
                                int64_t nrnk, float* tol, int64_t k,
                                float* reig, float* imeig,
                                float* z, int64_t ldz, float* res,
                                float* b, int64_t ldb, float* w,
                                int64_t ldw, float* s, int64_t lds,
                                float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_dgedmd_work_64( int matrix_layout, char jobs, char jobz,
                                char jobr, char jobf, int64_t whtsvd,
                                int64_t m, int64_t n, double* x,
                                int64_t ldx, double* y, int64_t ldy,
                                int64_t nrnk, double* tol, int64_t k,
                                double* reig, double *imeig,
                                double* z, int64_t ldz, double* res,
                                double* b, int64_t ldb, double* w,
                                int64_t ldw, double* s, int64_t lds,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_cgedmd_work_64( int matrix_layout, char jobs, char jobz,
                                char jobr, char jobf, int64_t whtsvd,
                                int64_t m, int64_t n,
                                lapack_complex_float* x, int64_t ldx,
                                lapack_complex_float* y, int64_t ldy,
                                int64_t nrnk, float* tol, int64_t k,
                                lapack_complex_float* eigs,
                                lapack_complex_float* z, int64_t ldz,
                                float* res,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* w, int64_t ldw,
                                lapack_complex_float* s, int64_t lds,
                                lapack_complex_float* zwork, int64_t lzwork,
                                float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_zgedmd_work_64( int matrix_layout, char jobs, char jobz,
                                char jobr, char jobf, int64_t whtsvd,
                                int64_t m, int64_t n,
                                lapack_complex_double* x, int64_t ldx,
                                lapack_complex_double* y, int64_t ldy,
                                int64_t nrnk, double* tol, int64_t k, 
                                lapack_complex_double* eigs,
                                lapack_complex_double* z, int64_t ldz,
                                double* res,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* w, int64_t ldw,
                                lapack_complex_double* s, int64_t lds,
                                lapack_complex_double* zwork, int64_t lzwork,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_sgedmdq_work_64( int matrix_layout, char jobs, char jobz,
                                 char jobr, char jobq, char jobt, char jobf,
                                 int64_t whtsvd, int64_t m, int64_t n,
                                 float* f, int64_t ldf, float* x,
                                 int64_t ldx, float* y, int64_t ldy,
                                 int64_t nrnk, float* tol, int64_t k,
                                 float* reig, float *imeig, float* z,
                                 int64_t ldz, float* res, float* b,
                                 int64_t ldb, float* v, int64_t ldv,
                                 float* s, int64_t lds, float* work,
                                 int64_t lwork, int64_t* iwork,
                                 int64_t liwork );

int64_t LAPACKE_dgedmdq_work_64( int matrix_layout, char jobs, char jobz,
                                 char jobr, char jobq, char jobt, char jobf,
                                 int64_t whtsvd, int64_t m, int64_t n,
                                 double* f, int64_t ldf, double* x,
                                 int64_t ldx, double* y, int64_t ldy,
                                 int64_t nrnk, double* tol, int64_t k,
                                 double* reig, double* imeig, double* z,
                                 int64_t ldz, double* res, double* b,
                                 int64_t ldb, double* v, int64_t ldv,
                                 double* s, int64_t lds, double* work,
                                 int64_t lwork, int64_t* iwork,
                                 int64_t liwork );

int64_t LAPACKE_cgedmdq_work_64( int matrix_layout, char jobs, char jobz,
                                 char jobr, char jobq, char jobt, char jobf,
                                 int64_t whtsvd, int64_t m, int64_t n,
                                 lapack_complex_float* f, int64_t ldf,
                                 lapack_complex_float* x, int64_t ldx,
                                 lapack_complex_float* y, int64_t ldy,
                                 int64_t nrnk, float* tol, int64_t k,
                                 lapack_complex_float* eigs,
                                 lapack_complex_float* z, int64_t ldz,
                                 float* res,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* v, int64_t ldv,
                                 lapack_complex_float* s, int64_t lds,
                                 lapack_complex_float* zwork, int64_t lzwork,
                                 float* work, int64_t lwork,
                                 int64_t* iwork, int64_t liwork);

int64_t LAPACKE_zgedmdq_work_64( int matrix_layout, char jobs, char jobz,
                                 char jobr, char jobq, char jobt, char jobf,
                                 int64_t whtsvd, int64_t m, int64_t n,
                                 lapack_complex_double* f, int64_t ldf,
                                 lapack_complex_double* x, int64_t ldx,
                                 lapack_complex_double* y, int64_t ldy,
                                 int64_t nrnk, double* tol, int64_t k,
                                 lapack_complex_double* eigs,
                                 lapack_complex_double* z, int64_t ldz,
                                 double* res,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* v, int64_t ldv,
                                 lapack_complex_double* s, int64_t lds,
                                 lapack_complex_double* zwork, int64_t lzwork,
                                 double* work, int64_t lwork,
                                 int64_t* iwork, int64_t liwork);


int64_t LAPACKE_sgesv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               float* a, int64_t lda, int64_t* ipiv,
                               float* b, int64_t ldb );
int64_t LAPACKE_dgesv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               double* a, int64_t lda, int64_t* ipiv,
                               double* b, int64_t ldb );
int64_t LAPACKE_cgesv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               lapack_complex_float* a, int64_t lda,
                               int64_t* ipiv, lapack_complex_float* b,
                               int64_t ldb );
int64_t LAPACKE_zgesv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               lapack_complex_double* a, int64_t lda,
                               int64_t* ipiv, lapack_complex_double* b,
                               int64_t ldb );
int64_t LAPACKE_dsgesv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                                double* a, int64_t lda, int64_t* ipiv,
                                double* b, int64_t ldb, double* x,
                                int64_t ldx, double* work, float* swork,
                                int64_t* iter );
int64_t LAPACKE_zcgesv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_double* b,
                                int64_t ldb, lapack_complex_double* x,
                                int64_t ldx, lapack_complex_double* work,
                                lapack_complex_float* swork, double* rwork,
                                int64_t* iter );

int64_t LAPACKE_sgesvd_work_64( int matrix_layout, char jobu, char jobvt,
                                int64_t m, int64_t n, float* a,
                                int64_t lda, float* s, float* u,
                                int64_t ldu, float* vt, int64_t ldvt,
                                float* work, int64_t lwork );
int64_t LAPACKE_dgesvd_work_64( int matrix_layout, char jobu, char jobvt,
                                int64_t m, int64_t n, double* a,
                                int64_t lda, double* s, double* u,
                                int64_t ldu, double* vt, int64_t ldvt,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgesvd_work_64( int matrix_layout, char jobu, char jobvt,
                                int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                float* s, lapack_complex_float* u,
                                int64_t ldu, lapack_complex_float* vt,
                                int64_t ldvt, lapack_complex_float* work,
                                int64_t lwork, float* rwork );
int64_t LAPACKE_zgesvd_work_64( int matrix_layout, char jobu, char jobvt,
                                int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                double* s, lapack_complex_double* u,
                                int64_t ldu, lapack_complex_double* vt,
                                int64_t ldvt, lapack_complex_double* work,
                                int64_t lwork, double* rwork );

int64_t LAPACKE_sgesvdx_work_64( int matrix_layout, char jobu, char jobvt, char range,
                                 int64_t m, int64_t n, float* a,
                                 int64_t lda, float vl, float vu,
                                 int64_t il, int64_t iu, int64_t* ns,
                                 float* s, float* u, int64_t ldu,
                                 float* vt, int64_t ldvt,
                                 float* work, int64_t lwork, int64_t* iwork );
int64_t LAPACKE_dgesvdx_work_64( int matrix_layout, char jobu, char jobvt, char range,
                                 int64_t m, int64_t n, double* a,
                                 int64_t lda, double vl, double vu,
                                 int64_t il, int64_t iu, int64_t* ns,
                                 double* s, double* u, int64_t ldu,
                                 double* vt, int64_t ldvt,
                                 double* work, int64_t lwork, int64_t* iwork );
int64_t LAPACKE_cgesvdx_work_64( int matrix_layout, char jobu, char jobvt, char range,
                                 int64_t m, int64_t n, lapack_complex_float* a,
                                 int64_t lda, float vl, float vu,
                                 int64_t il, int64_t iu, int64_t* ns,
                                 float* s, lapack_complex_float* u, int64_t ldu,
                                 lapack_complex_float* vt, int64_t ldvt,
                                 lapack_complex_float* work, int64_t lwork,
                                 float* rwork, int64_t* iwork );
int64_t LAPACKE_zgesvdx_work_64( int matrix_layout, char jobu, char jobvt, char range,
                                 int64_t m, int64_t n, lapack_complex_double* a,
                                 int64_t lda, double vl, double vu,
                                 int64_t il, int64_t iu, int64_t* ns,
                                 double* s, lapack_complex_double* u, int64_t ldu,
                                 lapack_complex_double* vt, int64_t ldvt,
                                 lapack_complex_double* work, int64_t lwork,
                                 double* rwork, int64_t* iwork );

int64_t LAPACKE_sgesvdq_work_64( int matrix_layout, char joba, char jobp,
                                char jobr, char jobu, char jobv,
                                int64_t m, int64_t n, float* a,
                                int64_t lda, float* s, float* u,
                                int64_t ldu, float* v, int64_t ldv,
                                int64_t* numrank,
                                int64_t* iwork, int64_t liwork,
                                float* work, int64_t lwork,
                                float* rwork, int64_t lrwork);
int64_t LAPACKE_dgesvdq_work_64( int matrix_layout, char joba, char jobp,
                                char jobr, char jobu, char jobv,
                                int64_t m, int64_t n, double* a,
                                int64_t lda, double* s, double* u,
                                int64_t ldu, double* v, int64_t ldv,
                                int64_t* numrank,
                                int64_t* iwork, int64_t liwork,
                                double* work, int64_t lwork,
                                double* rwork, int64_t lrwork);
int64_t LAPACKE_cgesvdq_work_64( int matrix_layout, char joba, char jobp,
                                char jobr, char jobu, char jobv,
                                int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                float* s, lapack_complex_float* u,
                                int64_t ldu, lapack_complex_float* v,
                                int64_t ldv, int64_t* numrank,
                                int64_t* iwork, int64_t liwork,
                                lapack_complex_float* cwork, int64_t lcwork,
                                float* rwork, int64_t lrwork);
int64_t LAPACKE_zgesvdq_work_64( int matrix_layout, char joba, char jobp,
                                char jobr, char jobu, char jobv,
                                int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                double* s, lapack_complex_double* u,
                                int64_t ldu, lapack_complex_double* v,
                                int64_t ldv, int64_t* numrank,
                                int64_t* iwork, int64_t liwork,
                                lapack_complex_double* cwork, int64_t lcwork,
                                double* rwork, int64_t lrwork);

int64_t LAPACKE_sgesvj_work_64( int matrix_layout, char joba, char jobu,
                                char jobv, int64_t m, int64_t n, float* a,
                                int64_t lda, float* sva, int64_t mv,
                                float* v, int64_t ldv, float* work,
                                int64_t lwork );
int64_t LAPACKE_dgesvj_work_64( int matrix_layout, char joba, char jobu,
                                char jobv, int64_t m, int64_t n,
                                double* a, int64_t lda, double* sva,
                                int64_t mv, double* v, int64_t ldv,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgesvj_work_64( int matrix_layout, char joba, char jobu,
                                char jobv, int64_t m, int64_t n, lapack_complex_float* a,
                                int64_t lda, float* sva, int64_t mv,
                                lapack_complex_float* v, int64_t ldv,
                                lapack_complex_float* cwork, int64_t lwork,
                                float* rwork,int64_t lrwork );
int64_t LAPACKE_zgesvj_work_64( int matrix_layout, char joba, char jobu,
                                char jobv, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda, double* sva,
                                int64_t mv, lapack_complex_double* v, int64_t ldv,
                                lapack_complex_double* cwork, int64_t lwork,
                                double* rwork, int64_t lrwork );

int64_t LAPACKE_sgesvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t nrhs, float* a,
                                int64_t lda, float* af, int64_t ldaf,
                                int64_t* ipiv, char* equed, float* r,
                                float* c, float* b, int64_t ldb, float* x,
                                int64_t ldx, float* rcond, float* ferr,
                                float* berr, float* work, int64_t* iwork );
int64_t LAPACKE_dgesvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t nrhs, double* a,
                                int64_t lda, double* af, int64_t ldaf,
                                int64_t* ipiv, char* equed, double* r,
                                double* c, double* b, int64_t ldb, double* x,
                                int64_t ldx, double* rcond, double* ferr,
                                double* berr, double* work, int64_t* iwork );
int64_t LAPACKE_cgesvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t nrhs,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* af, int64_t ldaf,
                                int64_t* ipiv, char* equed, float* r,
                                float* c, lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* x,
                                int64_t ldx, float* rcond, float* ferr,
                                float* berr, lapack_complex_float* work,
                                float* rwork );
int64_t LAPACKE_zgesvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t nrhs,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* af, int64_t ldaf,
                                int64_t* ipiv, char* equed, double* r,
                                double* c, lapack_complex_double* b,
                                int64_t ldb, lapack_complex_double* x,
                                int64_t ldx, double* rcond, double* ferr,
                                double* berr, lapack_complex_double* work,
                                double* rwork );

int64_t LAPACKE_sgesvxx_work_64( int matrix_layout, char fact, char trans,
                                 int64_t n, int64_t nrhs, float* a,
                                 int64_t lda, float* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, float* r,
                                 float* c, float* b, int64_t ldb, float* x,
                                 int64_t ldx, float* rcond, float* rpvgrw,
                                 float* berr, int64_t n_err_bnds,
                                 float* err_bnds_norm, float* err_bnds_comp,
                                 int64_t nparams, float* params, float* work,
                                 int64_t* iwork );
int64_t LAPACKE_dgesvxx_work_64( int matrix_layout, char fact, char trans,
                                 int64_t n, int64_t nrhs, double* a,
                                 int64_t lda, double* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, double* r,
                                 double* c, double* b, int64_t ldb,
                                 double* x, int64_t ldx, double* rcond,
                                 double* rpvgrw, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, double* work,
                                 int64_t* iwork );
int64_t LAPACKE_cgesvxx_work_64( int matrix_layout, char fact, char trans,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, float* r,
                                 float* c, lapack_complex_float* b,
                                 int64_t ldb, lapack_complex_float* x,
                                 int64_t ldx, float* rcond, float* rpvgrw,
                                 float* berr, int64_t n_err_bnds,
                                 float* err_bnds_norm, float* err_bnds_comp,
                                 int64_t nparams, float* params,
                                 lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zgesvxx_work_64( int matrix_layout, char fact, char trans,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, double* r,
                                 double* c, lapack_complex_double* b,
                                 int64_t ldb, lapack_complex_double* x,
                                 int64_t ldx, double* rcond, double* rpvgrw,
                                 double* berr, int64_t n_err_bnds,
                                 double* err_bnds_norm, double* err_bnds_comp,
                                 int64_t nparams, double* params,
                                 lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sgetf2_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dgetf2_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_cgetf2_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ipiv );
int64_t LAPACKE_zgetf2_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv );

int64_t LAPACKE_sgetrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dgetrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_cgetrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ipiv );
int64_t LAPACKE_zgetrf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv );

int64_t LAPACKE_sgetrf2_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dgetrf2_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, int64_t* ipiv );
int64_t LAPACKE_cgetrf2_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ipiv );
int64_t LAPACKE_zgetrf2_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv );

int64_t LAPACKE_sgetri_work_64( int matrix_layout, int64_t n, float* a,
                                int64_t lda, const int64_t* ipiv,
                                float* work, int64_t lwork );
int64_t LAPACKE_dgetri_work_64( int matrix_layout, int64_t n, double* a,
                                int64_t lda, const int64_t* ipiv,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgetri_work_64( int matrix_layout, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                const int64_t* ipiv,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgetri_work_64( int matrix_layout, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                const int64_t* ipiv,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgetrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const int64_t* ipiv, float* b,
                                int64_t ldb );
int64_t LAPACKE_dgetrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const double* a,
                                int64_t lda, const int64_t* ipiv,
                                double* b, int64_t ldb );
int64_t LAPACKE_cgetrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zgetrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sggbak_work_64( int matrix_layout, char job, char side,
                                int64_t n, int64_t ilo, int64_t ihi,
                                const float* lscale, const float* rscale,
                                int64_t m, float* v, int64_t ldv );
int64_t LAPACKE_dggbak_work_64( int matrix_layout, char job, char side,
                                int64_t n, int64_t ilo, int64_t ihi,
                                const double* lscale, const double* rscale,
                                int64_t m, double* v, int64_t ldv );
int64_t LAPACKE_cggbak_work_64( int matrix_layout, char job, char side,
                                int64_t n, int64_t ilo, int64_t ihi,
                                const float* lscale, const float* rscale,
                                int64_t m, lapack_complex_float* v,
                                int64_t ldv );
int64_t LAPACKE_zggbak_work_64( int matrix_layout, char job, char side,
                                int64_t n, int64_t ilo, int64_t ihi,
                                const double* lscale, const double* rscale,
                                int64_t m, lapack_complex_double* v,
                                int64_t ldv );

int64_t LAPACKE_sggbal_work_64( int matrix_layout, char job, int64_t n,
                                float* a, int64_t lda, float* b,
                                int64_t ldb, int64_t* ilo,
                                int64_t* ihi, float* lscale, float* rscale,
                                float* work );
int64_t LAPACKE_dggbal_work_64( int matrix_layout, char job, int64_t n,
                                double* a, int64_t lda, double* b,
                                int64_t ldb, int64_t* ilo,
                                int64_t* ihi, double* lscale, double* rscale,
                                double* work );
int64_t LAPACKE_cggbal_work_64( int matrix_layout, char job, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                int64_t* ilo, int64_t* ihi, float* lscale,
                                float* rscale, float* work );
int64_t LAPACKE_zggbal_work_64( int matrix_layout, char job, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                int64_t* ilo, int64_t* ihi,
                                double* lscale, double* rscale, double* work );

int64_t LAPACKE_sgges_work_64( int matrix_layout, char jobvsl, char jobvsr,
                               char sort, LAPACK_S_SELECT3 selctg, int64_t n,
                               float* a, int64_t lda, float* b,
                               int64_t ldb, int64_t* sdim, float* alphar,
                               float* alphai, float* beta, float* vsl,
                               int64_t ldvsl, float* vsr, int64_t ldvsr,
                               float* work, int64_t lwork,
                               lapack_logical* bwork );
int64_t LAPACKE_dgges_work_64( int matrix_layout, char jobvsl, char jobvsr,
                               char sort, LAPACK_D_SELECT3 selctg, int64_t n,
                               double* a, int64_t lda, double* b,
                               int64_t ldb, int64_t* sdim, double* alphar,
                               double* alphai, double* beta, double* vsl,
                               int64_t ldvsl, double* vsr, int64_t ldvsr,
                               double* work, int64_t lwork,
                               lapack_logical* bwork );
int64_t LAPACKE_cgges_work_64( int matrix_layout, char jobvsl, char jobvsr,
                               char sort, LAPACK_C_SELECT2 selctg, int64_t n,
                               lapack_complex_float* a, int64_t lda,
                               lapack_complex_float* b, int64_t ldb,
                               int64_t* sdim, lapack_complex_float* alpha,
                               lapack_complex_float* beta,
                               lapack_complex_float* vsl, int64_t ldvsl,
                               lapack_complex_float* vsr, int64_t ldvsr,
                               lapack_complex_float* work, int64_t lwork,
                               float* rwork, lapack_logical* bwork );
int64_t LAPACKE_zgges_work_64( int matrix_layout, char jobvsl, char jobvsr,
                               char sort, LAPACK_Z_SELECT2 selctg, int64_t n,
                               lapack_complex_double* a, int64_t lda,
                               lapack_complex_double* b, int64_t ldb,
                               int64_t* sdim, lapack_complex_double* alpha,
                               lapack_complex_double* beta,
                               lapack_complex_double* vsl, int64_t ldvsl,
                               lapack_complex_double* vsr, int64_t ldvsr,
                               lapack_complex_double* work, int64_t lwork,
                               double* rwork, lapack_logical* bwork );

int64_t LAPACKE_sgges3_work_64( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_S_SELECT3 selctg,
                                int64_t n,
                                float* a, int64_t lda,
                                float* b, int64_t ldb, int64_t* sdim,
                                float* alphar, float* alphai, float* beta,
                                float* vsl, int64_t ldvsl,
                                float* vsr, int64_t ldvsr,
                                float* work, int64_t lwork,
                                lapack_logical* bwork );
int64_t LAPACKE_dgges3_work_64( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_D_SELECT3 selctg,
                                int64_t n,
                                double* a, int64_t lda,
                                double* b, int64_t ldb, int64_t* sdim,
                                double* alphar, double* alphai, double* beta,
                                double* vsl, int64_t ldvsl,
                                double* vsr, int64_t ldvsr,
                                double* work, int64_t lwork,
                                lapack_logical* bwork );
int64_t LAPACKE_cgges3_work_64( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_C_SELECT2 selctg,
                                int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                int64_t* sdim, lapack_complex_float* alpha,
                                lapack_complex_float* beta,
                                lapack_complex_float* vsl, int64_t ldvsl,
                                lapack_complex_float* vsr, int64_t ldvsr,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, lapack_logical* bwork );
int64_t LAPACKE_zgges3_work_64( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_Z_SELECT2 selctg,
                                int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                int64_t* sdim, lapack_complex_double* alpha,
                                lapack_complex_double* beta,
                                lapack_complex_double* vsl, int64_t ldvsl,
                                lapack_complex_double* vsr, int64_t ldvsr,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, lapack_logical* bwork );

int64_t LAPACKE_sggesx_work_64( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_S_SELECT3 selctg, char sense,
                                int64_t n, float* a, int64_t lda,
                                float* b, int64_t ldb, int64_t* sdim,
                                float* alphar, float* alphai, float* beta,
                                float* vsl, int64_t ldvsl, float* vsr,
                                int64_t ldvsr, float* rconde, float* rcondv,
                                float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork,
                                lapack_logical* bwork );
int64_t LAPACKE_dggesx_work_64( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_D_SELECT3 selctg, char sense,
                                int64_t n, double* a, int64_t lda,
                                double* b, int64_t ldb, int64_t* sdim,
                                double* alphar, double* alphai, double* beta,
                                double* vsl, int64_t ldvsl, double* vsr,
                                int64_t ldvsr, double* rconde,
                                double* rcondv, double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork,
                                lapack_logical* bwork );
int64_t LAPACKE_cggesx_work_64( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_C_SELECT2 selctg, char sense,
                                int64_t n, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* b,
                                int64_t ldb, int64_t* sdim,
                                lapack_complex_float* alpha,
                                lapack_complex_float* beta,
                                lapack_complex_float* vsl, int64_t ldvsl,
                                lapack_complex_float* vsr, int64_t ldvsr,
                                float* rconde, float* rcondv,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, int64_t* iwork,
                                int64_t liwork, lapack_logical* bwork );
int64_t LAPACKE_zggesx_work_64( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_Z_SELECT2 selctg, char sense,
                                int64_t n, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb, int64_t* sdim,
                                lapack_complex_double* alpha,
                                lapack_complex_double* beta,
                                lapack_complex_double* vsl, int64_t ldvsl,
                                lapack_complex_double* vsr, int64_t ldvsr,
                                double* rconde, double* rcondv,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, int64_t* iwork,
                                int64_t liwork, lapack_logical* bwork );

int64_t LAPACKE_sggev_work_64( int matrix_layout, char jobvl, char jobvr,
                               int64_t n, float* a, int64_t lda, float* b,
                               int64_t ldb, float* alphar, float* alphai,
                               float* beta, float* vl, int64_t ldvl,
                               float* vr, int64_t ldvr, float* work,
                               int64_t lwork );
int64_t LAPACKE_dggev_work_64( int matrix_layout, char jobvl, char jobvr,
                               int64_t n, double* a, int64_t lda,
                               double* b, int64_t ldb, double* alphar,
                               double* alphai, double* beta, double* vl,
                               int64_t ldvl, double* vr, int64_t ldvr,
                               double* work, int64_t lwork );
int64_t LAPACKE_cggev_work_64( int matrix_layout, char jobvl, char jobvr,
                               int64_t n, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* b,
                               int64_t ldb, lapack_complex_float* alpha,
                               lapack_complex_float* beta,
                               lapack_complex_float* vl, int64_t ldvl,
                               lapack_complex_float* vr, int64_t ldvr,
                               lapack_complex_float* work, int64_t lwork,
                               float* rwork );
int64_t LAPACKE_zggev_work_64( int matrix_layout, char jobvl, char jobvr,
                               int64_t n, lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* b,
                               int64_t ldb, lapack_complex_double* alpha,
                               lapack_complex_double* beta,
                               lapack_complex_double* vl, int64_t ldvl,
                               lapack_complex_double* vr, int64_t ldvr,
                               lapack_complex_double* work, int64_t lwork,
                               double* rwork );

int64_t LAPACKE_sggev3_work_64( int matrix_layout, char jobvl, char jobvr,
                                int64_t n,
                                float* a, int64_t lda,
                                float* b, int64_t ldb,
                                float* alphar, float* alphai, float* beta,
                                float* vl, int64_t ldvl,
                                float* vr, int64_t ldvr,
                                float* work, int64_t lwork );
int64_t LAPACKE_dggev3_work_64( int matrix_layout, char jobvl, char jobvr,
                                int64_t n,
                                double* a, int64_t lda,
                                double* b, int64_t ldb,
                                double* alphar, double* alphai, double* beta,
                                double* vl, int64_t ldvl,
                                double* vr, int64_t ldvr,
                                double* work, int64_t lwork );
int64_t LAPACKE_cggev3_work_64( int matrix_layout, char jobvl, char jobvr,
                                int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* alpha,
                                lapack_complex_float* beta,
                                lapack_complex_float* vl, int64_t ldvl,
                                lapack_complex_float* vr, int64_t ldvr,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork );
int64_t LAPACKE_zggev3_work_64( int matrix_layout, char jobvl, char jobvr,
                                int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* alpha,
                                lapack_complex_double* beta,
                                lapack_complex_double* vl, int64_t ldvl,
                                lapack_complex_double* vr, int64_t ldvr,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork );

int64_t LAPACKE_sggevx_work_64( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, int64_t n, float* a,
                                int64_t lda, float* b, int64_t ldb,
                                float* alphar, float* alphai, float* beta,
                                float* vl, int64_t ldvl, float* vr,
                                int64_t ldvr, int64_t* ilo,
                                int64_t* ihi, float* lscale, float* rscale,
                                float* abnrm, float* bbnrm, float* rconde,
                                float* rcondv, float* work, int64_t lwork,
                                int64_t* iwork, lapack_logical* bwork );
int64_t LAPACKE_dggevx_work_64( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, int64_t n, double* a,
                                int64_t lda, double* b, int64_t ldb,
                                double* alphar, double* alphai, double* beta,
                                double* vl, int64_t ldvl, double* vr,
                                int64_t ldvr, int64_t* ilo,
                                int64_t* ihi, double* lscale, double* rscale,
                                double* abnrm, double* bbnrm, double* rconde,
                                double* rcondv, double* work, int64_t lwork,
                                int64_t* iwork, lapack_logical* bwork );
int64_t LAPACKE_cggevx_work_64( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* alpha,
                                lapack_complex_float* beta,
                                lapack_complex_float* vl, int64_t ldvl,
                                lapack_complex_float* vr, int64_t ldvr,
                                int64_t* ilo, int64_t* ihi, float* lscale,
                                float* rscale, float* abnrm, float* bbnrm,
                                float* rconde, float* rcondv,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, int64_t* iwork,
                                lapack_logical* bwork );
int64_t LAPACKE_zggevx_work_64( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* alpha,
                                lapack_complex_double* beta,
                                lapack_complex_double* vl, int64_t ldvl,
                                lapack_complex_double* vr, int64_t ldvr,
                                int64_t* ilo, int64_t* ihi,
                                double* lscale, double* rscale, double* abnrm,
                                double* bbnrm, double* rconde, double* rcondv,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, int64_t* iwork,
                                lapack_logical* bwork );

int64_t LAPACKE_sggglm_work_64( int matrix_layout, int64_t n, int64_t m,
                                int64_t p, float* a, int64_t lda,
                                float* b, int64_t ldb, float* d, float* x,
                                float* y, float* work, int64_t lwork );
int64_t LAPACKE_dggglm_work_64( int matrix_layout, int64_t n, int64_t m,
                                int64_t p, double* a, int64_t lda,
                                double* b, int64_t ldb, double* d, double* x,
                                double* y, double* work, int64_t lwork );
int64_t LAPACKE_cggglm_work_64( int matrix_layout, int64_t n, int64_t m,
                                int64_t p, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* d,
                                lapack_complex_float* x,
                                lapack_complex_float* y,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zggglm_work_64( int matrix_layout, int64_t n, int64_t m,
                                int64_t p, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb, lapack_complex_double* d,
                                lapack_complex_double* x,
                                lapack_complex_double* y,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgghrd_work_64( int matrix_layout, char compq, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                float* a, int64_t lda, float* b,
                                int64_t ldb, float* q, int64_t ldq,
                                float* z, int64_t ldz );
int64_t LAPACKE_dgghrd_work_64( int matrix_layout, char compq, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                double* a, int64_t lda, double* b,
                                int64_t ldb, double* q, int64_t ldq,
                                double* z, int64_t ldz );
int64_t LAPACKE_cgghrd_work_64( int matrix_layout, char compq, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* q, int64_t ldq,
                                lapack_complex_float* z, int64_t ldz );
int64_t LAPACKE_zgghrd_work_64( int matrix_layout, char compq, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* q, int64_t ldq,
                                lapack_complex_double* z, int64_t ldz );

int64_t LAPACKE_sgghd3_work_64( int matrix_layout, char compq, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                float* a, int64_t lda,
                                float* b, int64_t ldb,
                                float* q, int64_t ldq,
                                float* z, int64_t ldz,
                                float* work, int64_t lwork );
int64_t LAPACKE_dgghd3_work_64( int matrix_layout, char compq, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                double* a, int64_t lda,
                                double* b, int64_t ldb,
                                double* q, int64_t ldq,
                                double* z, int64_t ldz,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgghd3_work_64( int matrix_layout, char compq, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* q, int64_t ldq,
                                lapack_complex_float* z, int64_t ldz,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgghd3_work_64( int matrix_layout, char compq, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* q, int64_t ldq,
                                lapack_complex_double* z, int64_t ldz,
                                lapack_complex_double* work,
                                int64_t lwork );

int64_t LAPACKE_sgglse_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t p, float* a, int64_t lda,
                                float* b, int64_t ldb, float* c, float* d,
                                float* x, float* work, int64_t lwork );
int64_t LAPACKE_dgglse_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t p, double* a, int64_t lda,
                                double* b, int64_t ldb, double* c, double* d,
                                double* x, double* work, int64_t lwork );
int64_t LAPACKE_cgglse_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t p, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* c,
                                lapack_complex_float* d,
                                lapack_complex_float* x,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgglse_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t p, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb, lapack_complex_double* c,
                                lapack_complex_double* d,
                                lapack_complex_double* x,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sggqrf_work_64( int matrix_layout, int64_t n, int64_t m,
                                int64_t p, float* a, int64_t lda,
                                float* taua, float* b, int64_t ldb,
                                float* taub, float* work, int64_t lwork );
int64_t LAPACKE_dggqrf_work_64( int matrix_layout, int64_t n, int64_t m,
                                int64_t p, double* a, int64_t lda,
                                double* taua, double* b, int64_t ldb,
                                double* taub, double* work, int64_t lwork );
int64_t LAPACKE_cggqrf_work_64( int matrix_layout, int64_t n, int64_t m,
                                int64_t p, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* taua,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* taub,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zggqrf_work_64( int matrix_layout, int64_t n, int64_t m,
                                int64_t p, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* taua,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* taub,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sggrqf_work_64( int matrix_layout, int64_t m, int64_t p,
                                int64_t n, float* a, int64_t lda,
                                float* taua, float* b, int64_t ldb,
                                float* taub, float* work, int64_t lwork );
int64_t LAPACKE_dggrqf_work_64( int matrix_layout, int64_t m, int64_t p,
                                int64_t n, double* a, int64_t lda,
                                double* taua, double* b, int64_t ldb,
                                double* taub, double* work, int64_t lwork );
int64_t LAPACKE_cggrqf_work_64( int matrix_layout, int64_t m, int64_t p,
                                int64_t n, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* taua,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* taub,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zggrqf_work_64( int matrix_layout, int64_t m, int64_t p,
                                int64_t n, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* taua,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* taub,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sggsvd_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t n,
                                int64_t p, int64_t* k, int64_t* l,
                                float* a, int64_t lda, float* b,
                                int64_t ldb, float* alpha, float* beta,
                                float* u, int64_t ldu, float* v,
                                int64_t ldv, float* q, int64_t ldq,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dggsvd_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t n,
                                int64_t p, int64_t* k, int64_t* l,
                                double* a, int64_t lda, double* b,
                                int64_t ldb, double* alpha, double* beta,
                                double* u, int64_t ldu, double* v,
                                int64_t ldv, double* q, int64_t ldq,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cggsvd_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t n,
                                int64_t p, int64_t* k, int64_t* l,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                float* alpha, float* beta,
                                lapack_complex_float* u, int64_t ldu,
                                lapack_complex_float* v, int64_t ldv,
                                lapack_complex_float* q, int64_t ldq,
                                lapack_complex_float* work, float* rwork,
                                int64_t* iwork );
int64_t LAPACKE_zggsvd_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t n,
                                int64_t p, int64_t* k, int64_t* l,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                double* alpha, double* beta,
                                lapack_complex_double* u, int64_t ldu,
                                lapack_complex_double* v, int64_t ldv,
                                lapack_complex_double* q, int64_t ldq,
                                lapack_complex_double* work, double* rwork,
                                int64_t* iwork );

int64_t LAPACKE_sggsvd3_work_64( int matrix_layout, char jobu, char jobv,
                                 char jobq, int64_t m, int64_t n,
                                 int64_t p, int64_t* k, int64_t* l,
                                 float* a, int64_t lda, float* b,
                                 int64_t ldb, float* alpha, float* beta,
                                 float* u, int64_t ldu, float* v,
                                 int64_t ldv, float* q, int64_t ldq,
                                 float* work, int64_t lwork,
                                 int64_t* iwork );
int64_t LAPACKE_dggsvd3_work_64( int matrix_layout, char jobu, char jobv,
                                 char jobq, int64_t m, int64_t n,
                                 int64_t p, int64_t* k, int64_t* l,
                                 double* a, int64_t lda, double* b,
                                 int64_t ldb, double* alpha, double* beta,
                                 double* u, int64_t ldu, double* v,
                                 int64_t ldv, double* q, int64_t ldq,
                                 double* work, int64_t lwork,
                                 int64_t* iwork );
int64_t LAPACKE_cggsvd3_work_64( int matrix_layout, char jobu, char jobv,
                                 char jobq, int64_t m, int64_t n,
                                 int64_t p, int64_t* k, int64_t* l,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* b, int64_t ldb,
                                 float* alpha, float* beta,
                                 lapack_complex_float* u, int64_t ldu,
                                 lapack_complex_float* v, int64_t ldv,
                                 lapack_complex_float* q, int64_t ldq,
                                 lapack_complex_float* work, int64_t lwork,
                                 float* rwork, int64_t* iwork );
int64_t LAPACKE_zggsvd3_work_64( int matrix_layout, char jobu, char jobv,
                                 char jobq, int64_t m, int64_t n,
                                 int64_t p, int64_t* k, int64_t* l,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* b, int64_t ldb,
                                 double* alpha, double* beta,
                                 lapack_complex_double* u, int64_t ldu,
                                 lapack_complex_double* v, int64_t ldv,
                                 lapack_complex_double* q, int64_t ldq,
                                 lapack_complex_double* work, int64_t lwork,
                                 double* rwork, int64_t* iwork );

int64_t LAPACKE_sggsvp_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t p,
                                int64_t n, float* a, int64_t lda,
                                float* b, int64_t ldb, float tola,
                                float tolb, int64_t* k, int64_t* l,
                                float* u, int64_t ldu, float* v,
                                int64_t ldv, float* q, int64_t ldq,
                                int64_t* iwork, float* tau, float* work );
int64_t LAPACKE_dggsvp_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t p,
                                int64_t n, double* a, int64_t lda,
                                double* b, int64_t ldb, double tola,
                                double tolb, int64_t* k, int64_t* l,
                                double* u, int64_t ldu, double* v,
                                int64_t ldv, double* q, int64_t ldq,
                                int64_t* iwork, double* tau, double* work );
int64_t LAPACKE_cggsvp_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t p,
                                int64_t n, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* b,
                                int64_t ldb, float tola, float tolb,
                                int64_t* k, int64_t* l,
                                lapack_complex_float* u, int64_t ldu,
                                lapack_complex_float* v, int64_t ldv,
                                lapack_complex_float* q, int64_t ldq,
                                int64_t* iwork, float* rwork,
                                lapack_complex_float* tau,
                                lapack_complex_float* work );
int64_t LAPACKE_zggsvp_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t p,
                                int64_t n, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb, double tola, double tolb,
                                int64_t* k, int64_t* l,
                                lapack_complex_double* u, int64_t ldu,
                                lapack_complex_double* v, int64_t ldv,
                                lapack_complex_double* q, int64_t ldq,
                                int64_t* iwork, double* rwork,
                                lapack_complex_double* tau,
                                lapack_complex_double* work );

int64_t LAPACKE_sggsvp3_work_64( int matrix_layout, char jobu, char jobv,
                                 char jobq, int64_t m, int64_t p,
                                 int64_t n, float* a, int64_t lda,
                                 float* b, int64_t ldb, float tola,
                                 float tolb, int64_t* k, int64_t* l,
                                 float* u, int64_t ldu, float* v,
                                 int64_t ldv, float* q, int64_t ldq,
                                 int64_t* iwork, float* tau,
                                 float* work, int64_t lwork );
int64_t LAPACKE_dggsvp3_work_64( int matrix_layout, char jobu, char jobv,
                                 char jobq, int64_t m, int64_t p,
                                 int64_t n, double* a, int64_t lda,
                                 double* b, int64_t ldb, double tola,
                                 double tolb, int64_t* k, int64_t* l,
                                 double* u, int64_t ldu, double* v,
                                 int64_t ldv, double* q, int64_t ldq,
                                 int64_t* iwork, double* tau, double* work,
                                 int64_t lwork );
int64_t LAPACKE_cggsvp3_work_64( int matrix_layout, char jobu, char jobv,
                                 char jobq, int64_t m, int64_t p,
                                 int64_t n, lapack_complex_float* a,
                                 int64_t lda, lapack_complex_float* b,
                                 int64_t ldb, float tola, float tolb,
                                 int64_t* k, int64_t* l,
                                 lapack_complex_float* u, int64_t ldu,
                                 lapack_complex_float* v, int64_t ldv,
                                 lapack_complex_float* q, int64_t ldq,
                                 int64_t* iwork, float* rwork,
                                 lapack_complex_float* tau,
                                 lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zggsvp3_work_64( int matrix_layout, char jobu, char jobv,
                                 char jobq, int64_t m, int64_t p,
                                 int64_t n, lapack_complex_double* a,
                                 int64_t lda, lapack_complex_double* b,
                                 int64_t ldb, double tola, double tolb,
                                 int64_t* k, int64_t* l,
                                 lapack_complex_double* u, int64_t ldu,
                                 lapack_complex_double* v, int64_t ldv,
                                 lapack_complex_double* q, int64_t ldq,
                                 int64_t* iwork, double* rwork,
                                 lapack_complex_double* tau,
                                 lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgtcon_work_64( char norm, int64_t n, const float* dl,
                                const float* d, const float* du,
                                const float* du2, const int64_t* ipiv,
                                float anorm, float* rcond, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dgtcon_work_64( char norm, int64_t n, const double* dl,
                                const double* d, const double* du,
                                const double* du2, const int64_t* ipiv,
                                double anorm, double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_cgtcon_work_64( char norm, int64_t n,
                                const lapack_complex_float* dl,
                                const lapack_complex_float* d,
                                const lapack_complex_float* du,
                                const lapack_complex_float* du2,
                                const int64_t* ipiv, float anorm,
                                float* rcond, lapack_complex_float* work );
int64_t LAPACKE_zgtcon_work_64( char norm, int64_t n,
                                const lapack_complex_double* dl,
                                const lapack_complex_double* d,
                                const lapack_complex_double* du,
                                const lapack_complex_double* du2,
                                const int64_t* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

int64_t LAPACKE_sgtrfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const float* dl,
                                const float* d, const float* du,
                                const float* dlf, const float* df,
                                const float* duf, const float* du2,
                                const int64_t* ipiv, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* ferr, float* berr, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dgtrfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const double* dl,
                                const double* d, const double* du,
                                const double* dlf, const double* df,
                                const double* duf, const double* du2,
                                const int64_t* ipiv, const double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* ferr, double* berr, double* work,
                                int64_t* iwork );
int64_t LAPACKE_cgtrfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const lapack_complex_float* dl,
                                const lapack_complex_float* d,
                                const lapack_complex_float* du,
                                const lapack_complex_float* dlf,
                                const lapack_complex_float* df,
                                const lapack_complex_float* duf,
                                const lapack_complex_float* du2,
                                const int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zgtrfs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs,
                                const lapack_complex_double* dl,
                                const lapack_complex_double* d,
                                const lapack_complex_double* du,
                                const lapack_complex_double* dlf,
                                const lapack_complex_double* df,
                                const lapack_complex_double* duf,
                                const lapack_complex_double* du2,
                                const int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sgtsv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               float* dl, float* d, float* du, float* b,
                               int64_t ldb );
int64_t LAPACKE_dgtsv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               double* dl, double* d, double* du, double* b,
                               int64_t ldb );
int64_t LAPACKE_cgtsv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               lapack_complex_float* dl,
                               lapack_complex_float* d,
                               lapack_complex_float* du,
                               lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zgtsv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               lapack_complex_double* dl,
                               lapack_complex_double* d,
                               lapack_complex_double* du,
                               lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sgtsvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t nrhs, const float* dl,
                                const float* d, const float* du, float* dlf,
                                float* df, float* duf, float* du2,
                                int64_t* ipiv, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dgtsvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t nrhs, const double* dl,
                                const double* d, const double* du, double* dlf,
                                double* df, double* duf, double* du2,
                                int64_t* ipiv, const double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cgtsvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_float* dl,
                                const lapack_complex_float* d,
                                const lapack_complex_float* du,
                                lapack_complex_float* dlf,
                                lapack_complex_float* df,
                                lapack_complex_float* duf,
                                lapack_complex_float* du2, int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zgtsvx_work_64( int matrix_layout, char fact, char trans,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_double* dl,
                                const lapack_complex_double* d,
                                const lapack_complex_double* du,
                                lapack_complex_double* dlf,
                                lapack_complex_double* df,
                                lapack_complex_double* duf,
                                lapack_complex_double* du2, int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sgttrf_work_64( int64_t n, float* dl, float* d, float* du,
                                float* du2, int64_t* ipiv );
int64_t LAPACKE_dgttrf_work_64( int64_t n, double* dl, double* d, double* du,
                                double* du2, int64_t* ipiv );
int64_t LAPACKE_cgttrf_work_64( int64_t n, lapack_complex_float* dl,
                                lapack_complex_float* d,
                                lapack_complex_float* du,
                                lapack_complex_float* du2, int64_t* ipiv );
int64_t LAPACKE_zgttrf_work_64( int64_t n, lapack_complex_double* dl,
                                lapack_complex_double* d,
                                lapack_complex_double* du,
                                lapack_complex_double* du2, int64_t* ipiv );

int64_t LAPACKE_sgttrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const float* dl,
                                const float* d, const float* du,
                                const float* du2, const int64_t* ipiv,
                                float* b, int64_t ldb );
int64_t LAPACKE_dgttrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const double* dl,
                                const double* d, const double* du,
                                const double* du2, const int64_t* ipiv,
                                double* b, int64_t ldb );
int64_t LAPACKE_cgttrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs, const lapack_complex_float* dl,
                                const lapack_complex_float* d,
                                const lapack_complex_float* du,
                                const lapack_complex_float* du2,
                                const int64_t* ipiv, lapack_complex_float* b,
                                int64_t ldb );
int64_t LAPACKE_zgttrs_work_64( int matrix_layout, char trans, int64_t n,
                                int64_t nrhs,
                                const lapack_complex_double* dl,
                                const lapack_complex_double* d,
                                const lapack_complex_double* du,
                                const lapack_complex_double* du2,
                                const int64_t* ipiv,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_chbev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t kd,
                               lapack_complex_float* ab, int64_t ldab,
                               float* w, lapack_complex_float* z,
                               int64_t ldz, lapack_complex_float* work,
                               float* rwork );
int64_t LAPACKE_zhbev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t kd,
                               lapack_complex_double* ab, int64_t ldab,
                               double* w, lapack_complex_double* z,
                               int64_t ldz, lapack_complex_double* work,
                               double* rwork );

int64_t LAPACKE_chbevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t kd,
                                lapack_complex_float* ab, int64_t ldab,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_zhbevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t kd,
                                lapack_complex_double* ab, int64_t ldab,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );

int64_t LAPACKE_chbevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t kd,
                                lapack_complex_float* ab, int64_t ldab,
                                lapack_complex_float* q, int64_t ldq,
                                float vl, float vu, int64_t il,
                                int64_t iu, float abstol, int64_t* m,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                float* rwork, int64_t* iwork,
                                int64_t* ifail );
int64_t LAPACKE_zhbevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t kd,
                                lapack_complex_double* ab, int64_t ldab,
                                lapack_complex_double* q, int64_t ldq,
                                double vl, double vu, int64_t il,
                                int64_t iu, double abstol, int64_t* m,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                double* rwork, int64_t* iwork,
                                int64_t* ifail );

int64_t LAPACKE_chbgst_work_64( int matrix_layout, char vect, char uplo,
                                int64_t n, int64_t ka, int64_t kb,
                                lapack_complex_float* ab, int64_t ldab,
                                const lapack_complex_float* bb, int64_t ldbb,
                                lapack_complex_float* x, int64_t ldx,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zhbgst_work_64( int matrix_layout, char vect, char uplo,
                                int64_t n, int64_t ka, int64_t kb,
                                lapack_complex_double* ab, int64_t ldab,
                                const lapack_complex_double* bb,
                                int64_t ldbb, lapack_complex_double* x,
                                int64_t ldx, lapack_complex_double* work,
                                double* rwork );

int64_t LAPACKE_chbgv_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t ka, int64_t kb,
                               lapack_complex_float* ab, int64_t ldab,
                               lapack_complex_float* bb, int64_t ldbb,
                               float* w, lapack_complex_float* z,
                               int64_t ldz, lapack_complex_float* work,
                               float* rwork );
int64_t LAPACKE_zhbgv_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t ka, int64_t kb,
                               lapack_complex_double* ab, int64_t ldab,
                               lapack_complex_double* bb, int64_t ldbb,
                               double* w, lapack_complex_double* z,
                               int64_t ldz, lapack_complex_double* work,
                               double* rwork );

int64_t LAPACKE_chbgvd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t ka, int64_t kb,
                                lapack_complex_float* ab, int64_t ldab,
                                lapack_complex_float* bb, int64_t ldbb,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_zhbgvd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t ka, int64_t kb,
                                lapack_complex_double* ab, int64_t ldab,
                                lapack_complex_double* bb, int64_t ldbb,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );

int64_t LAPACKE_chbgvx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t ka,
                                int64_t kb, lapack_complex_float* ab,
                                int64_t ldab, lapack_complex_float* bb,
                                int64_t ldbb, lapack_complex_float* q,
                                int64_t ldq, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w,
                                lapack_complex_float* z, int64_t ldz,
                                lapack_complex_float* work, float* rwork,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_zhbgvx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t ka,
                                int64_t kb, lapack_complex_double* ab,
                                int64_t ldab, lapack_complex_double* bb,
                                int64_t ldbb, lapack_complex_double* q,
                                int64_t ldq, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w,
                                lapack_complex_double* z, int64_t ldz,
                                lapack_complex_double* work, double* rwork,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_chbtrd_work_64( int matrix_layout, char vect, char uplo,
                                int64_t n, int64_t kd,
                                lapack_complex_float* ab, int64_t ldab,
                                float* d, float* e, lapack_complex_float* q,
                                int64_t ldq, lapack_complex_float* work );
int64_t LAPACKE_zhbtrd_work_64( int matrix_layout, char vect, char uplo,
                                int64_t n, int64_t kd,
                                lapack_complex_double* ab, int64_t ldab,
                                double* d, double* e, lapack_complex_double* q,
                                int64_t ldq, lapack_complex_double* work );

int64_t LAPACKE_checon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                const int64_t* ipiv, float anorm,
                                float* rcond, lapack_complex_float* work );
int64_t LAPACKE_zhecon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                const int64_t* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

int64_t LAPACKE_cheequb_work_64( int matrix_layout, char uplo, int64_t n,
                                 const lapack_complex_float* a, int64_t lda,
                                 float* s, float* scond, float* amax,
                                 lapack_complex_float* work );
int64_t LAPACKE_zheequb_work_64( int matrix_layout, char uplo, int64_t n,
                                 const lapack_complex_double* a, int64_t lda,
                                 double* s, double* scond, double* amax,
                                 lapack_complex_double* work );

int64_t LAPACKE_cheev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, lapack_complex_float* a,
                               int64_t lda, float* w,
                               lapack_complex_float* work, int64_t lwork,
                               float* rwork );
int64_t LAPACKE_zheev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, lapack_complex_double* a,
                               int64_t lda, double* w,
                               lapack_complex_double* work, int64_t lwork,
                               double* rwork );

int64_t LAPACKE_cheevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, lapack_complex_float* a,
                                int64_t lda, float* w,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_zheevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, lapack_complex_double* a,
                                int64_t lda, double* w,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_cheevr_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                float vl, float vu, int64_t il,
                                int64_t iu, float abstol, int64_t* m,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, int64_t* isuppz,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_zheevr_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                double vl, double vu, int64_t il,
                                int64_t iu, double abstol, int64_t* m,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, int64_t* isuppz,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_cheevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                float vl, float vu, int64_t il,
                                int64_t iu, float abstol, int64_t* m,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_zheevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                double vl, double vu, int64_t il,
                                int64_t iu, double abstol, int64_t* m,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_chegst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* b,
                                int64_t ldb );
int64_t LAPACKE_zhegst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, lapack_complex_double* a,
                                int64_t lda, const lapack_complex_double* b,
                                int64_t ldb );

int64_t LAPACKE_chegv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* b,
                               int64_t ldb, float* w,
                               lapack_complex_float* work, int64_t lwork,
                               float* rwork );
int64_t LAPACKE_zhegv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n,
                               lapack_complex_double* a, int64_t lda,
                               lapack_complex_double* b, int64_t ldb,
                               double* w, lapack_complex_double* work,
                               int64_t lwork, double* rwork );

int64_t LAPACKE_chegvd_work_64( int matrix_layout, int64_t itype, char jobz,
                                char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                float* w, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_zhegvd_work_64( int matrix_layout, int64_t itype, char jobz,
                                char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                double* w, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );

int64_t LAPACKE_chegvx_work_64( int matrix_layout, int64_t itype, char jobz,
                                char range, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                float vl, float vu, int64_t il,
                                int64_t iu, float abstol, int64_t* m,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_zhegvx_work_64( int matrix_layout, int64_t itype, char jobz,
                                char range, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                double vl, double vu, int64_t il,
                                int64_t iu, double abstol, int64_t* m,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_cherfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* af,
                                int64_t ldaf, const int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zherfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, const lapack_complex_double* af,
                                int64_t ldaf, const int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_cherfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs,
                                 const lapack_complex_float* a, int64_t lda,
                                 const lapack_complex_float* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const float* s, const lapack_complex_float* b,
                                 int64_t ldb, lapack_complex_float* x,
                                 int64_t ldx, float* rcond, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, lapack_complex_float* work,
                                 float* rwork );
int64_t LAPACKE_zherfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs,
                                 const lapack_complex_double* a, int64_t lda,
                                 const lapack_complex_double* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const double* s,
                                 const lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* x, int64_t ldx,
                                 double* rcond, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, lapack_complex_double* work,
                                 double* rwork );

int64_t LAPACKE_chesv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_float* b, int64_t ldb,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zhesv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_double* b, int64_t ldb,
                               lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_chesvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* af, int64_t ldaf,
                                int64_t* ipiv, const lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* x,
                                int64_t ldx, float* rcond, float* ferr,
                                float* berr, lapack_complex_float* work,
                                int64_t lwork, float* rwork );
int64_t LAPACKE_zhesvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* af, int64_t ldaf,
                                int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork );

int64_t LAPACKE_chesvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, float* s,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* x, int64_t ldx,
                                 float* rcond, float* rpvgrw, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, lapack_complex_float* work,
                                 float* rwork );
int64_t LAPACKE_zhesvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, double* s,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* x, int64_t ldx,
                                 double* rcond, double* rpvgrw, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, lapack_complex_double* work,
                                 double* rwork );

int64_t LAPACKE_chetrd_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                float* d, float* e, lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zhetrd_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                double* d, double* e,
                                lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_chetrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_float* work,
                                int64_t lwork );
int64_t LAPACKE_zhetrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_double* work,
                                int64_t lwork );

int64_t LAPACKE_chetri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                const int64_t* ipiv,
                                lapack_complex_float* work );
int64_t LAPACKE_zhetri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                const int64_t* ipiv,
                                lapack_complex_double* work );

int64_t LAPACKE_chetrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zhetrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_chfrk_work_64( int matrix_layout, char transr, char uplo,
                               char trans, int64_t n, int64_t k,
                               float alpha, const lapack_complex_float* a,
                               int64_t lda, float beta,
                               lapack_complex_float* c );
int64_t LAPACKE_zhfrk_work_64( int matrix_layout, char transr, char uplo,
                               char trans, int64_t n, int64_t k,
                               double alpha, const lapack_complex_double* a,
                               int64_t lda, double beta,
                               lapack_complex_double* c );

int64_t LAPACKE_shgeqz_work_64( int matrix_layout, char job, char compq,
                                char compz, int64_t n, int64_t ilo,
                                int64_t ihi, float* h, int64_t ldh,
                                float* t, int64_t ldt, float* alphar,
                                float* alphai, float* beta, float* q,
                                int64_t ldq, float* z, int64_t ldz,
                                float* work, int64_t lwork );
int64_t LAPACKE_dhgeqz_work_64( int matrix_layout, char job, char compq,
                                char compz, int64_t n, int64_t ilo,
                                int64_t ihi, double* h, int64_t ldh,
                                double* t, int64_t ldt, double* alphar,
                                double* alphai, double* beta, double* q,
                                int64_t ldq, double* z, int64_t ldz,
                                double* work, int64_t lwork );
int64_t LAPACKE_chgeqz_work_64( int matrix_layout, char job, char compq,
                                char compz, int64_t n, int64_t ilo,
                                int64_t ihi, lapack_complex_float* h,
                                int64_t ldh, lapack_complex_float* t,
                                int64_t ldt, lapack_complex_float* alpha,
                                lapack_complex_float* beta,
                                lapack_complex_float* q, int64_t ldq,
                                lapack_complex_float* z, int64_t ldz,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork );
int64_t LAPACKE_zhgeqz_work_64( int matrix_layout, char job, char compq,
                                char compz, int64_t n, int64_t ilo,
                                int64_t ihi, lapack_complex_double* h,
                                int64_t ldh, lapack_complex_double* t,
                                int64_t ldt, lapack_complex_double* alpha,
                                lapack_complex_double* beta,
                                lapack_complex_double* q, int64_t ldq,
                                lapack_complex_double* z, int64_t ldz,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork );

int64_t LAPACKE_chpcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* ap,
                                const int64_t* ipiv, float anorm,
                                float* rcond, lapack_complex_float* work );
int64_t LAPACKE_zhpcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* ap,
                                const int64_t* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

int64_t LAPACKE_chpev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, lapack_complex_float* ap, float* w,
                               lapack_complex_float* z, int64_t ldz,
                               lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zhpev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, lapack_complex_double* ap,
                               double* w, lapack_complex_double* z,
                               int64_t ldz, lapack_complex_double* work,
                               double* rwork );

int64_t LAPACKE_chpevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, lapack_complex_float* ap,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_zhpevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, lapack_complex_double* ap,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );

int64_t LAPACKE_chpevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_float* ap, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w,
                                lapack_complex_float* z, int64_t ldz,
                                lapack_complex_float* work, float* rwork,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_zhpevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_double* ap, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w,
                                lapack_complex_double* z, int64_t ldz,
                                lapack_complex_double* work, double* rwork,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_chpgst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, lapack_complex_float* ap,
                                const lapack_complex_float* bp );
int64_t LAPACKE_zhpgst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, lapack_complex_double* ap,
                                const lapack_complex_double* bp );

int64_t LAPACKE_chpgv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n,
                               lapack_complex_float* ap,
                               lapack_complex_float* bp, float* w,
                               lapack_complex_float* z, int64_t ldz,
                               lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zhpgv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n,
                               lapack_complex_double* ap,
                               lapack_complex_double* bp, double* w,
                               lapack_complex_double* z, int64_t ldz,
                               lapack_complex_double* work, double* rwork );

int64_t LAPACKE_chpgvd_work_64( int matrix_layout, int64_t itype, char jobz,
                                char uplo, int64_t n,
                                lapack_complex_float* ap,
                                lapack_complex_float* bp, float* w,
                                lapack_complex_float* z, int64_t ldz,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_zhpgvd_work_64( int matrix_layout, int64_t itype, char jobz,
                                char uplo, int64_t n,
                                lapack_complex_double* ap,
                                lapack_complex_double* bp, double* w,
                                lapack_complex_double* z, int64_t ldz,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_chpgvx_work_64( int matrix_layout, int64_t itype, char jobz,
                                char range, char uplo, int64_t n,
                                lapack_complex_float* ap,
                                lapack_complex_float* bp, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w,
                                lapack_complex_float* z, int64_t ldz,
                                lapack_complex_float* work, float* rwork,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_zhpgvx_work_64( int matrix_layout, int64_t itype, char jobz,
                                char range, char uplo, int64_t n,
                                lapack_complex_double* ap,
                                lapack_complex_double* bp, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w,
                                lapack_complex_double* z, int64_t ldz,
                                lapack_complex_double* work, double* rwork,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_chprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* ap,
                                const lapack_complex_float* afp,
                                const int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zhprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs,
                                const lapack_complex_double* ap,
                                const lapack_complex_double* afp,
                                const int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_chpsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* ap,
                               int64_t* ipiv, lapack_complex_float* b,
                               int64_t ldb );
int64_t LAPACKE_zhpsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* ap,
                               int64_t* ipiv, lapack_complex_double* b,
                               int64_t ldb );

int64_t LAPACKE_chpsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_float* ap,
                                lapack_complex_float* afp, int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zhpsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_double* ap,
                                lapack_complex_double* afp, int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_chptrd_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* ap, float* d, float* e,
                                lapack_complex_float* tau );
int64_t LAPACKE_zhptrd_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* ap, double* d, double* e,
                                lapack_complex_double* tau );

int64_t LAPACKE_chptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* ap, int64_t* ipiv );
int64_t LAPACKE_zhptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* ap, int64_t* ipiv );

int64_t LAPACKE_chptri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* ap,
                                const int64_t* ipiv,
                                lapack_complex_float* work );
int64_t LAPACKE_zhptri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* ap,
                                const int64_t* ipiv,
                                lapack_complex_double* work );

int64_t LAPACKE_chptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* ap,
                                const int64_t* ipiv, lapack_complex_float* b,
                                int64_t ldb );
int64_t LAPACKE_zhptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs,
                                const lapack_complex_double* ap,
                                const int64_t* ipiv,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_shsein_work_64( int matrix_layout, char job, char eigsrc,
                                char initv, lapack_logical* select,
                                int64_t n, const float* h, int64_t ldh,
                                float* wr, const float* wi, float* vl,
                                int64_t ldvl, float* vr, int64_t ldvr,
                                int64_t mm, int64_t* m, float* work,
                                int64_t* ifaill, int64_t* ifailr );
int64_t LAPACKE_dhsein_work_64( int matrix_layout, char job, char eigsrc,
                                char initv, lapack_logical* select,
                                int64_t n, const double* h, int64_t ldh,
                                double* wr, const double* wi, double* vl,
                                int64_t ldvl, double* vr, int64_t ldvr,
                                int64_t mm, int64_t* m, double* work,
                                int64_t* ifaill, int64_t* ifailr );
int64_t LAPACKE_chsein_work_64( int matrix_layout, char job, char eigsrc,
                                char initv, const lapack_logical* select,
                                int64_t n, const lapack_complex_float* h,
                                int64_t ldh, lapack_complex_float* w,
                                lapack_complex_float* vl, int64_t ldvl,
                                lapack_complex_float* vr, int64_t ldvr,
                                int64_t mm, int64_t* m,
                                lapack_complex_float* work, float* rwork,
                                int64_t* ifaill, int64_t* ifailr );
int64_t LAPACKE_zhsein_work_64( int matrix_layout, char job, char eigsrc,
                                char initv, const lapack_logical* select,
                                int64_t n, const lapack_complex_double* h,
                                int64_t ldh, lapack_complex_double* w,
                                lapack_complex_double* vl, int64_t ldvl,
                                lapack_complex_double* vr, int64_t ldvr,
                                int64_t mm, int64_t* m,
                                lapack_complex_double* work, double* rwork,
                                int64_t* ifaill, int64_t* ifailr );

int64_t LAPACKE_shseqr_work_64( int matrix_layout, char job, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                float* h, int64_t ldh, float* wr, float* wi,
                                float* z, int64_t ldz, float* work,
                                int64_t lwork );
int64_t LAPACKE_dhseqr_work_64( int matrix_layout, char job, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                double* h, int64_t ldh, double* wr,
                                double* wi, double* z, int64_t ldz,
                                double* work, int64_t lwork );
int64_t LAPACKE_chseqr_work_64( int matrix_layout, char job, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                lapack_complex_float* h, int64_t ldh,
                                lapack_complex_float* w,
                                lapack_complex_float* z, int64_t ldz,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zhseqr_work_64( int matrix_layout, char job, char compz,
                                int64_t n, int64_t ilo, int64_t ihi,
                                lapack_complex_double* h, int64_t ldh,
                                lapack_complex_double* w,
                                lapack_complex_double* z, int64_t ldz,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_clacgv_work_64( int64_t n, lapack_complex_float* x,
                                int64_t incx );
int64_t LAPACKE_zlacgv_work_64( int64_t n, lapack_complex_double* x,
                                int64_t incx );

int64_t LAPACKE_slacn2_work_64( int64_t n, float* v, float* x,
                                int64_t* isgn, float* est, int64_t* kase,
                                int64_t* isave );
int64_t LAPACKE_dlacn2_work_64( int64_t n, double* v, double* x,
                                int64_t* isgn, double* est, int64_t* kase,
                                int64_t* isave );
int64_t LAPACKE_clacn2_work_64( int64_t n, lapack_complex_float* v,
                                lapack_complex_float* x,
                                float* est, int64_t* kase,
                                int64_t* isave );
int64_t LAPACKE_zlacn2_work_64( int64_t n, lapack_complex_double* v,
                                lapack_complex_double* x,
                                double* est, int64_t* kase,
                                int64_t* isave );

int64_t LAPACKE_slacpy_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, const float* a, int64_t lda,
                                float* b, int64_t ldb );
int64_t LAPACKE_dlacpy_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, const double* a, int64_t lda,
                                double* b, int64_t ldb );
int64_t LAPACKE_clacpy_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, const lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* b,
                                int64_t ldb );
int64_t LAPACKE_zlacpy_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, const lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb );

int64_t LAPACKE_clacp2_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, const float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zlacp2_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, const double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_zlag2c_work_64( int matrix_layout, int64_t m, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                lapack_complex_float* sa, int64_t ldsa );

int64_t LAPACKE_slag2d_work_64( int matrix_layout, int64_t m, int64_t n,
                                const float* sa, int64_t ldsa, double* a,
                                int64_t lda );

int64_t LAPACKE_dlag2s_work_64( int matrix_layout, int64_t m, int64_t n,
                                const double* a, int64_t lda, float* sa,
                                int64_t ldsa );

int64_t LAPACKE_clag2z_work_64( int matrix_layout, int64_t m, int64_t n,
                                const lapack_complex_float* sa, int64_t ldsa,
                                lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_slagge_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku, const float* d,
                                float* a, int64_t lda, int64_t* iseed,
                                float* work );
int64_t LAPACKE_dlagge_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku, const double* d,
                                double* a, int64_t lda, int64_t* iseed,
                                double* work );
int64_t LAPACKE_clagge_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku, const float* d,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* iseed, lapack_complex_float* work );
int64_t LAPACKE_zlagge_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t kl, int64_t ku, const double* d,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* iseed,
                                lapack_complex_double* work );

int64_t LAPACKE_claghe_work_64( int matrix_layout, int64_t n, int64_t k,
                                const float* d, lapack_complex_float* a,
                                int64_t lda, int64_t* iseed,
                                lapack_complex_float* work );
int64_t LAPACKE_zlaghe_work_64( int matrix_layout, int64_t n, int64_t k,
                                const double* d, lapack_complex_double* a,
                                int64_t lda, int64_t* iseed,
                                lapack_complex_double* work );

int64_t LAPACKE_slagsy_work_64( int matrix_layout, int64_t n, int64_t k,
                                const float* d, float* a, int64_t lda,
                                int64_t* iseed, float* work );
int64_t LAPACKE_dlagsy_work_64( int matrix_layout, int64_t n, int64_t k,
                                const double* d, double* a, int64_t lda,
                                int64_t* iseed, double* work );
int64_t LAPACKE_clagsy_work_64( int matrix_layout, int64_t n, int64_t k,
                                const float* d, lapack_complex_float* a,
                                int64_t lda, int64_t* iseed,
                                lapack_complex_float* work );
int64_t LAPACKE_zlagsy_work_64( int matrix_layout, int64_t n, int64_t k,
                                const double* d, lapack_complex_double* a,
                                int64_t lda, int64_t* iseed,
                                lapack_complex_double* work );

int64_t LAPACKE_slapmr_work_64( int matrix_layout, lapack_logical forwrd,
                                int64_t m, int64_t n, float* x,
                                int64_t ldx, int64_t* k );
int64_t LAPACKE_dlapmr_work_64( int matrix_layout, lapack_logical forwrd,
                                int64_t m, int64_t n, double* x,
                                int64_t ldx, int64_t* k );
int64_t LAPACKE_clapmr_work_64( int matrix_layout, lapack_logical forwrd,
                                int64_t m, int64_t n,
                                lapack_complex_float* x, int64_t ldx,
                                int64_t* k );
int64_t LAPACKE_zlapmr_work_64( int matrix_layout, lapack_logical forwrd,
                                int64_t m, int64_t n,
                                lapack_complex_double* x, int64_t ldx,
                                int64_t* k );

int64_t LAPACKE_slapmt_work_64( int matrix_layout, lapack_logical forwrd,
                                int64_t m, int64_t n, float* x,
                                int64_t ldx, int64_t* k );
int64_t LAPACKE_dlapmt_work_64( int matrix_layout, lapack_logical forwrd,
                                int64_t m, int64_t n, double* x,
                                int64_t ldx, int64_t* k );
int64_t LAPACKE_clapmt_work_64( int matrix_layout, lapack_logical forwrd,
                                int64_t m, int64_t n,
                                lapack_complex_float* x, int64_t ldx,
                                int64_t* k );
int64_t LAPACKE_zlapmt_work_64( int matrix_layout, lapack_logical forwrd,
                                int64_t m, int64_t n,
                                lapack_complex_double* x, int64_t ldx,
                                int64_t* k );

int64_t LAPACKE_slartgp_work_64( float f, float g, float* cs, float* sn,
                                 float* r );
int64_t LAPACKE_dlartgp_work_64( double f, double g, double* cs, double* sn,
                                 double* r );

int64_t LAPACKE_slartgs_work_64( float x, float y, float sigma, float* cs,
                                 float* sn );
int64_t LAPACKE_dlartgs_work_64( double x, double y, double sigma, double* cs,
                                 double* sn );

float LAPACKE_slapy2_work_64( float x, float y );
double LAPACKE_dlapy2_work_64( double x, double y );

float LAPACKE_slapy3_work_64( float x, float y, float z );
double LAPACKE_dlapy3_work_64( double x, double y, double z );

float LAPACKE_slamch_work_64( char cmach );
double LAPACKE_dlamch_work_64( char cmach );

float LAPACKE_slangb_work_64( int matrix_layout, char norm, int64_t n,
                           int64_t kl, int64_t ku, const float* ab,
                           int64_t ldab, float* work );
double LAPACKE_dlangb_work_64( int matrix_layout, char norm, int64_t n,
                            int64_t kl, int64_t ku, const double* ab,
                            int64_t ldab, double* work );
float LAPACKE_clangb_work_64( int matrix_layout, char norm, int64_t n,
                           int64_t kl, int64_t ku,
                           const lapack_complex_float* ab, int64_t ldab,
                           float* work );
double LAPACKE_zlangb_work_64( int matrix_layout, char norm, int64_t n,
                            int64_t kl, int64_t ku,
                            const lapack_complex_double* ab, int64_t ldab,
                            double* work );

float LAPACKE_slange_work_64( int matrix_layout, char norm, int64_t m,
                                int64_t n, const float* a, int64_t lda,
                                float* work );
double LAPACKE_dlange_work_64( int matrix_layout, char norm, int64_t m,
                                int64_t n, const double* a, int64_t lda,
                                double* work );
float LAPACKE_clange_work_64( int matrix_layout, char norm, int64_t m,
                                int64_t n, const lapack_complex_float* a,
                                int64_t lda, float* work );
double LAPACKE_zlange_work_64( int matrix_layout, char norm, int64_t m,
                                int64_t n, const lapack_complex_double* a,
                                int64_t lda, double* work );

float LAPACKE_clanhe_work_64( int matrix_layout, char norm, char uplo,
                                int64_t n, const lapack_complex_float* a,
                                int64_t lda, float* work );
double LAPACKE_zlanhe_work_64( int matrix_layout, char norm, char uplo,
                                int64_t n, const lapack_complex_double* a,
                                int64_t lda, double* work );

int64_t LAPACKE_clacrm_work_64( int matrix_layout, int64_t m, int64_t n,
                                const lapack_complex_float* a,
                                int64_t lda, const float* b,
                                int64_t ldb, lapack_complex_float* c,
                                int64_t ldc, float* work );
int64_t LAPACKE_zlacrm_work_64( int matrix_layout, int64_t m, int64_t n,
                                const lapack_complex_double* a,
                                int64_t lda, const double* b,
                                int64_t ldb, lapack_complex_double* c,
                                int64_t ldc, double* work );

int64_t LAPACKE_clarcm_work_64( int matrix_layout, int64_t m, int64_t n,
                                const float* a, int64_t lda,
                                const lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* c,
                                int64_t ldc, float* work );
int64_t LAPACKE_zlarcm_work_64( int matrix_layout, int64_t m, int64_t n,
                                const double* a, int64_t lda,
                                const lapack_complex_double* b,
                                int64_t ldb, lapack_complex_double* c,
                                int64_t ldc, double* work );

float LAPACKE_slansy_work_64( int matrix_layout, char norm, char uplo,
                                int64_t n, const float* a, int64_t lda,
                                float* work );
double LAPACKE_dlansy_work_64( int matrix_layout, char norm, char uplo,
                                int64_t n, const double* a, int64_t lda,
                                double* work );
float LAPACKE_clansy_work_64( int matrix_layout, char norm, char uplo,
                                int64_t n, const lapack_complex_float* a,
                                int64_t lda, float* work );
double LAPACKE_zlansy_work_64( int matrix_layout, char norm, char uplo,
                                int64_t n, const lapack_complex_double* a,
                                int64_t lda, double* work );

float LAPACKE_slanky_work_64( int matrix_layout, char norm, char uplo,
                                int64_t n, const float* a, int64_t lda,
                                float* work );
double LAPACKE_dlanky_work_64( int matrix_layout, char norm, char uplo,
                                int64_t n, const double* a, int64_t lda,
                                double* work );

float LAPACKE_slantr_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t m, int64_t n, const float* a,
                                int64_t lda, float* work );
double LAPACKE_dlantr_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t m, int64_t n,
                                const double* a, int64_t lda, double* work );
float LAPACKE_clantr_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t m, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                float* work );
double LAPACKE_zlantr_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t m, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                double* work );

int64_t LAPACKE_slarfb_work_64( int matrix_layout, char side, char trans,
                                char direct, char storev, int64_t m,
                                int64_t n, int64_t k, const float* v,
                                int64_t ldv, const float* t, int64_t ldt,
                                float* c, int64_t ldc, float* work,
                                int64_t ldwork );
int64_t LAPACKE_dlarfb_work_64( int matrix_layout, char side, char trans,
                                char direct, char storev, int64_t m,
                                int64_t n, int64_t k, const double* v,
                                int64_t ldv, const double* t, int64_t ldt,
                                double* c, int64_t ldc, double* work,
                                int64_t ldwork );
int64_t LAPACKE_clarfb_work_64( int matrix_layout, char side, char trans,
                                char direct, char storev, int64_t m,
                                int64_t n, int64_t k,
                                const lapack_complex_float* v, int64_t ldv,
                                const lapack_complex_float* t, int64_t ldt,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t ldwork );
int64_t LAPACKE_zlarfb_work_64( int matrix_layout, char side, char trans,
                                char direct, char storev, int64_t m,
                                int64_t n, int64_t k,
                                const lapack_complex_double* v, int64_t ldv,
                                const lapack_complex_double* t, int64_t ldt,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work,
                                int64_t ldwork );

int64_t LAPACKE_slarfg_work_64( int64_t n, float* alpha, float* x,
                                int64_t incx, float* tau );
int64_t LAPACKE_dlarfg_work_64( int64_t n, double* alpha, double* x,
                                int64_t incx, double* tau );
int64_t LAPACKE_clarfg_work_64( int64_t n, lapack_complex_float* alpha,
                                lapack_complex_float* x, int64_t incx,
                                lapack_complex_float* tau );
int64_t LAPACKE_zlarfg_work_64( int64_t n, lapack_complex_double* alpha,
                                lapack_complex_double* x, int64_t incx,
                                lapack_complex_double* tau );

int64_t LAPACKE_slarft_work_64( int matrix_layout, char direct, char storev,
                                int64_t n, int64_t k, const float* v,
                                int64_t ldv, const float* tau, float* t,
                                int64_t ldt );
int64_t LAPACKE_dlarft_work_64( int matrix_layout, char direct, char storev,
                                int64_t n, int64_t k, const double* v,
                                int64_t ldv, const double* tau, double* t,
                                int64_t ldt );
int64_t LAPACKE_clarft_work_64( int matrix_layout, char direct, char storev,
                                int64_t n, int64_t k,
                                const lapack_complex_float* v, int64_t ldv,
                                const lapack_complex_float* tau,
                                lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_zlarft_work_64( int matrix_layout, char direct, char storev,
                                int64_t n, int64_t k,
                                const lapack_complex_double* v, int64_t ldv,
                                const lapack_complex_double* tau,
                                lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_slarfx_work_64( int matrix_layout, char side, int64_t m,
                                int64_t n, const float* v, float tau,
                                float* c, int64_t ldc, float* work );
int64_t LAPACKE_dlarfx_work_64( int matrix_layout, char side, int64_t m,
                                int64_t n, const double* v, double tau,
                                double* c, int64_t ldc, double* work );
int64_t LAPACKE_clarfx_work_64( int matrix_layout, char side, int64_t m,
                                int64_t n, const lapack_complex_float* v,
                                lapack_complex_float tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work );
int64_t LAPACKE_zlarfx_work_64( int matrix_layout, char side, int64_t m,
                                int64_t n, const lapack_complex_double* v,
                                lapack_complex_double tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work );

int64_t LAPACKE_slarnv_work_64( int64_t idist, int64_t* iseed,
                                int64_t n, float* x );
int64_t LAPACKE_dlarnv_work_64( int64_t idist, int64_t* iseed,
                                int64_t n, double* x );
int64_t LAPACKE_clarnv_work_64( int64_t idist, int64_t* iseed,
                                int64_t n, lapack_complex_float* x );
int64_t LAPACKE_zlarnv_work_64( int64_t idist, int64_t* iseed,
                                int64_t n, lapack_complex_double* x );


int64_t LAPACKE_slascl_work_64( int matrix_layout, char type, int64_t kl,
                           int64_t ku, float cfrom, float cto,
                           int64_t m, int64_t n, float* a,
                           int64_t lda );
int64_t LAPACKE_dlascl_work_64( int matrix_layout, char type, int64_t kl,
                           int64_t ku, double cfrom, double cto,
                           int64_t m, int64_t n, double* a,
                           int64_t lda );
int64_t LAPACKE_clascl_work_64( int matrix_layout, char type, int64_t kl,
                           int64_t ku, float cfrom, float cto,
                           int64_t m, int64_t n, lapack_complex_float* a,
                           int64_t lda );
int64_t LAPACKE_zlascl_work_64( int matrix_layout, char type, int64_t kl,
                           int64_t ku, double cfrom, double cto,
                           int64_t m, int64_t n, lapack_complex_double* a,
                           int64_t lda );

int64_t LAPACKE_slaset_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, float alpha, float beta, float* a,
                                int64_t lda );
int64_t LAPACKE_dlaset_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, double alpha, double beta,
                                double* a, int64_t lda );
int64_t LAPACKE_claset_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, lapack_complex_float alpha,
                                lapack_complex_float beta,
                                lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zlaset_work_64( int matrix_layout, char uplo, int64_t m,
                                int64_t n, lapack_complex_double alpha,
                                lapack_complex_double beta,
                                lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_slasrt_work_64( char id, int64_t n, float* d );
int64_t LAPACKE_dlasrt_work_64( char id, int64_t n, double* d );

int64_t LAPACKE_slassq_work_64( int64_t n,                 float* x, int64_t incx,  float* scale,  float* sumsq );
int64_t LAPACKE_dlassq_work_64( int64_t n,                double* x, int64_t incx, double* scale, double* sumsq );
int64_t LAPACKE_classq_work_64( int64_t n,  lapack_complex_float* x, int64_t incx,  float* scale,  float* sumsq );
int64_t LAPACKE_zlassq_work_64( int64_t n, lapack_complex_double* x, int64_t incx, double* scale, double* sumsq );

int64_t LAPACKE_slaswp_work_64( int matrix_layout, int64_t n, float* a,
                                int64_t lda, int64_t k1, int64_t k2,
                                const int64_t* ipiv, int64_t incx );
int64_t LAPACKE_dlaswp_work_64( int matrix_layout, int64_t n, double* a,
                                int64_t lda, int64_t k1, int64_t k2,
                                const int64_t* ipiv, int64_t incx );
int64_t LAPACKE_claswp_work_64( int matrix_layout, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t k1, int64_t k2,
                                const int64_t* ipiv, int64_t incx );
int64_t LAPACKE_zlaswp_work_64( int matrix_layout, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t k1, int64_t k2,
                                const int64_t* ipiv, int64_t incx );

int64_t LAPACKE_slatms_work_64( int matrix_layout, int64_t m, int64_t n,
                                char dist, int64_t* iseed, char sym,
                                float* d, int64_t mode, float cond,
                                float dmax, int64_t kl, int64_t ku,
                                char pack, float* a, int64_t lda,
                                float* work );
int64_t LAPACKE_dlatms_work_64( int matrix_layout, int64_t m, int64_t n,
                                char dist, int64_t* iseed, char sym,
                                double* d, int64_t mode, double cond,
                                double dmax, int64_t kl, int64_t ku,
                                char pack, double* a, int64_t lda,
                                double* work );
int64_t LAPACKE_clatms_work_64( int matrix_layout, int64_t m, int64_t n,
                                char dist, int64_t* iseed, char sym,
                                float* d, int64_t mode, float cond,
                                float dmax, int64_t kl, int64_t ku,
                                char pack, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* work );
int64_t LAPACKE_zlatms_work_64( int matrix_layout, int64_t m, int64_t n,
                                char dist, int64_t* iseed, char sym,
                                double* d, int64_t mode, double cond,
                                double dmax, int64_t kl, int64_t ku,
                                char pack, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* work );

int64_t LAPACKE_slauum_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda );
int64_t LAPACKE_dlauum_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda );
int64_t LAPACKE_clauum_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zlauum_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_sopgtr_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* ap, const float* tau, float* q,
                                int64_t ldq, float* work );
int64_t LAPACKE_dopgtr_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* ap, const double* tau, double* q,
                                int64_t ldq, double* work );

int64_t LAPACKE_sopmtr_work_64( int matrix_layout, char side, char uplo,
                                char trans, int64_t m, int64_t n,
                                const float* ap, const float* tau, float* c,
                                int64_t ldc, float* work );
int64_t LAPACKE_dopmtr_work_64( int matrix_layout, char side, char uplo,
                                char trans, int64_t m, int64_t n,
                                const double* ap, const double* tau, double* c,
                                int64_t ldc, double* work );

int64_t LAPACKE_sorgbr_work_64( int matrix_layout, char vect, int64_t m,
                                int64_t n, int64_t k, float* a,
                                int64_t lda, const float* tau, float* work,
                                int64_t lwork );
int64_t LAPACKE_dorgbr_work_64( int matrix_layout, char vect, int64_t m,
                                int64_t n, int64_t k, double* a,
                                int64_t lda, const double* tau, double* work,
                                int64_t lwork );

int64_t LAPACKE_sorghr_work_64( int matrix_layout, int64_t n, int64_t ilo,
                                int64_t ihi, float* a, int64_t lda,
                                const float* tau, float* work,
                                int64_t lwork );
int64_t LAPACKE_dorghr_work_64( int matrix_layout, int64_t n, int64_t ilo,
                                int64_t ihi, double* a, int64_t lda,
                                const double* tau, double* work,
                                int64_t lwork );

int64_t LAPACKE_sorglq_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, float* a, int64_t lda,
                                const float* tau, float* work,
                                int64_t lwork );
int64_t LAPACKE_dorglq_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, double* a, int64_t lda,
                                const double* tau, double* work,
                                int64_t lwork );

int64_t LAPACKE_sorgql_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, float* a, int64_t lda,
                                const float* tau, float* work,
                                int64_t lwork );
int64_t LAPACKE_dorgql_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, double* a, int64_t lda,
                                const double* tau, double* work,
                                int64_t lwork );

int64_t LAPACKE_sorgqr_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, float* a, int64_t lda,
                                const float* tau, float* work,
                                int64_t lwork );
int64_t LAPACKE_dorgqr_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, double* a, int64_t lda,
                                const double* tau, double* work,
                                int64_t lwork );

int64_t LAPACKE_sorgrq_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, float* a, int64_t lda,
                                const float* tau, float* work,
                                int64_t lwork );
int64_t LAPACKE_dorgrq_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, double* a, int64_t lda,
                                const double* tau, double* work,
                                int64_t lwork );

int64_t LAPACKE_sorgtr_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda, const float* tau,
                                float* work, int64_t lwork );
int64_t LAPACKE_dorgtr_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda, const double* tau,
                                double* work, int64_t lwork );

int64_t LAPACKE_sorgtsqr_row_work_64( int matrix_layout,
                                      int64_t m, int64_t n,
                                      int64_t mb, int64_t nb,
                                      float* a, int64_t lda,
                                      const float* t, int64_t ldt,
                                      float* work, int64_t lwork );
int64_t LAPACKE_dorgtsqr_row_work_64( int matrix_layout,
                                      int64_t m, int64_t n,
                                      int64_t mb, int64_t nb,
                                      double* a, int64_t lda,
                                      const double* t, int64_t ldt,
                                      double* work, int64_t lwork );

int64_t LAPACKE_sormbr_work_64( int matrix_layout, char vect, char side,
                                char trans, int64_t m, int64_t n,
                                int64_t k, const float* a, int64_t lda,
                                const float* tau, float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dormbr_work_64( int matrix_layout, char vect, char side,
                                char trans, int64_t m, int64_t n,
                                int64_t k, const double* a, int64_t lda,
                                const double* tau, double* c, int64_t ldc,
                                double* work, int64_t lwork );

int64_t LAPACKE_sormhr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t ilo,
                                int64_t ihi, const float* a, int64_t lda,
                                const float* tau, float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dormhr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t ilo,
                                int64_t ihi, const double* a, int64_t lda,
                                const double* tau, double* c, int64_t ldc,
                                double* work, int64_t lwork );

int64_t LAPACKE_sormlq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const float* a, int64_t lda,
                                const float* tau, float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dormlq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const double* a, int64_t lda,
                                const double* tau, double* c, int64_t ldc,
                                double* work, int64_t lwork );

int64_t LAPACKE_sormql_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const float* a, int64_t lda,
                                const float* tau, float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dormql_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const double* a, int64_t lda,
                                const double* tau, double* c, int64_t ldc,
                                double* work, int64_t lwork );

int64_t LAPACKE_sormqr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const float* a, int64_t lda,
                                const float* tau, float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dormqr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const double* a, int64_t lda,
                                const double* tau, double* c, int64_t ldc,
                                double* work, int64_t lwork );

int64_t LAPACKE_sormrq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const float* a, int64_t lda,
                                const float* tau, float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dormrq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const double* a, int64_t lda,
                                const double* tau, double* c, int64_t ldc,
                                double* work, int64_t lwork );

int64_t LAPACKE_sormrz_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                int64_t l, const float* a, int64_t lda,
                                const float* tau, float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dormrz_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                int64_t l, const double* a, int64_t lda,
                                const double* tau, double* c, int64_t ldc,
                                double* work, int64_t lwork );

int64_t LAPACKE_sormtr_work_64( int matrix_layout, char side, char uplo,
                                char trans, int64_t m, int64_t n,
                                const float* a, int64_t lda,
                                const float* tau, float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dormtr_work_64( int matrix_layout, char side, char uplo,
                                char trans, int64_t m, int64_t n,
                                const double* a, int64_t lda,
                                const double* tau, double* c, int64_t ldc,
                                double* work, int64_t lwork );

int64_t LAPACKE_spbcon_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, const float* ab, int64_t ldab,
                                float anorm, float* rcond, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dpbcon_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, const double* ab,
                                int64_t ldab, double anorm, double* rcond,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cpbcon_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, const lapack_complex_float* ab,
                                int64_t ldab, float anorm, float* rcond,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zpbcon_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, const lapack_complex_double* ab,
                                int64_t ldab, double anorm, double* rcond,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_spbequ_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, const float* ab, int64_t ldab,
                                float* s, float* scond, float* amax );
int64_t LAPACKE_dpbequ_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, const double* ab,
                                int64_t ldab, double* s, double* scond,
                                double* amax );
int64_t LAPACKE_cpbequ_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, const lapack_complex_float* ab,
                                int64_t ldab, float* s, float* scond,
                                float* amax );
int64_t LAPACKE_zpbequ_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, const lapack_complex_double* ab,
                                int64_t ldab, double* s, double* scond,
                                double* amax );

int64_t LAPACKE_spbrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, int64_t nrhs, const float* ab,
                                int64_t ldab, const float* afb,
                                int64_t ldafb, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* ferr, float* berr, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dpbrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, int64_t nrhs,
                                const double* ab, int64_t ldab,
                                const double* afb, int64_t ldafb,
                                const double* b, int64_t ldb, double* x,
                                int64_t ldx, double* ferr, double* berr,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cpbrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, int64_t nrhs,
                                const lapack_complex_float* ab, int64_t ldab,
                                const lapack_complex_float* afb,
                                int64_t ldafb, const lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* x,
                                int64_t ldx, float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zpbrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, int64_t nrhs,
                                const lapack_complex_double* ab,
                                int64_t ldab,
                                const lapack_complex_double* afb,
                                int64_t ldafb,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_spbstf_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kb, float* bb, int64_t ldbb );
int64_t LAPACKE_dpbstf_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kb, double* bb, int64_t ldbb );
int64_t LAPACKE_cpbstf_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kb, lapack_complex_float* bb,
                                int64_t ldbb );
int64_t LAPACKE_zpbstf_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kb, lapack_complex_double* bb,
                                int64_t ldbb );

int64_t LAPACKE_spbsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t kd, int64_t nrhs, float* ab,
                               int64_t ldab, float* b, int64_t ldb );
int64_t LAPACKE_dpbsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t kd, int64_t nrhs, double* ab,
                               int64_t ldab, double* b, int64_t ldb );
int64_t LAPACKE_cpbsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t kd, int64_t nrhs,
                               lapack_complex_float* ab, int64_t ldab,
                               lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpbsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t kd, int64_t nrhs,
                               lapack_complex_double* ab, int64_t ldab,
                               lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_spbsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t kd, int64_t nrhs,
                                float* ab, int64_t ldab, float* afb,
                                int64_t ldafb, char* equed, float* s,
                                float* b, int64_t ldb, float* x,
                                int64_t ldx, float* rcond, float* ferr,
                                float* berr, float* work, int64_t* iwork );
int64_t LAPACKE_dpbsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t kd, int64_t nrhs,
                                double* ab, int64_t ldab, double* afb,
                                int64_t ldafb, char* equed, double* s,
                                double* b, int64_t ldb, double* x,
                                int64_t ldx, double* rcond, double* ferr,
                                double* berr, double* work, int64_t* iwork );
int64_t LAPACKE_cpbsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t kd, int64_t nrhs,
                                lapack_complex_float* ab, int64_t ldab,
                                lapack_complex_float* afb, int64_t ldafb,
                                char* equed, float* s, lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* x,
                                int64_t ldx, float* rcond, float* ferr,
                                float* berr, lapack_complex_float* work,
                                float* rwork );
int64_t LAPACKE_zpbsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t kd, int64_t nrhs,
                                lapack_complex_double* ab, int64_t ldab,
                                lapack_complex_double* afb, int64_t ldafb,
                                char* equed, double* s,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_spbtrf_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, float* ab, int64_t ldab );
int64_t LAPACKE_dpbtrf_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, double* ab, int64_t ldab );
int64_t LAPACKE_cpbtrf_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, lapack_complex_float* ab,
                                int64_t ldab );
int64_t LAPACKE_zpbtrf_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, lapack_complex_double* ab,
                                int64_t ldab );

int64_t LAPACKE_spbtrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, int64_t nrhs, const float* ab,
                                int64_t ldab, float* b, int64_t ldb );
int64_t LAPACKE_dpbtrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, int64_t nrhs,
                                const double* ab, int64_t ldab, double* b,
                                int64_t ldb );
int64_t LAPACKE_cpbtrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, int64_t nrhs,
                                const lapack_complex_float* ab, int64_t ldab,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpbtrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t kd, int64_t nrhs,
                                const lapack_complex_double* ab,
                                int64_t ldab, lapack_complex_double* b,
                                int64_t ldb );

int64_t LAPACKE_spftrf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, float* a );
int64_t LAPACKE_dpftrf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, double* a );
int64_t LAPACKE_cpftrf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, lapack_complex_float* a );
int64_t LAPACKE_zpftrf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, lapack_complex_double* a );

int64_t LAPACKE_spftri_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, float* a );
int64_t LAPACKE_dpftri_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, double* a );
int64_t LAPACKE_cpftri_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, lapack_complex_float* a );
int64_t LAPACKE_zpftri_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, lapack_complex_double* a );

int64_t LAPACKE_spftrs_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, int64_t nrhs, const float* a,
                                float* b, int64_t ldb );
int64_t LAPACKE_dpftrs_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, int64_t nrhs, const double* a,
                                double* b, int64_t ldb );
int64_t LAPACKE_cpftrs_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_float* a,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpftrs_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_double* a,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_spocon_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* a, int64_t lda, float anorm,
                                float* rcond, float* work, int64_t* iwork );
int64_t LAPACKE_dpocon_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* a, int64_t lda, double anorm,
                                double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_cpocon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                float anorm, float* rcond,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zpocon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                double anorm, double* rcond,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_spoequ_work_64( int matrix_layout, int64_t n, const float* a,
                                int64_t lda, float* s, float* scond,
                                float* amax );
int64_t LAPACKE_dpoequ_work_64( int matrix_layout, int64_t n, const double* a,
                                int64_t lda, double* s, double* scond,
                                double* amax );
int64_t LAPACKE_cpoequ_work_64( int matrix_layout, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                float* s, float* scond, float* amax );
int64_t LAPACKE_zpoequ_work_64( int matrix_layout, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                double* s, double* scond, double* amax );

int64_t LAPACKE_spoequb_work_64( int matrix_layout, int64_t n, const float* a,
                                 int64_t lda, float* s, float* scond,
                                 float* amax );
int64_t LAPACKE_dpoequb_work_64( int matrix_layout, int64_t n,
                                 const double* a, int64_t lda, double* s,
                                 double* scond, double* amax );
int64_t LAPACKE_cpoequb_work_64( int matrix_layout, int64_t n,
                                 const lapack_complex_float* a, int64_t lda,
                                 float* s, float* scond, float* amax );
int64_t LAPACKE_zpoequb_work_64( int matrix_layout, int64_t n,
                                 const lapack_complex_double* a, int64_t lda,
                                 double* s, double* scond, double* amax );

int64_t LAPACKE_sporfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const float* af, int64_t ldaf,
                                const float* b, int64_t ldb, float* x,
                                int64_t ldx, float* ferr, float* berr,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dporfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* a,
                                int64_t lda, const double* af,
                                int64_t ldaf, const double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* ferr, double* berr, double* work,
                                int64_t* iwork );
int64_t LAPACKE_cporfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* af,
                                int64_t ldaf, const lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* x,
                                int64_t ldx, float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zporfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, const lapack_complex_double* af,
                                int64_t ldaf, const lapack_complex_double* b,
                                int64_t ldb, lapack_complex_double* x,
                                int64_t ldx, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sporfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs, const float* a,
                                 int64_t lda, const float* af,
                                 int64_t ldaf, const float* s,
                                 const float* b, int64_t ldb, float* x,
                                 int64_t ldx, float* rcond, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, float* work,
                                 int64_t* iwork );
int64_t LAPACKE_dporfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs, const double* a,
                                 int64_t lda, const double* af,
                                 int64_t ldaf, const double* s,
                                 const double* b, int64_t ldb, double* x,
                                 int64_t ldx, double* rcond, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, double* work,
                                 int64_t* iwork );
int64_t LAPACKE_cporfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs,
                                 const lapack_complex_float* a, int64_t lda,
                                 const lapack_complex_float* af,
                                 int64_t ldaf, const float* s,
                                 const lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* x, int64_t ldx,
                                 float* rcond, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, lapack_complex_float* work,
                                 float* rwork );
int64_t LAPACKE_zporfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs,
                                 const lapack_complex_double* a, int64_t lda,
                                 const lapack_complex_double* af,
                                 int64_t ldaf, const double* s,
                                 const lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* x, int64_t ldx,
                                 double* rcond, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, lapack_complex_double* work,
                                 double* rwork );

int64_t LAPACKE_sposv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* a, int64_t lda,
                               float* b, int64_t ldb );
int64_t LAPACKE_dposv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* a, int64_t lda,
                               double* b, int64_t ldb );
int64_t LAPACKE_cposv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* b,
                               int64_t ldb );
int64_t LAPACKE_zposv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* b,
                               int64_t ldb );
int64_t LAPACKE_dsposv_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, double* a, int64_t lda,
                                double* b, int64_t ldb, double* x,
                                int64_t ldx, double* work, float* swork,
                                int64_t* iter );
int64_t LAPACKE_zcposv_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb, lapack_complex_double* x,
                                int64_t ldx, lapack_complex_double* work,
                                lapack_complex_float* swork, double* rwork,
                                int64_t* iter );

int64_t LAPACKE_sposvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, float* a,
                                int64_t lda, float* af, int64_t ldaf,
                                char* equed, float* s, float* b, int64_t ldb,
                                float* x, int64_t ldx, float* rcond,
                                float* ferr, float* berr, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dposvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, double* a,
                                int64_t lda, double* af, int64_t ldaf,
                                char* equed, double* s, double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cposvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* af, int64_t ldaf,
                                char* equed, float* s, lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* x,
                                int64_t ldx, float* rcond, float* ferr,
                                float* berr, lapack_complex_float* work,
                                float* rwork );
int64_t LAPACKE_zposvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* af, int64_t ldaf,
                                char* equed, double* s,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sposvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs, float* a,
                                 int64_t lda, float* af, int64_t ldaf,
                                 char* equed, float* s, float* b,
                                 int64_t ldb, float* x, int64_t ldx,
                                 float* rcond, float* rpvgrw, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, float* work,
                                 int64_t* iwork );
int64_t LAPACKE_dposvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs, double* a,
                                 int64_t lda, double* af, int64_t ldaf,
                                 char* equed, double* s, double* b,
                                 int64_t ldb, double* x, int64_t ldx,
                                 double* rcond, double* rpvgrw, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, double* work,
                                 int64_t* iwork );
int64_t LAPACKE_cposvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* af, int64_t ldaf,
                                 char* equed, float* s, lapack_complex_float* b,
                                 int64_t ldb, lapack_complex_float* x,
                                 int64_t ldx, float* rcond, float* rpvgrw,
                                 float* berr, int64_t n_err_bnds,
                                 float* err_bnds_norm, float* err_bnds_comp,
                                 int64_t nparams, float* params,
                                 lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zposvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* af, int64_t ldaf,
                                 char* equed, double* s,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* x, int64_t ldx,
                                 double* rcond, double* rpvgrw, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, lapack_complex_double* work,
                                 double* rwork );

int64_t LAPACKE_spotrf2_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda );
int64_t LAPACKE_dpotrf2_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda );
int64_t LAPACKE_cpotrf2_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zpotrf2_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_spotrf_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda );
int64_t LAPACKE_dpotrf_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda );
int64_t LAPACKE_cpotrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zpotrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_spotri_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda );
int64_t LAPACKE_dpotri_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda );
int64_t LAPACKE_cpotri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zpotri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_spotrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                float* b, int64_t ldb );
int64_t LAPACKE_dpotrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* a,
                                int64_t lda, double* b, int64_t ldb );
int64_t LAPACKE_cpotrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* b,
                                int64_t ldb );
int64_t LAPACKE_zpotrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* b,
                                int64_t ldb );

int64_t LAPACKE_sppcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* ap, float anorm, float* rcond,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dppcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* ap, double anorm, double* rcond,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cppcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* ap, float anorm,
                                float* rcond, lapack_complex_float* work,
                                float* rwork );
int64_t LAPACKE_zppcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* ap, double anorm,
                                double* rcond, lapack_complex_double* work,
                                double* rwork );

int64_t LAPACKE_sppequ_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* ap, float* s, float* scond,
                                float* amax );
int64_t LAPACKE_dppequ_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* ap, double* s, double* scond,
                                double* amax );
int64_t LAPACKE_cppequ_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* ap, float* s,
                                float* scond, float* amax );
int64_t LAPACKE_zppequ_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* ap, double* s,
                                double* scond, double* amax );

int64_t LAPACKE_spprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* ap,
                                const float* afp, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* ferr, float* berr, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dpprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* ap,
                                const double* afp, const double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* ferr, double* berr, double* work,
                                int64_t* iwork );
int64_t LAPACKE_cpprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* ap,
                                const lapack_complex_float* afp,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zpprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs,
                                const lapack_complex_double* ap,
                                const lapack_complex_double* afp,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sppsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* ap, float* b,
                               int64_t ldb );
int64_t LAPACKE_dppsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* ap, double* b,
                               int64_t ldb );
int64_t LAPACKE_cppsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* ap,
                               lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zppsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* ap,
                               lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sppsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, float* ap,
                                float* afp, char* equed, float* s, float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dppsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, double* ap,
                                double* afp, char* equed, double* s, double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cppsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                lapack_complex_float* ap,
                                lapack_complex_float* afp, char* equed,
                                float* s, lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* x,
                                int64_t ldx, float* rcond, float* ferr,
                                float* berr, lapack_complex_float* work,
                                float* rwork );
int64_t LAPACKE_zppsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                lapack_complex_double* ap,
                                lapack_complex_double* afp, char* equed,
                                double* s, lapack_complex_double* b,
                                int64_t ldb, lapack_complex_double* x,
                                int64_t ldx, double* rcond, double* ferr,
                                double* berr, lapack_complex_double* work,
                                double* rwork );

int64_t LAPACKE_spptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                float* ap );
int64_t LAPACKE_dpptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                double* ap );
int64_t LAPACKE_cpptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* ap );
int64_t LAPACKE_zpptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* ap );

int64_t LAPACKE_spptri_work_64( int matrix_layout, char uplo, int64_t n,
                                float* ap );
int64_t LAPACKE_dpptri_work_64( int matrix_layout, char uplo, int64_t n,
                                double* ap );
int64_t LAPACKE_cpptri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* ap );
int64_t LAPACKE_zpptri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* ap );

int64_t LAPACKE_spptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* ap, float* b,
                                int64_t ldb );
int64_t LAPACKE_dpptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* ap, double* b,
                                int64_t ldb );
int64_t LAPACKE_cpptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* ap,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs,
                                const lapack_complex_double* ap,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_spstrf_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda, int64_t* piv,
                                int64_t* rank, float tol, float* work );
int64_t LAPACKE_dpstrf_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda, int64_t* piv,
                                int64_t* rank, double tol, double* work );
int64_t LAPACKE_cpstrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* piv, int64_t* rank, float tol,
                                float* work );
int64_t LAPACKE_zpstrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* piv, int64_t* rank, double tol,
                                double* work );

int64_t LAPACKE_sptcon_work_64( int64_t n, const float* d, const float* e,
                                float anorm, float* rcond, float* work );
int64_t LAPACKE_dptcon_work_64( int64_t n, const double* d, const double* e,
                                double anorm, double* rcond, double* work );
int64_t LAPACKE_cptcon_work_64( int64_t n, const float* d,
                                const lapack_complex_float* e, float anorm,
                                float* rcond, float* work );
int64_t LAPACKE_zptcon_work_64( int64_t n, const double* d,
                                const lapack_complex_double* e, double anorm,
                                double* rcond, double* work );

int64_t LAPACKE_spteqr_work_64( int matrix_layout, char compz, int64_t n,
                                float* d, float* e, float* z, int64_t ldz,
                                float* work );
int64_t LAPACKE_dpteqr_work_64( int matrix_layout, char compz, int64_t n,
                                double* d, double* e, double* z, int64_t ldz,
                                double* work );
int64_t LAPACKE_cpteqr_work_64( int matrix_layout, char compz, int64_t n,
                                float* d, float* e, lapack_complex_float* z,
                                int64_t ldz, float* work );
int64_t LAPACKE_zpteqr_work_64( int matrix_layout, char compz, int64_t n,
                                double* d, double* e, lapack_complex_double* z,
                                int64_t ldz, double* work );

int64_t LAPACKE_sptrfs_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                                const float* d, const float* e, const float* df,
                                const float* ef, const float* b, int64_t ldb,
                                float* x, int64_t ldx, float* ferr,
                                float* berr, float* work );
int64_t LAPACKE_dptrfs_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                                const double* d, const double* e,
                                const double* df, const double* ef,
                                const double* b, int64_t ldb, double* x,
                                int64_t ldx, double* ferr, double* berr,
                                double* work );
int64_t LAPACKE_cptrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* d,
                                const lapack_complex_float* e, const float* df,
                                const lapack_complex_float* ef,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zptrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* d,
                                const lapack_complex_double* e,
                                const double* df,
                                const lapack_complex_double* ef,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sptsv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               float* d, float* e, float* b, int64_t ldb );
int64_t LAPACKE_dptsv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               double* d, double* e, double* b,
                               int64_t ldb );
int64_t LAPACKE_cptsv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               float* d, lapack_complex_float* e,
                               lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zptsv_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                               double* d, lapack_complex_double* e,
                               lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sptsvx_work_64( int matrix_layout, char fact, int64_t n,
                                int64_t nrhs, const float* d, const float* e,
                                float* df, float* ef, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                float* work );
int64_t LAPACKE_dptsvx_work_64( int matrix_layout, char fact, int64_t n,
                                int64_t nrhs, const double* d,
                                const double* e, double* df, double* ef,
                                const double* b, int64_t ldb, double* x,
                                int64_t ldx, double* rcond, double* ferr,
                                double* berr, double* work );
int64_t LAPACKE_cptsvx_work_64( int matrix_layout, char fact, int64_t n,
                                int64_t nrhs, const float* d,
                                const lapack_complex_float* e, float* df,
                                lapack_complex_float* ef,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zptsvx_work_64( int matrix_layout, char fact, int64_t n,
                                int64_t nrhs, const double* d,
                                const lapack_complex_double* e, double* df,
                                lapack_complex_double* ef,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_spttrf_work_64( int64_t n, float* d, float* e );
int64_t LAPACKE_dpttrf_work_64( int64_t n, double* d, double* e );
int64_t LAPACKE_cpttrf_work_64( int64_t n, float* d,
                                lapack_complex_float* e );
int64_t LAPACKE_zpttrf_work_64( int64_t n, double* d,
                                lapack_complex_double* e );

int64_t LAPACKE_spttrs_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                                const float* d, const float* e, float* b,
                                int64_t ldb );
int64_t LAPACKE_dpttrs_work_64( int matrix_layout, int64_t n, int64_t nrhs,
                                const double* d, const double* e, double* b,
                                int64_t ldb );
int64_t LAPACKE_cpttrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* d,
                                const lapack_complex_float* e,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zpttrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* d,
                                const lapack_complex_double* e,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_ssbev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t kd, float* ab,
                               int64_t ldab, float* w, float* z,
                               int64_t ldz, float* work );
int64_t LAPACKE_dsbev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t kd, double* ab,
                               int64_t ldab, double* w, double* z,
                               int64_t ldz, double* work );

int64_t LAPACKE_ssbevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t kd, float* ab,
                                int64_t ldab, float* w, float* z,
                                int64_t ldz, float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dsbevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t kd, double* ab,
                                int64_t ldab, double* w, double* z,
                                int64_t ldz, double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_ssbevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t kd,
                                float* ab, int64_t ldab, float* q,
                                int64_t ldq, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w, float* z,
                                int64_t ldz, float* work,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_dsbevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t kd,
                                double* ab, int64_t ldab, double* q,
                                int64_t ldq, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w, double* z,
                                int64_t ldz, double* work,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_ssbgst_work_64( int matrix_layout, char vect, char uplo,
                                int64_t n, int64_t ka, int64_t kb,
                                float* ab, int64_t ldab, const float* bb,
                                int64_t ldbb, float* x, int64_t ldx,
                                float* work );
int64_t LAPACKE_dsbgst_work_64( int matrix_layout, char vect, char uplo,
                                int64_t n, int64_t ka, int64_t kb,
                                double* ab, int64_t ldab, const double* bb,
                                int64_t ldbb, double* x, int64_t ldx,
                                double* work );

int64_t LAPACKE_ssbgv_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t ka, int64_t kb,
                               float* ab, int64_t ldab, float* bb,
                               int64_t ldbb, float* w, float* z,
                               int64_t ldz, float* work );
int64_t LAPACKE_dsbgv_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t ka, int64_t kb,
                               double* ab, int64_t ldab, double* bb,
                               int64_t ldbb, double* w, double* z,
                               int64_t ldz, double* work );

int64_t LAPACKE_ssbgvd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t ka, int64_t kb,
                                float* ab, int64_t ldab, float* bb,
                                int64_t ldbb, float* w, float* z,
                                int64_t ldz, float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dsbgvd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t ka, int64_t kb,
                                double* ab, int64_t ldab, double* bb,
                                int64_t ldbb, double* w, double* z,
                                int64_t ldz, double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_ssbgvx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t ka,
                                int64_t kb, float* ab, int64_t ldab,
                                float* bb, int64_t ldbb, float* q,
                                int64_t ldq, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w, float* z,
                                int64_t ldz, float* work, int64_t* iwork,
                                int64_t* ifail );
int64_t LAPACKE_dsbgvx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t ka,
                                int64_t kb, double* ab, int64_t ldab,
                                double* bb, int64_t ldbb, double* q,
                                int64_t ldq, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w, double* z,
                                int64_t ldz, double* work, int64_t* iwork,
                                int64_t* ifail );

int64_t LAPACKE_ssbtrd_work_64( int matrix_layout, char vect, char uplo,
                                int64_t n, int64_t kd, float* ab,
                                int64_t ldab, float* d, float* e, float* q,
                                int64_t ldq, float* work );
int64_t LAPACKE_dsbtrd_work_64( int matrix_layout, char vect, char uplo,
                                int64_t n, int64_t kd, double* ab,
                                int64_t ldab, double* d, double* e,
                                double* q, int64_t ldq, double* work );

int64_t LAPACKE_ssfrk_work_64( int matrix_layout, char transr, char uplo,
                               char trans, int64_t n, int64_t k,
                               float alpha, const float* a, int64_t lda,
                               float beta, float* c );
int64_t LAPACKE_dsfrk_work_64( int matrix_layout, char transr, char uplo,
                               char trans, int64_t n, int64_t k,
                               double alpha, const double* a, int64_t lda,
                               double beta, double* c );

int64_t LAPACKE_sspcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* ap, const int64_t* ipiv,
                                float anorm, float* rcond, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dspcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* ap, const int64_t* ipiv,
                                double anorm, double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_cspcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* ap,
                                const int64_t* ipiv, float anorm,
                                float* rcond, lapack_complex_float* work );
int64_t LAPACKE_zspcon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* ap,
                                const int64_t* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

int64_t LAPACKE_sspev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, float* ap, float* w, float* z,
                               int64_t ldz, float* work );
int64_t LAPACKE_dspev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, double* ap, double* w, double* z,
                               int64_t ldz, double* work );

int64_t LAPACKE_sspevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, float* ap, float* w, float* z,
                                int64_t ldz, float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dspevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, double* ap, double* w, double* z,
                                int64_t ldz, double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_sspevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, float* ap, float vl,
                                float vu, int64_t il, int64_t iu,
                                float abstol, int64_t* m, float* w, float* z,
                                int64_t ldz, float* work, int64_t* iwork,
                                int64_t* ifail );
int64_t LAPACKE_dspevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, double* ap, double vl,
                                double vu, int64_t il, int64_t iu,
                                double abstol, int64_t* m, double* w,
                                double* z, int64_t ldz, double* work,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_sspgst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, float* ap, const float* bp );
int64_t LAPACKE_dspgst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, double* ap, const double* bp );

int64_t LAPACKE_sspgv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, float* ap, float* bp,
                               float* w, float* z, int64_t ldz,
                               float* work );
int64_t LAPACKE_dspgv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, double* ap, double* bp,
                               double* w, double* z, int64_t ldz,
                               double* work );

int64_t LAPACKE_sspgvd_work_64( int matrix_layout, int64_t itype, char jobz,
                                char uplo, int64_t n, float* ap, float* bp,
                                float* w, float* z, int64_t ldz, float* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_dspgvd_work_64( int matrix_layout, int64_t itype, char jobz,
                                char uplo, int64_t n, double* ap, double* bp,
                                double* w, double* z, int64_t ldz,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_sspgvx_work_64( int matrix_layout, int64_t itype, char jobz,
                                char range, char uplo, int64_t n, float* ap,
                                float* bp, float vl, float vu, int64_t il,
                                int64_t iu, float abstol, int64_t* m,
                                float* w, float* z, int64_t ldz, float* work,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_dspgvx_work_64( int matrix_layout, int64_t itype, char jobz,
                                char range, char uplo, int64_t n, double* ap,
                                double* bp, double vl, double vu, int64_t il,
                                int64_t iu, double abstol, int64_t* m,
                                double* w, double* z, int64_t ldz,
                                double* work, int64_t* iwork,
                                int64_t* ifail );

int64_t LAPACKE_ssprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* ap,
                                const float* afp, const int64_t* ipiv,
                                const float* b, int64_t ldb, float* x,
                                int64_t ldx, float* ferr, float* berr,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dsprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* ap,
                                const double* afp, const int64_t* ipiv,
                                const double* b, int64_t ldb, double* x,
                                int64_t ldx, double* ferr, double* berr,
                                double* work, int64_t* iwork );
int64_t LAPACKE_csprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* ap,
                                const lapack_complex_float* afp,
                                const int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zsprfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs,
                                const lapack_complex_double* ap,
                                const lapack_complex_double* afp,
                                const int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_sspsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* ap, int64_t* ipiv,
                               float* b, int64_t ldb );
int64_t LAPACKE_dspsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* ap, int64_t* ipiv,
                               double* b, int64_t ldb );
int64_t LAPACKE_cspsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* ap,
                               int64_t* ipiv, lapack_complex_float* b,
                               int64_t ldb );
int64_t LAPACKE_zspsv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* ap,
                               int64_t* ipiv, lapack_complex_double* b,
                               int64_t ldb );

int64_t LAPACKE_sspsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, const float* ap,
                                float* afp, int64_t* ipiv, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dspsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, const double* ap,
                                double* afp, int64_t* ipiv, const double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, int64_t* iwork );
int64_t LAPACKE_cspsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_float* ap,
                                lapack_complex_float* afp, int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zspsvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_double* ap,
                                lapack_complex_double* afp, int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_ssptrd_work_64( int matrix_layout, char uplo, int64_t n,
                                float* ap, float* d, float* e, float* tau );
int64_t LAPACKE_dsptrd_work_64( int matrix_layout, char uplo, int64_t n,
                                double* ap, double* d, double* e, double* tau );

int64_t LAPACKE_ssptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                float* ap, int64_t* ipiv );
int64_t LAPACKE_dsptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                double* ap, int64_t* ipiv );
int64_t LAPACKE_csptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* ap, int64_t* ipiv );
int64_t LAPACKE_zsptrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* ap, int64_t* ipiv );

int64_t LAPACKE_ssptri_work_64( int matrix_layout, char uplo, int64_t n,
                                float* ap, const int64_t* ipiv,
                                float* work );
int64_t LAPACKE_dsptri_work_64( int matrix_layout, char uplo, int64_t n,
                                double* ap, const int64_t* ipiv,
                                double* work );
int64_t LAPACKE_csptri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* ap,
                                const int64_t* ipiv,
                                lapack_complex_float* work );
int64_t LAPACKE_zsptri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* ap,
                                const int64_t* ipiv,
                                lapack_complex_double* work );

int64_t LAPACKE_ssptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* ap,
                                const int64_t* ipiv, float* b,
                                int64_t ldb );
int64_t LAPACKE_dsptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* ap,
                                const int64_t* ipiv, double* b,
                                int64_t ldb );
int64_t LAPACKE_csptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* ap,
                                const int64_t* ipiv, lapack_complex_float* b,
                                int64_t ldb );
int64_t LAPACKE_zsptrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs,
                                const lapack_complex_double* ap,
                                const int64_t* ipiv,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sstebz_work_64( char range, char order, int64_t n, float vl,
                                float vu, int64_t il, int64_t iu,
                                float abstol, const float* d, const float* e,
                                int64_t* m, int64_t* nsplit, float* w,
                                int64_t* iblock, int64_t* isplit,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dstebz_work_64( char range, char order, int64_t n, double vl,
                                double vu, int64_t il, int64_t iu,
                                double abstol, const double* d, const double* e,
                                int64_t* m, int64_t* nsplit, double* w,
                                int64_t* iblock, int64_t* isplit,
                                double* work, int64_t* iwork );

int64_t LAPACKE_sstedc_work_64( int matrix_layout, char compz, int64_t n,
                                float* d, float* e, float* z, int64_t ldz,
                                float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dstedc_work_64( int matrix_layout, char compz, int64_t n,
                                double* d, double* e, double* z, int64_t ldz,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_cstedc_work_64( int matrix_layout, char compz, int64_t n,
                                float* d, float* e, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_zstedc_work_64( int matrix_layout, char compz, int64_t n,
                                double* d, double* e, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );

int64_t LAPACKE_sstegr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, float* d, float* e, float vl,
                                float vu, int64_t il, int64_t iu,
                                float abstol, int64_t* m, float* w, float* z,
                                int64_t ldz, int64_t* isuppz, float* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_dstegr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, double* d, double* e, double vl,
                                double vu, int64_t il, int64_t iu,
                                double abstol, int64_t* m, double* w,
                                double* z, int64_t ldz, int64_t* isuppz,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_cstegr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, float* d, float* e, float vl,
                                float vu, int64_t il, int64_t iu,
                                float abstol, int64_t* m, float* w,
                                lapack_complex_float* z, int64_t ldz,
                                int64_t* isuppz, float* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_zstegr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, double* d, double* e, double vl,
                                double vu, int64_t il, int64_t iu,
                                double abstol, int64_t* m, double* w,
                                lapack_complex_double* z, int64_t ldz,
                                int64_t* isuppz, double* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );

int64_t LAPACKE_sstein_work_64( int matrix_layout, int64_t n, const float* d,
                                const float* e, int64_t m, const float* w,
                                const int64_t* iblock,
                                const int64_t* isplit, float* z,
                                int64_t ldz, float* work, int64_t* iwork,
                                int64_t* ifailv );
int64_t LAPACKE_dstein_work_64( int matrix_layout, int64_t n, const double* d,
                                const double* e, int64_t m, const double* w,
                                const int64_t* iblock,
                                const int64_t* isplit, double* z,
                                int64_t ldz, double* work, int64_t* iwork,
                                int64_t* ifailv );
int64_t LAPACKE_cstein_work_64( int matrix_layout, int64_t n, const float* d,
                                const float* e, int64_t m, const float* w,
                                const int64_t* iblock,
                                const int64_t* isplit,
                                lapack_complex_float* z, int64_t ldz,
                                float* work, int64_t* iwork,
                                int64_t* ifailv );
int64_t LAPACKE_zstein_work_64( int matrix_layout, int64_t n, const double* d,
                                const double* e, int64_t m, const double* w,
                                const int64_t* iblock,
                                const int64_t* isplit,
                                lapack_complex_double* z, int64_t ldz,
                                double* work, int64_t* iwork,
                                int64_t* ifailv );

int64_t LAPACKE_sstemr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, float* d, float* e, float vl,
                                float vu, int64_t il, int64_t iu,
                                int64_t* m, float* w, float* z,
                                int64_t ldz, int64_t nzc,
                                int64_t* isuppz, lapack_logical* tryrac,
                                float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dstemr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, double* d, double* e, double vl,
                                double vu, int64_t il, int64_t iu,
                                int64_t* m, double* w, double* z,
                                int64_t ldz, int64_t nzc,
                                int64_t* isuppz, lapack_logical* tryrac,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_cstemr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, float* d, float* e, float vl,
                                float vu, int64_t il, int64_t iu,
                                int64_t* m, float* w,
                                lapack_complex_float* z, int64_t ldz,
                                int64_t nzc, int64_t* isuppz,
                                lapack_logical* tryrac, float* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_zstemr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, double* d, double* e, double vl,
                                double vu, int64_t il, int64_t iu,
                                int64_t* m, double* w,
                                lapack_complex_double* z, int64_t ldz,
                                int64_t nzc, int64_t* isuppz,
                                lapack_logical* tryrac, double* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );

int64_t LAPACKE_ssteqr_work_64( int matrix_layout, char compz, int64_t n,
                                float* d, float* e, float* z, int64_t ldz,
                                float* work );
int64_t LAPACKE_dsteqr_work_64( int matrix_layout, char compz, int64_t n,
                                double* d, double* e, double* z, int64_t ldz,
                                double* work );
int64_t LAPACKE_csteqr_work_64( int matrix_layout, char compz, int64_t n,
                                float* d, float* e, lapack_complex_float* z,
                                int64_t ldz, float* work );
int64_t LAPACKE_zsteqr_work_64( int matrix_layout, char compz, int64_t n,
                                double* d, double* e, lapack_complex_double* z,
                                int64_t ldz, double* work );

int64_t LAPACKE_skteqr_work_64( int matrix_layout, char compz, int64_t n,
                                float* e, float* z, int64_t ldz,
                                float* work );
int64_t LAPACKE_dkteqr_work_64( int matrix_layout, char compz, int64_t n,
                                double* e, double* z, int64_t ldz,
                                double* work );

int64_t LAPACKE_ssterf_work_64( int64_t n, float* d, float* e );
int64_t LAPACKE_dsterf_work_64( int64_t n, double* d, double* e );

int64_t LAPACKE_sstev_work_64( int matrix_layout, char jobz, int64_t n,
                               float* d, float* e, float* z, int64_t ldz,
                               float* work );
int64_t LAPACKE_dstev_work_64( int matrix_layout, char jobz, int64_t n,
                               double* d, double* e, double* z, int64_t ldz,
                               double* work );

int64_t LAPACKE_sktev_work_64( int matrix_layout, char jobz, int64_t n,
                               float* d, float* e, float* z, int64_t ldz,
                               float* work );
int64_t LAPACKE_dktev_work_64( int matrix_layout, char jobz, int64_t n,
                               double* d, double* e, double* z, int64_t ldz,
                               double* work );

int64_t LAPACKE_sstevd_work_64( int matrix_layout, char jobz, int64_t n,
                                float* d, float* e, float* z, int64_t ldz,
                                float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dstevd_work_64( int matrix_layout, char jobz, int64_t n,
                                double* d, double* e, double* z, int64_t ldz,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_sstevr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, float* d, float* e, float vl,
                                float vu, int64_t il, int64_t iu,
                                float abstol, int64_t* m, float* w, float* z,
                                int64_t ldz, int64_t* isuppz, float* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_dstevr_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, double* d, double* e, double vl,
                                double vu, int64_t il, int64_t iu,
                                double abstol, int64_t* m, double* w,
                                double* z, int64_t ldz, int64_t* isuppz,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_sstevx_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, float* d, float* e, float vl,
                                float vu, int64_t il, int64_t iu,
                                float abstol, int64_t* m, float* w, float* z,
                                int64_t ldz, float* work, int64_t* iwork,
                                int64_t* ifail );
int64_t LAPACKE_dstevx_work_64( int matrix_layout, char jobz, char range,
                                int64_t n, double* d, double* e, double vl,
                                double vu, int64_t il, int64_t iu,
                                double abstol, int64_t* m, double* w,
                                double* z, int64_t ldz, double* work,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_ssycon_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* a, int64_t lda,
                                const int64_t* ipiv, float anorm,
                                float* rcond, float* work, int64_t* iwork );
int64_t LAPACKE_dsycon_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* a, int64_t lda,
                                const int64_t* ipiv, double anorm,
                                double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_csycon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                const int64_t* ipiv, float anorm,
                                float* rcond, lapack_complex_float* work );
int64_t LAPACKE_zsycon_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                const int64_t* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

int64_t LAPACKE_skycon_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* a, int64_t lda,
                                const int64_t* ipiv, float anorm,
                                float* rcond, float* work, int64_t* iwork );
int64_t LAPACKE_dkycon_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* a, int64_t lda,
                                const int64_t* ipiv, double anorm,
                                double* rcond, double* work,
                                int64_t* iwork );

int64_t LAPACKE_ssyequb_work_64( int matrix_layout, char uplo, int64_t n,
                                 const float* a, int64_t lda, float* s,
                                 float* scond, float* amax, float* work );
int64_t LAPACKE_dsyequb_work_64( int matrix_layout, char uplo, int64_t n,
                                 const double* a, int64_t lda, double* s,
                                 double* scond, double* amax, double* work );
int64_t LAPACKE_csyequb_work_64( int matrix_layout, char uplo, int64_t n,
                                 const lapack_complex_float* a, int64_t lda,
                                 float* s, float* scond, float* amax,
                                 lapack_complex_float* work );
int64_t LAPACKE_zsyequb_work_64( int matrix_layout, char uplo, int64_t n,
                                 const lapack_complex_double* a, int64_t lda,
                                 double* s, double* scond, double* amax,
                                 lapack_complex_double* work );

int64_t LAPACKE_ssyev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, float* a, int64_t lda, float* w,
                               float* work, int64_t lwork );
int64_t LAPACKE_dsyev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, double* a, int64_t lda,
                               double* w, double* work, int64_t lwork );

int64_t LAPACKE_skyev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, float* a, int64_t lda, float* w,
                               float* work, int64_t lwork );
int64_t LAPACKE_dkyev_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, double* a, int64_t lda,
                               double* w, double* work, int64_t lwork );

int64_t LAPACKE_ssyevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, float* a, int64_t lda,
                                float* w, float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dsyevd_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, double* a, int64_t lda,
                                double* w, double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_ssyevr_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, float* a,
                                int64_t lda, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w, float* z,
                                int64_t ldz, int64_t* isuppz, float* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_dsyevr_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, double* a,
                                int64_t lda, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w, double* z,
                                int64_t ldz, int64_t* isuppz,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_ssyevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, float* a,
                                int64_t lda, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w, float* z,
                                int64_t ldz, float* work, int64_t lwork,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_dsyevx_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, double* a,
                                int64_t lda, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w, double* z,
                                int64_t ldz, double* work, int64_t lwork,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_ssygst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, float* a, int64_t lda,
                                const float* b, int64_t ldb );
int64_t LAPACKE_dsygst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, double* a, int64_t lda,
                                const double* b, int64_t ldb );

int64_t LAPACKE_skygst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, float* a, int64_t lda,
                                const float* b, int64_t ldb );
int64_t LAPACKE_dkygst_work_64( int matrix_layout, int64_t itype, char uplo,
                                int64_t n, double* a, int64_t lda,
                                const double* b, int64_t ldb );

int64_t LAPACKE_ssygv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, float* a,
                               int64_t lda, float* b, int64_t ldb,
                               float* w, float* work, int64_t lwork );
int64_t LAPACKE_dsygv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, double* a,
                               int64_t lda, double* b, int64_t ldb,
                               double* w, double* work, int64_t lwork );

int64_t LAPACKE_skygv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, float* a,
                               int64_t lda, float* b, int64_t ldb,
                               float* w, float* work, int64_t lwork );
int64_t LAPACKE_dkygv_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, double* a,
                               int64_t lda, double* b, int64_t ldb,
                               double* w, double* work, int64_t lwork );

int64_t LAPACKE_ssygvd_work_64( int matrix_layout, int64_t itype, char jobz,
                                char uplo, int64_t n, float* a,
                                int64_t lda, float* b, int64_t ldb,
                                float* w, float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dsygvd_work_64( int matrix_layout, int64_t itype, char jobz,
                                char uplo, int64_t n, double* a,
                                int64_t lda, double* b, int64_t ldb,
                                double* w, double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_ssygvx_work_64( int matrix_layout, int64_t itype, char jobz,
                                char range, char uplo, int64_t n, float* a,
                                int64_t lda, float* b, int64_t ldb,
                                float vl, float vu, int64_t il,
                                int64_t iu, float abstol, int64_t* m,
                                float* w, float* z, int64_t ldz, float* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t* ifail );
int64_t LAPACKE_dsygvx_work_64( int matrix_layout, int64_t itype, char jobz,
                                char range, char uplo, int64_t n, double* a,
                                int64_t lda, double* b, int64_t ldb,
                                double vl, double vu, int64_t il,
                                int64_t iu, double abstol, int64_t* m,
                                double* w, double* z, int64_t ldz,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_ssyrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const float* af, int64_t ldaf,
                                const int64_t* ipiv, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* ferr, float* berr, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dsyrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* a,
                                int64_t lda, const double* af,
                                int64_t ldaf, const int64_t* ipiv,
                                const double* b, int64_t ldb, double* x,
                                int64_t ldx, double* ferr, double* berr,
                                double* work, int64_t* iwork );
int64_t LAPACKE_csyrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* af,
                                int64_t ldaf, const int64_t* ipiv,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_zsyrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, const lapack_complex_double* af,
                                int64_t ldaf, const int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_skyrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const float* af, int64_t ldaf,
                                const int64_t* ipiv, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* ferr, float* berr, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dkyrfs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* a,
                                int64_t lda, const double* af,
                                int64_t ldaf, const int64_t* ipiv,
                                const double* b, int64_t ldb, double* x,
                                int64_t ldx, double* ferr, double* berr,
                                double* work, int64_t* iwork );

int64_t LAPACKE_ssyrfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs, const float* a,
                                 int64_t lda, const float* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const float* s, const float* b, int64_t ldb,
                                 float* x, int64_t ldx, float* rcond,
                                 float* berr, int64_t n_err_bnds,
                                 float* err_bnds_norm, float* err_bnds_comp,
                                 int64_t nparams, float* params, float* work,
                                 int64_t* iwork );
int64_t LAPACKE_dsyrfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs, const double* a,
                                 int64_t lda, const double* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const double* s, const double* b,
                                 int64_t ldb, double* x, int64_t ldx,
                                 double* rcond, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, double* work,
                                 int64_t* iwork );
int64_t LAPACKE_csyrfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs,
                                 const lapack_complex_float* a, int64_t lda,
                                 const lapack_complex_float* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const float* s, const lapack_complex_float* b,
                                 int64_t ldb, lapack_complex_float* x,
                                 int64_t ldx, float* rcond, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, lapack_complex_float* work,
                                 float* rwork );
int64_t LAPACKE_zsyrfsx_work_64( int matrix_layout, char uplo, char equed,
                                 int64_t n, int64_t nrhs,
                                 const lapack_complex_double* a, int64_t lda,
                                 const lapack_complex_double* af,
                                 int64_t ldaf, const int64_t* ipiv,
                                 const double* s,
                                 const lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* x, int64_t ldx,
                                 double* rcond, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, lapack_complex_double* work,
                                 double* rwork );

int64_t LAPACKE_ssysv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* a, int64_t lda,
                               int64_t* ipiv, float* b, int64_t ldb,
                               float* work, int64_t lwork );
int64_t LAPACKE_dsysv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* a, int64_t lda,
                               int64_t* ipiv, double* b, int64_t ldb,
                               double* work, int64_t lwork );
int64_t LAPACKE_csysv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_float* b, int64_t ldb,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zsysv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_double* b, int64_t ldb,
                               lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_skysv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* a, int64_t lda,
                               int64_t* ipiv, float* b, int64_t ldb,
                               float* work, int64_t lwork );
int64_t LAPACKE_dkysv_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* a, int64_t lda,
                               int64_t* ipiv, double* b, int64_t ldb,
                               double* work, int64_t lwork );

int64_t LAPACKE_ssysvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, const float* a,
                                int64_t lda, float* af, int64_t ldaf,
                                int64_t* ipiv, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                float* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_dsysvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, const double* a,
                                int64_t lda, double* af, int64_t ldaf,
                                int64_t* ipiv, const double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_csysvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* af, int64_t ldaf,
                                int64_t* ipiv, const lapack_complex_float* b,
                                int64_t ldb, lapack_complex_float* x,
                                int64_t ldx, float* rcond, float* ferr,
                                float* berr, lapack_complex_float* work,
                                int64_t lwork, float* rwork );
int64_t LAPACKE_zsysvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs,
                                const lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* af, int64_t ldaf,
                                int64_t* ipiv,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork );

int64_t LAPACKE_skysvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, const float* a,
                                int64_t lda, float* af, int64_t ldaf,
                                int64_t* ipiv, const float* b,
                                int64_t ldb, float* x, int64_t ldx,
                                float* rcond, float* ferr, float* berr,
                                float* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_dkysvx_work_64( int matrix_layout, char fact, char uplo,
                                int64_t n, int64_t nrhs, const double* a,
                                int64_t lda, double* af, int64_t ldaf,
                                int64_t* ipiv, const double* b,
                                int64_t ldb, double* x, int64_t ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, int64_t lwork,
                                int64_t* iwork );

int64_t LAPACKE_ssysvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs, float* a,
                                 int64_t lda, float* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, float* s,
                                 float* b, int64_t ldb, float* x,
                                 int64_t ldx, float* rcond, float* rpvgrw,
                                 float* berr, int64_t n_err_bnds,
                                 float* err_bnds_norm, float* err_bnds_comp,
                                 int64_t nparams, float* params, float* work,
                                 int64_t* iwork );
int64_t LAPACKE_dsysvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs, double* a,
                                 int64_t lda, double* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, double* s,
                                 double* b, int64_t ldb, double* x,
                                 int64_t ldx, double* rcond, double* rpvgrw,
                                 double* berr, int64_t n_err_bnds,
                                 double* err_bnds_norm, double* err_bnds_comp,
                                 int64_t nparams, double* params,
                                 double* work, int64_t* iwork );
int64_t LAPACKE_csysvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, float* s,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* x, int64_t ldx,
                                 float* rcond, float* rpvgrw, float* berr,
                                 int64_t n_err_bnds, float* err_bnds_norm,
                                 float* err_bnds_comp, int64_t nparams,
                                 float* params, lapack_complex_float* work,
                                 float* rwork );
int64_t LAPACKE_zsysvxx_work_64( int matrix_layout, char fact, char uplo,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* af, int64_t ldaf,
                                 int64_t* ipiv, char* equed, double* s,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* x, int64_t ldx,
                                 double* rcond, double* rpvgrw, double* berr,
                                 int64_t n_err_bnds, double* err_bnds_norm,
                                 double* err_bnds_comp, int64_t nparams,
                                 double* params, lapack_complex_double* work,
                                 double* rwork );

int64_t LAPACKE_ssytrd_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda, float* d, float* e,
                                float* tau, float* work, int64_t lwork );
int64_t LAPACKE_dsytrd_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda, double* d, double* e,
                                double* tau, double* work, int64_t lwork );

int64_t LAPACKE_skytrd_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda, float* e,
                                float* tau, float* work, int64_t lwork );
int64_t LAPACKE_dkytrd_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda, double* e,
                                double* tau, double* work, int64_t lwork );

int64_t LAPACKE_ssytrf_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda, int64_t* ipiv,
                                float* work, int64_t lwork );
int64_t LAPACKE_dsytrf_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda, int64_t* ipiv,
                                double* work, int64_t lwork );
int64_t LAPACKE_csytrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_float* work,
                                int64_t lwork );
int64_t LAPACKE_zsytrf_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_double* work,
                                int64_t lwork );

int64_t LAPACKE_skytrf_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda, int64_t* ipiv,
                                float* work, int64_t lwork );
int64_t LAPACKE_dkytrf_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda, int64_t* ipiv,
                                double* work, int64_t lwork );

int64_t LAPACKE_ssytri_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda,
                                const int64_t* ipiv, float* work );
int64_t LAPACKE_dsytri_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda,
                                const int64_t* ipiv, double* work );
int64_t LAPACKE_csytri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                const int64_t* ipiv,
                                lapack_complex_float* work );
int64_t LAPACKE_zsytri_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                const int64_t* ipiv,
                                lapack_complex_double* work );

int64_t LAPACKE_skytri_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda,
                                const int64_t* ipiv, float* work );
int64_t LAPACKE_dkytri_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda,
                                const int64_t* ipiv, double* work );

int64_t LAPACKE_ssytrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const int64_t* ipiv, float* b,
                                int64_t ldb );
int64_t LAPACKE_dsytrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* a,
                                int64_t lda, const int64_t* ipiv,
                                double* b, int64_t ldb );
int64_t LAPACKE_csytrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zsytrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_skytrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const int64_t* ipiv, float* b,
                                int64_t ldb );
int64_t LAPACKE_dkytrs_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* a,
                                int64_t lda, const int64_t* ipiv,
                                double* b, int64_t ldb );

int64_t LAPACKE_stbcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n, int64_t kd,
                                const float* ab, int64_t ldab, float* rcond,
                                float* work, int64_t* iwork );
int64_t LAPACKE_dtbcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n, int64_t kd,
                                const double* ab, int64_t ldab,
                                double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_ctbcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n, int64_t kd,
                                const lapack_complex_float* ab, int64_t ldab,
                                float* rcond, lapack_complex_float* work,
                                float* rwork );
int64_t LAPACKE_ztbcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n, int64_t kd,
                                const lapack_complex_double* ab,
                                int64_t ldab, double* rcond,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_stbrfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t kd,
                                int64_t nrhs, const float* ab,
                                int64_t ldab, const float* b, int64_t ldb,
                                const float* x, int64_t ldx, float* ferr,
                                float* berr, float* work, int64_t* iwork );
int64_t LAPACKE_dtbrfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t kd,
                                int64_t nrhs, const double* ab,
                                int64_t ldab, const double* b,
                                int64_t ldb, const double* x, int64_t ldx,
                                double* ferr, double* berr, double* work,
                                int64_t* iwork );
int64_t LAPACKE_ctbrfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t kd,
                                int64_t nrhs, const lapack_complex_float* ab,
                                int64_t ldab, const lapack_complex_float* b,
                                int64_t ldb, const lapack_complex_float* x,
                                int64_t ldx, float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_ztbrfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t kd,
                                int64_t nrhs,
                                const lapack_complex_double* ab,
                                int64_t ldab, const lapack_complex_double* b,
                                int64_t ldb, const lapack_complex_double* x,
                                int64_t ldx, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_stbtrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t kd,
                                int64_t nrhs, const float* ab,
                                int64_t ldab, float* b, int64_t ldb );
int64_t LAPACKE_dtbtrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t kd,
                                int64_t nrhs, const double* ab,
                                int64_t ldab, double* b, int64_t ldb );
int64_t LAPACKE_ctbtrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t kd,
                                int64_t nrhs, const lapack_complex_float* ab,
                                int64_t ldab, lapack_complex_float* b,
                                int64_t ldb );
int64_t LAPACKE_ztbtrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t kd,
                                int64_t nrhs,
                                const lapack_complex_double* ab,
                                int64_t ldab, lapack_complex_double* b,
                                int64_t ldb );

int64_t LAPACKE_stfsm_work_64( int matrix_layout, char transr, char side,
                               char uplo, char trans, char diag, int64_t m,
                               int64_t n, float alpha, const float* a,
                               float* b, int64_t ldb );
int64_t LAPACKE_dtfsm_work_64( int matrix_layout, char transr, char side,
                               char uplo, char trans, char diag, int64_t m,
                               int64_t n, double alpha, const double* a,
                               double* b, int64_t ldb );
int64_t LAPACKE_ctfsm_work_64( int matrix_layout, char transr, char side,
                               char uplo, char trans, char diag, int64_t m,
                               int64_t n, lapack_complex_float alpha,
                               const lapack_complex_float* a,
                               lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_ztfsm_work_64( int matrix_layout, char transr, char side,
                               char uplo, char trans, char diag, int64_t m,
                               int64_t n, lapack_complex_double alpha,
                               const lapack_complex_double* a,
                               lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_stftri_work_64( int matrix_layout, char transr, char uplo,
                                char diag, int64_t n, float* a );
int64_t LAPACKE_dtftri_work_64( int matrix_layout, char transr, char uplo,
                                char diag, int64_t n, double* a );
int64_t LAPACKE_ctftri_work_64( int matrix_layout, char transr, char uplo,
                                char diag, int64_t n,
                                lapack_complex_float* a );
int64_t LAPACKE_ztftri_work_64( int matrix_layout, char transr, char uplo,
                                char diag, int64_t n,
                                lapack_complex_double* a );

int64_t LAPACKE_stfttp_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const float* arf, float* ap );
int64_t LAPACKE_dtfttp_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const double* arf, double* ap );
int64_t LAPACKE_ctfttp_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const lapack_complex_float* arf,
                                lapack_complex_float* ap );
int64_t LAPACKE_ztfttp_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const lapack_complex_double* arf,
                                lapack_complex_double* ap );

int64_t LAPACKE_stfttr_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const float* arf, float* a,
                                int64_t lda );
int64_t LAPACKE_dtfttr_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const double* arf, double* a,
                                int64_t lda );
int64_t LAPACKE_ctfttr_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const lapack_complex_float* arf,
                                lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_ztfttr_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const lapack_complex_double* arf,
                                lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_stgevc_work_64( int matrix_layout, char side, char howmny,
                                const lapack_logical* select, int64_t n,
                                const float* s, int64_t lds, const float* p,
                                int64_t ldp, float* vl, int64_t ldvl,
                                float* vr, int64_t ldvr, int64_t mm,
                                int64_t* m, float* work );
int64_t LAPACKE_dtgevc_work_64( int matrix_layout, char side, char howmny,
                                const lapack_logical* select, int64_t n,
                                const double* s, int64_t lds,
                                const double* p, int64_t ldp, double* vl,
                                int64_t ldvl, double* vr, int64_t ldvr,
                                int64_t mm, int64_t* m, double* work );
int64_t LAPACKE_ctgevc_work_64( int matrix_layout, char side, char howmny,
                                const lapack_logical* select, int64_t n,
                                const lapack_complex_float* s, int64_t lds,
                                const lapack_complex_float* p, int64_t ldp,
                                lapack_complex_float* vl, int64_t ldvl,
                                lapack_complex_float* vr, int64_t ldvr,
                                int64_t mm, int64_t* m,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_ztgevc_work_64( int matrix_layout, char side, char howmny,
                                const lapack_logical* select, int64_t n,
                                const lapack_complex_double* s, int64_t lds,
                                const lapack_complex_double* p, int64_t ldp,
                                lapack_complex_double* vl, int64_t ldvl,
                                lapack_complex_double* vr, int64_t ldvr,
                                int64_t mm, int64_t* m,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_stgexc_work_64( int matrix_layout, lapack_logical wantq,
                                lapack_logical wantz, int64_t n, float* a,
                                int64_t lda, float* b, int64_t ldb,
                                float* q, int64_t ldq, float* z,
                                int64_t ldz, int64_t* ifst,
                                int64_t* ilst, float* work,
                                int64_t lwork );
int64_t LAPACKE_dtgexc_work_64( int matrix_layout, lapack_logical wantq,
                                lapack_logical wantz, int64_t n, double* a,
                                int64_t lda, double* b, int64_t ldb,
                                double* q, int64_t ldq, double* z,
                                int64_t ldz, int64_t* ifst,
                                int64_t* ilst, double* work,
                                int64_t lwork );
int64_t LAPACKE_ctgexc_work_64( int matrix_layout, lapack_logical wantq,
                                lapack_logical wantz, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* q, int64_t ldq,
                                lapack_complex_float* z, int64_t ldz,
                                int64_t ifst, int64_t ilst );
int64_t LAPACKE_ztgexc_work_64( int matrix_layout, lapack_logical wantq,
                                lapack_logical wantz, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* q, int64_t ldq,
                                lapack_complex_double* z, int64_t ldz,
                                int64_t ifst, int64_t ilst );

int64_t LAPACKE_stgsen_work_64( int matrix_layout, int64_t ijob,
                                lapack_logical wantq, lapack_logical wantz,
                                const lapack_logical* select, int64_t n,
                                float* a, int64_t lda, float* b,
                                int64_t ldb, float* alphar, float* alphai,
                                float* beta, float* q, int64_t ldq, float* z,
                                int64_t ldz, int64_t* m, float* pl,
                                float* pr, float* dif, float* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_dtgsen_work_64( int matrix_layout, int64_t ijob,
                                lapack_logical wantq, lapack_logical wantz,
                                const lapack_logical* select, int64_t n,
                                double* a, int64_t lda, double* b,
                                int64_t ldb, double* alphar, double* alphai,
                                double* beta, double* q, int64_t ldq,
                                double* z, int64_t ldz, int64_t* m,
                                double* pl, double* pr, double* dif,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_ctgsen_work_64( int matrix_layout, int64_t ijob,
                                lapack_logical wantq, lapack_logical wantz,
                                const lapack_logical* select, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* alpha,
                                lapack_complex_float* beta,
                                lapack_complex_float* q, int64_t ldq,
                                lapack_complex_float* z, int64_t ldz,
                                int64_t* m, float* pl, float* pr, float* dif,
                                lapack_complex_float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_ztgsen_work_64( int matrix_layout, int64_t ijob,
                                lapack_logical wantq, lapack_logical wantz,
                                const lapack_logical* select, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* alpha,
                                lapack_complex_double* beta,
                                lapack_complex_double* q, int64_t ldq,
                                lapack_complex_double* z, int64_t ldz,
                                int64_t* m, double* pl, double* pr,
                                double* dif, lapack_complex_double* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );

int64_t LAPACKE_stgsja_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t p,
                                int64_t n, int64_t k, int64_t l,
                                float* a, int64_t lda, float* b,
                                int64_t ldb, float tola, float tolb,
                                float* alpha, float* beta, float* u,
                                int64_t ldu, float* v, int64_t ldv,
                                float* q, int64_t ldq, float* work,
                                int64_t* ncycle );
int64_t LAPACKE_dtgsja_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t p,
                                int64_t n, int64_t k, int64_t l,
                                double* a, int64_t lda, double* b,
                                int64_t ldb, double tola, double tolb,
                                double* alpha, double* beta, double* u,
                                int64_t ldu, double* v, int64_t ldv,
                                double* q, int64_t ldq, double* work,
                                int64_t* ncycle );
int64_t LAPACKE_ctgsja_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t p,
                                int64_t n, int64_t k, int64_t l,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                float tola, float tolb, float* alpha,
                                float* beta, lapack_complex_float* u,
                                int64_t ldu, lapack_complex_float* v,
                                int64_t ldv, lapack_complex_float* q,
                                int64_t ldq, lapack_complex_float* work,
                                int64_t* ncycle );
int64_t LAPACKE_ztgsja_work_64( int matrix_layout, char jobu, char jobv,
                                char jobq, int64_t m, int64_t p,
                                int64_t n, int64_t k, int64_t l,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                double tola, double tolb, double* alpha,
                                double* beta, lapack_complex_double* u,
                                int64_t ldu, lapack_complex_double* v,
                                int64_t ldv, lapack_complex_double* q,
                                int64_t ldq, lapack_complex_double* work,
                                int64_t* ncycle );

int64_t LAPACKE_stgsna_work_64( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, int64_t n,
                                const float* a, int64_t lda, const float* b,
                                int64_t ldb, const float* vl,
                                int64_t ldvl, const float* vr,
                                int64_t ldvr, float* s, float* dif,
                                int64_t mm, int64_t* m, float* work,
                                int64_t lwork, int64_t* iwork );
int64_t LAPACKE_dtgsna_work_64( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, int64_t n,
                                const double* a, int64_t lda,
                                const double* b, int64_t ldb,
                                const double* vl, int64_t ldvl,
                                const double* vr, int64_t ldvr, double* s,
                                double* dif, int64_t mm, int64_t* m,
                                double* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_ctgsna_work_64( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* b, int64_t ldb,
                                const lapack_complex_float* vl, int64_t ldvl,
                                const lapack_complex_float* vr, int64_t ldvr,
                                float* s, float* dif, int64_t mm,
                                int64_t* m, lapack_complex_float* work,
                                int64_t lwork, int64_t* iwork );
int64_t LAPACKE_ztgsna_work_64( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* b, int64_t ldb,
                                const lapack_complex_double* vl,
                                int64_t ldvl,
                                const lapack_complex_double* vr,
                                int64_t ldvr, double* s, double* dif,
                                int64_t mm, int64_t* m,
                                lapack_complex_double* work, int64_t lwork,
                                int64_t* iwork );

int64_t LAPACKE_stgsyl_work_64( int matrix_layout, char trans, int64_t ijob,
                                int64_t m, int64_t n, const float* a,
                                int64_t lda, const float* b, int64_t ldb,
                                float* c, int64_t ldc, const float* d,
                                int64_t ldd, const float* e, int64_t lde,
                                float* f, int64_t ldf, float* scale,
                                float* dif, float* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_dtgsyl_work_64( int matrix_layout, char trans, int64_t ijob,
                                int64_t m, int64_t n, const double* a,
                                int64_t lda, const double* b, int64_t ldb,
                                double* c, int64_t ldc, const double* d,
                                int64_t ldd, const double* e, int64_t lde,
                                double* f, int64_t ldf, double* scale,
                                double* dif, double* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_ctgsyl_work_64( int matrix_layout, char trans, int64_t ijob,
                                int64_t m, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* c, int64_t ldc,
                                const lapack_complex_float* d, int64_t ldd,
                                const lapack_complex_float* e, int64_t lde,
                                lapack_complex_float* f, int64_t ldf,
                                float* scale, float* dif,
                                lapack_complex_float* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_ztgsyl_work_64( int matrix_layout, char trans, int64_t ijob,
                                int64_t m, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* c, int64_t ldc,
                                const lapack_complex_double* d, int64_t ldd,
                                const lapack_complex_double* e, int64_t lde,
                                lapack_complex_double* f, int64_t ldf,
                                double* scale, double* dif,
                                lapack_complex_double* work, int64_t lwork,
                                int64_t* iwork );

int64_t LAPACKE_stpcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n, const float* ap,
                                float* rcond, float* work, int64_t* iwork );
int64_t LAPACKE_dtpcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n, const double* ap,
                                double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_ctpcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n,
                                const lapack_complex_float* ap, float* rcond,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_ztpcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n,
                                const lapack_complex_double* ap, double* rcond,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_stprfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const float* ap, const float* b, int64_t ldb,
                                const float* x, int64_t ldx, float* ferr,
                                float* berr, float* work, int64_t* iwork );
int64_t LAPACKE_dtprfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const double* ap, const double* b,
                                int64_t ldb, const double* x, int64_t ldx,
                                double* ferr, double* berr, double* work,
                                int64_t* iwork );
int64_t LAPACKE_ctprfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const lapack_complex_float* ap,
                                const lapack_complex_float* b, int64_t ldb,
                                const lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_ztprfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const lapack_complex_double* ap,
                                const lapack_complex_double* b, int64_t ldb,
                                const lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_stptri_work_64( int matrix_layout, char uplo, char diag,
                                int64_t n, float* ap );
int64_t LAPACKE_dtptri_work_64( int matrix_layout, char uplo, char diag,
                                int64_t n, double* ap );
int64_t LAPACKE_ctptri_work_64( int matrix_layout, char uplo, char diag,
                                int64_t n, lapack_complex_float* ap );
int64_t LAPACKE_ztptri_work_64( int matrix_layout, char uplo, char diag,
                                int64_t n, lapack_complex_double* ap );

int64_t LAPACKE_stptrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const float* ap, float* b, int64_t ldb );
int64_t LAPACKE_dtptrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const double* ap, double* b, int64_t ldb );
int64_t LAPACKE_ctptrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const lapack_complex_float* ap,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_ztptrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const lapack_complex_double* ap,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_stpttf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const float* ap, float* arf );
int64_t LAPACKE_dtpttf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const double* ap, double* arf );
int64_t LAPACKE_ctpttf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const lapack_complex_float* ap,
                                lapack_complex_float* arf );
int64_t LAPACKE_ztpttf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const lapack_complex_double* ap,
                                lapack_complex_double* arf );

int64_t LAPACKE_stpttr_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* ap, float* a, int64_t lda );
int64_t LAPACKE_dtpttr_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* ap, double* a, int64_t lda );
int64_t LAPACKE_ctpttr_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* ap,
                                lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_ztpttr_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* ap,
                                lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_strcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n, const float* a,
                                int64_t lda, float* rcond, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dtrcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n, const double* a,
                                int64_t lda, double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_ctrcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                float* rcond, lapack_complex_float* work,
                                float* rwork );
int64_t LAPACKE_ztrcon_work_64( int matrix_layout, char norm, char uplo,
                                char diag, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                double* rcond, lapack_complex_double* work,
                                double* rwork );

int64_t LAPACKE_strevc_work_64( int matrix_layout, char side, char howmny,
                                lapack_logical* select, int64_t n,
                                const float* t, int64_t ldt, float* vl,
                                int64_t ldvl, float* vr, int64_t ldvr,
                                int64_t mm, int64_t* m, float* work );
int64_t LAPACKE_dtrevc_work_64( int matrix_layout, char side, char howmny,
                                lapack_logical* select, int64_t n,
                                const double* t, int64_t ldt, double* vl,
                                int64_t ldvl, double* vr, int64_t ldvr,
                                int64_t mm, int64_t* m, double* work );
int64_t LAPACKE_ctrevc_work_64( int matrix_layout, char side, char howmny,
                                const lapack_logical* select, int64_t n,
                                lapack_complex_float* t, int64_t ldt,
                                lapack_complex_float* vl, int64_t ldvl,
                                lapack_complex_float* vr, int64_t ldvr,
                                int64_t mm, int64_t* m,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_ztrevc_work_64( int matrix_layout, char side, char howmny,
                                const lapack_logical* select, int64_t n,
                                lapack_complex_double* t, int64_t ldt,
                                lapack_complex_double* vl, int64_t ldvl,
                                lapack_complex_double* vr, int64_t ldvr,
                                int64_t mm, int64_t* m,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_strexc_work_64( int matrix_layout, char compq, int64_t n,
                                float* t, int64_t ldt, float* q,
                                int64_t ldq, int64_t* ifst,
                                int64_t* ilst, float* work );
int64_t LAPACKE_dtrexc_work_64( int matrix_layout, char compq, int64_t n,
                                double* t, int64_t ldt, double* q,
                                int64_t ldq, int64_t* ifst,
                                int64_t* ilst, double* work );
int64_t LAPACKE_ctrexc_work_64( int matrix_layout, char compq, int64_t n,
                                lapack_complex_float* t, int64_t ldt,
                                lapack_complex_float* q, int64_t ldq,
                                int64_t ifst, int64_t ilst );
int64_t LAPACKE_ztrexc_work_64( int matrix_layout, char compq, int64_t n,
                                lapack_complex_double* t, int64_t ldt,
                                lapack_complex_double* q, int64_t ldq,
                                int64_t ifst, int64_t ilst );

int64_t LAPACKE_strrfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const float* a, int64_t lda, const float* b,
                                int64_t ldb, const float* x, int64_t ldx,
                                float* ferr, float* berr, float* work,
                                int64_t* iwork );
int64_t LAPACKE_dtrrfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const double* a, int64_t lda,
                                const double* b, int64_t ldb,
                                const double* x, int64_t ldx, double* ferr,
                                double* berr, double* work, int64_t* iwork );
int64_t LAPACKE_ctrrfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* b, int64_t ldb,
                                const lapack_complex_float* x, int64_t ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork );
int64_t LAPACKE_ztrrfs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* b, int64_t ldb,
                                const lapack_complex_double* x, int64_t ldx,
                                double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

int64_t LAPACKE_strsen_work_64( int matrix_layout, char job, char compq,
                                const lapack_logical* select, int64_t n,
                                float* t, int64_t ldt, float* q,
                                int64_t ldq, float* wr, float* wi,
                                int64_t* m, float* s, float* sep,
                                float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dtrsen_work_64( int matrix_layout, char job, char compq,
                                const lapack_logical* select, int64_t n,
                                double* t, int64_t ldt, double* q,
                                int64_t ldq, double* wr, double* wi,
                                int64_t* m, double* s, double* sep,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_ctrsen_work_64( int matrix_layout, char job, char compq,
                                const lapack_logical* select, int64_t n,
                                lapack_complex_float* t, int64_t ldt,
                                lapack_complex_float* q, int64_t ldq,
                                lapack_complex_float* w, int64_t* m,
                                float* s, float* sep,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_ztrsen_work_64( int matrix_layout, char job, char compq,
                                const lapack_logical* select, int64_t n,
                                lapack_complex_double* t, int64_t ldt,
                                lapack_complex_double* q, int64_t ldq,
                                lapack_complex_double* w, int64_t* m,
                                double* s, double* sep,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_strsna_work_64( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, int64_t n,
                                const float* t, int64_t ldt, const float* vl,
                                int64_t ldvl, const float* vr,
                                int64_t ldvr, float* s, float* sep,
                                int64_t mm, int64_t* m, float* work,
                                int64_t ldwork, int64_t* iwork );
int64_t LAPACKE_dtrsna_work_64( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, int64_t n,
                                const double* t, int64_t ldt,
                                const double* vl, int64_t ldvl,
                                const double* vr, int64_t ldvr, double* s,
                                double* sep, int64_t mm, int64_t* m,
                                double* work, int64_t ldwork,
                                int64_t* iwork );
int64_t LAPACKE_ctrsna_work_64( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, int64_t n,
                                const lapack_complex_float* t, int64_t ldt,
                                const lapack_complex_float* vl, int64_t ldvl,
                                const lapack_complex_float* vr, int64_t ldvr,
                                float* s, float* sep, int64_t mm,
                                int64_t* m, lapack_complex_float* work,
                                int64_t ldwork, float* rwork );
int64_t LAPACKE_ztrsna_work_64( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, int64_t n,
                                const lapack_complex_double* t, int64_t ldt,
                                const lapack_complex_double* vl,
                                int64_t ldvl,
                                const lapack_complex_double* vr,
                                int64_t ldvr, double* s, double* sep,
                                int64_t mm, int64_t* m,
                                lapack_complex_double* work, int64_t ldwork,
                                double* rwork );

int64_t LAPACKE_strsyl_work_64( int matrix_layout, char trana, char tranb,
                                int64_t isgn, int64_t m, int64_t n,
                                const float* a, int64_t lda, const float* b,
                                int64_t ldb, float* c, int64_t ldc,
                                float* scale );
int64_t LAPACKE_dtrsyl_work_64( int matrix_layout, char trana, char tranb,
                                int64_t isgn, int64_t m, int64_t n,
                                const double* a, int64_t lda,
                                const double* b, int64_t ldb, double* c,
                                int64_t ldc, double* scale );
int64_t LAPACKE_ctrsyl_work_64( int matrix_layout, char trana, char tranb,
                                int64_t isgn, int64_t m, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* c, int64_t ldc,
                                float* scale );
int64_t LAPACKE_ztrsyl_work_64( int matrix_layout, char trana, char tranb,
                                int64_t isgn, int64_t m, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* c, int64_t ldc,
                                double* scale );

int64_t LAPACKE_strsyl3_work_64( int matrix_layout, char trana, char tranb,
                                 int64_t isgn, int64_t m, int64_t n,
                                 const float* a, int64_t lda,
                                 const float* b, int64_t ldb,
                                 float* c, int64_t ldc, float* scale,
                                 int64_t* iwork, int64_t liwork,
                                 float* swork, int64_t ldswork );
int64_t LAPACKE_dtrsyl3_work_64( int matrix_layout, char trana, char tranb,
                                 int64_t isgn, int64_t m, int64_t n,
                                 const double* a, int64_t lda,
                                 const double* b, int64_t ldb,
                                 double* c, int64_t ldc, double* scale,
                                 int64_t* iwork, int64_t liwork,
                                 double* swork, int64_t ldswork );
int64_t LAPACKE_ctrsyl3_work_64( int matrix_layout, char trana, char tranb,
                                 int64_t isgn, int64_t m, int64_t n,
                                 const lapack_complex_float* a, int64_t lda,
                                 const lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* c, int64_t ldc,
                                 float* scale, float* swork,
                                 int64_t ldswork );
int64_t LAPACKE_ztrsyl3_work_64( int matrix_layout, char trana, char tranb,
                                 int64_t isgn, int64_t m, int64_t n,
                                 const lapack_complex_double* a, int64_t lda,
                                 const lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* c, int64_t ldc,
                                 double* scale, double* swork,
                                 int64_t ldswork );

int64_t LAPACKE_strtri_work_64( int matrix_layout, char uplo, char diag,
                                int64_t n, float* a, int64_t lda );
int64_t LAPACKE_dtrtri_work_64( int matrix_layout, char uplo, char diag,
                                int64_t n, double* a, int64_t lda );
int64_t LAPACKE_ctrtri_work_64( int matrix_layout, char uplo, char diag,
                                int64_t n, lapack_complex_float* a,
                                int64_t lda );
int64_t LAPACKE_ztrtri_work_64( int matrix_layout, char uplo, char diag,
                                int64_t n, lapack_complex_double* a,
                                int64_t lda );

int64_t LAPACKE_strtrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const float* a, int64_t lda, float* b,
                                int64_t ldb );
int64_t LAPACKE_dtrtrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const double* a, int64_t lda, double* b,
                                int64_t ldb );
int64_t LAPACKE_ctrtrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_ztrtrs_work_64( int matrix_layout, char uplo, char trans,
                                char diag, int64_t n, int64_t nrhs,
                                const lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_strttf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const float* a, int64_t lda,
                                float* arf );
int64_t LAPACKE_dtrttf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const double* a, int64_t lda,
                                double* arf );
int64_t LAPACKE_ctrttf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* arf );
int64_t LAPACKE_ztrttf_work_64( int matrix_layout, char transr, char uplo,
                                int64_t n, const lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* arf );

int64_t LAPACKE_strttp_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* a, int64_t lda, float* ap );
int64_t LAPACKE_dtrttp_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* a, int64_t lda, double* ap );
int64_t LAPACKE_ctrttp_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* ap );
int64_t LAPACKE_ztrttp_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* ap );

int64_t LAPACKE_stzrzf_work_64( int matrix_layout, int64_t m, int64_t n,
                                float* a, int64_t lda, float* tau,
                                float* work, int64_t lwork );
int64_t LAPACKE_dtzrzf_work_64( int matrix_layout, int64_t m, int64_t n,
                                double* a, int64_t lda, double* tau,
                                double* work, int64_t lwork );
int64_t LAPACKE_ctzrzf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_ztzrzf_work_64( int matrix_layout, int64_t m, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cungbr_work_64( int matrix_layout, char vect, int64_t m,
                                int64_t n, int64_t k,
                                lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zungbr_work_64( int matrix_layout, char vect, int64_t m,
                                int64_t n, int64_t k,
                                lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunghr_work_64( int matrix_layout, int64_t n, int64_t ilo,
                                int64_t ihi, lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunghr_work_64( int matrix_layout, int64_t n, int64_t ilo,
                                int64_t ihi, lapack_complex_double* a,
                                int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunglq_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunglq_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, lapack_complex_double* a,
                                int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cungql_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zungql_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, lapack_complex_double* a,
                                int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cungqr_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zungqr_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, lapack_complex_double* a,
                                int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cungrq_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zungrq_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t k, lapack_complex_double* a,
                                int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cungtr_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* tau,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zungtr_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cungtsqr_row_work_64( int matrix_layout,
                                      int64_t m, int64_t n,
                                      int64_t mb, int64_t nb,
                                      lapack_complex_float* a, int64_t lda,
                                      const lapack_complex_float* t, int64_t ldt,
                                      lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zungtsqr_row_work_64( int matrix_layout,
                                      int64_t m, int64_t n,
                                      int64_t mb, int64_t nb,
                                      lapack_complex_double* a, int64_t lda,
                                      const lapack_complex_double* t, int64_t ldt,
                                      lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunmbr_work_64( int matrix_layout, char vect, char side,
                                char trans, int64_t m, int64_t n,
                                int64_t k, const lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunmbr_work_64( int matrix_layout, char vect, char side,
                                char trans, int64_t m, int64_t n,
                                int64_t k, const lapack_complex_double* a,
                                int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunmhr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t ilo,
                                int64_t ihi, const lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunmhr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t ilo,
                                int64_t ihi, const lapack_complex_double* a,
                                int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunmlq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunmlq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunmql_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunmql_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunmqr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunmqr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunmrq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunmrq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunmrz_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                int64_t l, const lapack_complex_float* a,
                                int64_t lda, const lapack_complex_float* tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunmrz_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                int64_t l, const lapack_complex_double* a,
                                int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cunmtr_work_64( int matrix_layout, char side, char uplo,
                                char trans, int64_t m, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zunmtr_work_64( int matrix_layout, char side, char uplo,
                                char trans, int64_t m, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_cupgtr_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* ap,
                                const lapack_complex_float* tau,
                                lapack_complex_float* q, int64_t ldq,
                                lapack_complex_float* work );
int64_t LAPACKE_zupgtr_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* ap,
                                const lapack_complex_double* tau,
                                lapack_complex_double* q, int64_t ldq,
                                lapack_complex_double* work );

int64_t LAPACKE_cupmtr_work_64( int matrix_layout, char side, char uplo,
                                char trans, int64_t m, int64_t n,
                                const lapack_complex_float* ap,
                                const lapack_complex_float* tau,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work );
int64_t LAPACKE_zupmtr_work_64( int matrix_layout, char side, char uplo,
                                char trans, int64_t m, int64_t n,
                                const lapack_complex_double* ap,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work );

int64_t LAPACKE_claghe_64( int matrix_layout, int64_t n, int64_t k,
                           const float* d, lapack_complex_float* a,
                           int64_t lda, int64_t* iseed );
int64_t LAPACKE_zlaghe_64( int matrix_layout, int64_t n, int64_t k,
                           const double* d, lapack_complex_double* a,
                           int64_t lda, int64_t* iseed );

int64_t LAPACKE_slagsy_64( int matrix_layout, int64_t n, int64_t k,
                           const float* d, float* a, int64_t lda,
                           int64_t* iseed );
int64_t LAPACKE_dlagsy_64( int matrix_layout, int64_t n, int64_t k,
                           const double* d, double* a, int64_t lda,
                           int64_t* iseed );
int64_t LAPACKE_clagsy_64( int matrix_layout, int64_t n, int64_t k,
                           const float* d, lapack_complex_float* a,
                           int64_t lda, int64_t* iseed );
int64_t LAPACKE_zlagsy_64( int matrix_layout, int64_t n, int64_t k,
                           const double* d, lapack_complex_double* a,
                           int64_t lda, int64_t* iseed );

int64_t LAPACKE_slapmr_64( int matrix_layout, lapack_logical forwrd,
                           int64_t m, int64_t n, float* x, int64_t ldx,
                           int64_t* k );
int64_t LAPACKE_dlapmr_64( int matrix_layout, lapack_logical forwrd,
                           int64_t m, int64_t n, double* x,
                           int64_t ldx, int64_t* k );
int64_t LAPACKE_clapmr_64( int matrix_layout, lapack_logical forwrd,
                           int64_t m, int64_t n, lapack_complex_float* x,
                           int64_t ldx, int64_t* k );
int64_t LAPACKE_zlapmr_64( int matrix_layout, lapack_logical forwrd,
                           int64_t m, int64_t n, lapack_complex_double* x,
                           int64_t ldx, int64_t* k );

int64_t LAPACKE_slapmt_64( int matrix_layout, lapack_logical forwrd,
                           int64_t m, int64_t n, float* x, int64_t ldx,
                           int64_t* k );
int64_t LAPACKE_dlapmt_64( int matrix_layout, lapack_logical forwrd,
                           int64_t m, int64_t n, double* x,
                           int64_t ldx, int64_t* k );
int64_t LAPACKE_clapmt_64( int matrix_layout, lapack_logical forwrd,
                           int64_t m, int64_t n, lapack_complex_float* x,
                           int64_t ldx, int64_t* k );
int64_t LAPACKE_zlapmt_64( int matrix_layout, lapack_logical forwrd,
                           int64_t m, int64_t n, lapack_complex_double* x,
                           int64_t ldx, int64_t* k );

float LAPACKE_slapy2_64( float x, float y );
double LAPACKE_dlapy2_64( double x, double y );

float LAPACKE_slapy3_64( float x, float y, float z );
double LAPACKE_dlapy3_64( double x, double y, double z );

int64_t LAPACKE_slartgp_64( float f, float g, float* cs, float* sn, float* r );
int64_t LAPACKE_dlartgp_64( double f, double g, double* cs, double* sn,
                            double* r );

int64_t LAPACKE_slartgs_64( float x, float y, float sigma, float* cs,
                            float* sn );
int64_t LAPACKE_dlartgs_64( double x, double y, double sigma, double* cs,
                            double* sn );


//LAPACK 3.3.0
int64_t LAPACKE_cbbcsd_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, char jobv2t, char trans, int64_t m,
                           int64_t p, int64_t q, float* theta, float* phi,
                           lapack_complex_float* u1, int64_t ldu1,
                           lapack_complex_float* u2, int64_t ldu2,
                           lapack_complex_float* v1t, int64_t ldv1t,
                           lapack_complex_float* v2t, int64_t ldv2t,
                           float* b11d, float* b11e, float* b12d, float* b12e,
                           float* b21d, float* b21e, float* b22d, float* b22e );
int64_t LAPACKE_cbbcsd_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, char jobv2t, char trans,
                                int64_t m, int64_t p, int64_t q,
                                float* theta, float* phi,
                                lapack_complex_float* u1, int64_t ldu1,
                                lapack_complex_float* u2, int64_t ldu2,
                                lapack_complex_float* v1t, int64_t ldv1t,
                                lapack_complex_float* v2t, int64_t ldv2t,
                                float* b11d, float* b11e, float* b12d,
                                float* b12e, float* b21d, float* b21e,
                                float* b22d, float* b22e, float* rwork,
                                int64_t lrwork );
int64_t LAPACKE_cheswapr_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_float* a, int64_t lda,
                             int64_t i1, int64_t i2 );
int64_t LAPACKE_cheswapr_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_float* a, int64_t lda,
                                  int64_t i1, int64_t i2 );
int64_t LAPACKE_chetri2_64( int matrix_layout, char uplo, int64_t n,
                            lapack_complex_float* a, int64_t lda,
                            const int64_t* ipiv );
int64_t LAPACKE_chetri2_work_64( int matrix_layout, char uplo, int64_t n,
                                 lapack_complex_float* a, int64_t lda,
                                 const int64_t* ipiv,
                                 lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_chetri2x_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_float* a, int64_t lda,
                             const int64_t* ipiv, int64_t nb );
int64_t LAPACKE_chetri2x_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_float* a, int64_t lda,
                                  const int64_t* ipiv,
                                  lapack_complex_float* work, int64_t nb );
int64_t LAPACKE_chetrs2_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_float* a,
                            int64_t lda, const int64_t* ipiv,
                            lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_chetrs2_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_float* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* work );
int64_t LAPACKE_csyconv_64( int matrix_layout, char uplo, char way, int64_t n,
                            lapack_complex_float* a, int64_t lda,
                            const int64_t* ipiv, lapack_complex_float* e  );
int64_t LAPACKE_csyconv_work_64( int matrix_layout, char uplo, char way,
                                 int64_t n, lapack_complex_float* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_float* e );
int64_t LAPACKE_csyswapr_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_float* a, int64_t lda,
                             int64_t i1, int64_t i2 );
int64_t LAPACKE_csyswapr_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_float* a, int64_t lda,
                                  int64_t i1, int64_t i2 );
int64_t LAPACKE_csytri2_64( int matrix_layout, char uplo, int64_t n,
                            lapack_complex_float* a, int64_t lda,
                            const int64_t* ipiv );
int64_t LAPACKE_csytri2_work_64( int matrix_layout, char uplo, int64_t n,
                                 lapack_complex_float* a, int64_t lda,
                                 const int64_t* ipiv,
                                 lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_csytri2x_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_float* a, int64_t lda,
                             const int64_t* ipiv, int64_t nb );
int64_t LAPACKE_csytri2x_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_float* a, int64_t lda,
                                  const int64_t* ipiv,
                                  lapack_complex_float* work, int64_t nb );
int64_t LAPACKE_csytrs2_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_float* a,
                            int64_t lda, const int64_t* ipiv,
                            lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_csytrs2_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_float* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* work );
int64_t LAPACKE_cunbdb_64( int matrix_layout, char trans, char signs,
                           int64_t m, int64_t p, int64_t q,
                           lapack_complex_float* x11, int64_t ldx11,
                           lapack_complex_float* x12, int64_t ldx12,
                           lapack_complex_float* x21, int64_t ldx21,
                           lapack_complex_float* x22, int64_t ldx22,
                           float* theta, float* phi,
                           lapack_complex_float* taup1,
                           lapack_complex_float* taup2,
                           lapack_complex_float* tauq1,
                           lapack_complex_float* tauq2 );
int64_t LAPACKE_cunbdb_work_64( int matrix_layout, char trans, char signs,
                                int64_t m, int64_t p, int64_t q,
                                lapack_complex_float* x11, int64_t ldx11,
                                lapack_complex_float* x12, int64_t ldx12,
                                lapack_complex_float* x21, int64_t ldx21,
                                lapack_complex_float* x22, int64_t ldx22,
                                float* theta, float* phi,
                                lapack_complex_float* taup1,
                                lapack_complex_float* taup2,
                                lapack_complex_float* tauq1,
                                lapack_complex_float* tauq2,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_cuncsd_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, char jobv2t, char trans, char signs,
                           int64_t m, int64_t p, int64_t q,
                           lapack_complex_float* x11, int64_t ldx11,
                           lapack_complex_float* x12, int64_t ldx12,
                           lapack_complex_float* x21, int64_t ldx21,
                           lapack_complex_float* x22, int64_t ldx22,
                           float* theta, lapack_complex_float* u1,
                           int64_t ldu1, lapack_complex_float* u2,
                           int64_t ldu2, lapack_complex_float* v1t,
                           int64_t ldv1t, lapack_complex_float* v2t,
                           int64_t ldv2t );
int64_t LAPACKE_cuncsd_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, char jobv2t, char trans,
                                char signs, int64_t m, int64_t p,
                                int64_t q, lapack_complex_float* x11,
                                int64_t ldx11, lapack_complex_float* x12,
                                int64_t ldx12, lapack_complex_float* x21,
                                int64_t ldx21, lapack_complex_float* x22,
                                int64_t ldx22, float* theta,
                                lapack_complex_float* u1, int64_t ldu1,
                                lapack_complex_float* u2, int64_t ldu2,
                                lapack_complex_float* v1t, int64_t ldv1t,
                                lapack_complex_float* v2t, int64_t ldv2t,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, int64_t lrwork,
                                int64_t* iwork );
int64_t LAPACKE_cuncsd2by1_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, int64_t m, int64_t p, int64_t q,
                           lapack_complex_float* x11, int64_t ldx11,
                           lapack_complex_float* x21, int64_t ldx21,
                           float* theta, lapack_complex_float* u1,
                           int64_t ldu1, lapack_complex_float* u2,
                           int64_t ldu2, lapack_complex_float* v1t, int64_t ldv1t );
int64_t LAPACKE_cuncsd2by1_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, int64_t m, int64_t p,
                                int64_t q, lapack_complex_float* x11, int64_t ldx11,
                                lapack_complex_float* x21, int64_t ldx21,
                                float* theta, lapack_complex_float* u1,
                                int64_t ldu1, lapack_complex_float* u2,
                                int64_t ldu2, lapack_complex_float* v1t,
                                int64_t ldv1t, lapack_complex_float* work,
                                int64_t lwork, float* rwork, int64_t lrwork,
                                int64_t* iwork );
int64_t LAPACKE_dbbcsd_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, char jobv2t, char trans, int64_t m,
                           int64_t p, int64_t q, double* theta,
                           double* phi, double* u1, int64_t ldu1, double* u2,
                           int64_t ldu2, double* v1t, int64_t ldv1t,
                           double* v2t, int64_t ldv2t, double* b11d,
                           double* b11e, double* b12d, double* b12e,
                           double* b21d, double* b21e, double* b22d,
                           double* b22e );
int64_t LAPACKE_dbbcsd_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, char jobv2t, char trans,
                                int64_t m, int64_t p, int64_t q,
                                double* theta, double* phi, double* u1,
                                int64_t ldu1, double* u2, int64_t ldu2,
                                double* v1t, int64_t ldv1t, double* v2t,
                                int64_t ldv2t, double* b11d, double* b11e,
                                double* b12d, double* b12e, double* b21d,
                                double* b21e, double* b22d, double* b22e,
                                double* work, int64_t lwork );
int64_t LAPACKE_dorbdb_64( int matrix_layout, char trans, char signs,
                           int64_t m, int64_t p, int64_t q,
                           double* x11, int64_t ldx11, double* x12,
                           int64_t ldx12, double* x21, int64_t ldx21,
                           double* x22, int64_t ldx22, double* theta,
                           double* phi, double* taup1, double* taup2,
                           double* tauq1, double* tauq2 );
int64_t LAPACKE_dorbdb_work_64( int matrix_layout, char trans, char signs,
                                int64_t m, int64_t p, int64_t q,
                                double* x11, int64_t ldx11, double* x12,
                                int64_t ldx12, double* x21, int64_t ldx21,
                                double* x22, int64_t ldx22, double* theta,
                                double* phi, double* taup1, double* taup2,
                                double* tauq1, double* tauq2, double* work,
                                int64_t lwork );
int64_t LAPACKE_dorcsd_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, char jobv2t, char trans, char signs,
                           int64_t m, int64_t p, int64_t q,
                           double* x11, int64_t ldx11, double* x12,
                           int64_t ldx12, double* x21, int64_t ldx21,
                           double* x22, int64_t ldx22, double* theta,
                           double* u1, int64_t ldu1, double* u2,
                           int64_t ldu2, double* v1t, int64_t ldv1t,
                           double* v2t, int64_t ldv2t );
int64_t LAPACKE_dorcsd_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, char jobv2t, char trans,
                                char signs, int64_t m, int64_t p,
                                int64_t q, double* x11, int64_t ldx11,
                                double* x12, int64_t ldx12, double* x21,
                                int64_t ldx21, double* x22, int64_t ldx22,
                                double* theta, double* u1, int64_t ldu1,
                                double* u2, int64_t ldu2, double* v1t,
                                int64_t ldv1t, double* v2t, int64_t ldv2t,
                                double* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_dorcsd2by1_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, int64_t m, int64_t p, int64_t q,
                           double* x11, int64_t ldx11, double* x21, int64_t ldx21,
                           double* theta, double* u1, int64_t ldu1, double* u2,
                           int64_t ldu2, double* v1t, int64_t ldv1t);
int64_t LAPACKE_dorcsd2by1_work_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, int64_t m, int64_t p, int64_t q,
                           double* x11, int64_t ldx11, double* x21, int64_t ldx21,
                           double* theta, double* u1, int64_t ldu1, double* u2,
                           int64_t ldu2, double* v1t, int64_t ldv1t,
                           double* work, int64_t lwork, int64_t* iwork );
int64_t LAPACKE_dsyconv_64( int matrix_layout, char uplo, char way, int64_t n,
                            double* a, int64_t lda, const int64_t* ipiv, double* e);
int64_t LAPACKE_dsyconv_work_64( int matrix_layout, char uplo, char way,
                                 int64_t n, double* a, int64_t lda,
                                 const int64_t* ipiv, double* e );
int64_t LAPACKE_dkyconv_64( int matrix_layout, char uplo, char way, int64_t n,
                            double* a, int64_t lda, const int64_t* ipiv, double* e);
int64_t LAPACKE_dkyconv_work_64( int matrix_layout, char uplo, char way,
                                 int64_t n, double* a, int64_t lda,
                                 const int64_t* ipiv, double* e );
int64_t LAPACKE_dsyswapr_64( int matrix_layout, char uplo, int64_t n,
                             double* a, int64_t lda, int64_t i1,
                             int64_t i2 );
int64_t LAPACKE_dsyswapr_work_64( int matrix_layout, char uplo, int64_t n,
                                  double* a, int64_t lda, int64_t i1,
                                  int64_t i2 );
int64_t LAPACKE_dkyswapr_64( int matrix_layout, char uplo, int64_t n,
                             double* a, int64_t lda, int64_t i1,
                             int64_t i2 );
int64_t LAPACKE_dkyswapr_work_64( int matrix_layout, char uplo, int64_t n,
                                  double* a, int64_t lda, int64_t i1,
                                  int64_t i2 );
int64_t LAPACKE_dsytri2_64( int matrix_layout, char uplo, int64_t n,
                            double* a, int64_t lda, const int64_t* ipiv );
int64_t LAPACKE_dsytri2_work_64( int matrix_layout, char uplo, int64_t n,
                                 double* a, int64_t lda,
                                 const int64_t* ipiv,
                                 double* work, int64_t lwork );
int64_t LAPACKE_dkytri2_64( int matrix_layout, char uplo, int64_t n,
                            double* a, int64_t lda, const int64_t* ipiv );
int64_t LAPACKE_dkytri2_work_64( int matrix_layout, char uplo, int64_t n,
                                 double* a, int64_t lda,
                                 const int64_t* ipiv,
                                 double* work, int64_t lwork );
int64_t LAPACKE_dsytri2x_64( int matrix_layout, char uplo, int64_t n,
                             double* a, int64_t lda, const int64_t* ipiv,
                             int64_t nb );
int64_t LAPACKE_dsytri2x_work_64( int matrix_layout, char uplo, int64_t n,
                                  double* a, int64_t lda,
                                  const int64_t* ipiv, double* work,
                                  int64_t nb );
int64_t LAPACKE_dkytri2x_64( int matrix_layout, char uplo, int64_t n,
                             double* a, int64_t lda, const int64_t* ipiv,
                             int64_t nb );
int64_t LAPACKE_dkytri2x_work_64( int matrix_layout, char uplo, int64_t n,
                                  double* a, int64_t lda,
                                  const int64_t* ipiv, double* work,
                                  int64_t nb );
int64_t LAPACKE_dsytrs2_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const double* a, int64_t lda,
                            const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_dsytrs2_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const double* a,
                                 int64_t lda, const int64_t* ipiv,
                                 double* b, int64_t ldb, double* work );
int64_t LAPACKE_dkytrs2_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const double* a, int64_t lda,
                            const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_dkytrs2_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const double* a,
                                 int64_t lda, const int64_t* ipiv,
                                 double* b, int64_t ldb, double* work );
int64_t LAPACKE_sbbcsd_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, char jobv2t, char trans, int64_t m,
                           int64_t p, int64_t q, float* theta, float* phi,
                           float* u1, int64_t ldu1, float* u2,
                           int64_t ldu2, float* v1t, int64_t ldv1t,
                           float* v2t, int64_t ldv2t, float* b11d,
                           float* b11e, float* b12d, float* b12e, float* b21d,
                           float* b21e, float* b22d, float* b22e );
int64_t LAPACKE_sbbcsd_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, char jobv2t, char trans,
                                int64_t m, int64_t p, int64_t q,
                                float* theta, float* phi, float* u1,
                                int64_t ldu1, float* u2, int64_t ldu2,
                                float* v1t, int64_t ldv1t, float* v2t,
                                int64_t ldv2t, float* b11d, float* b11e,
                                float* b12d, float* b12e, float* b21d,
                                float* b21e, float* b22d, float* b22e,
                                float* work, int64_t lwork );
int64_t LAPACKE_sorbdb_64( int matrix_layout, char trans, char signs,
                           int64_t m, int64_t p, int64_t q, float* x11,
                           int64_t ldx11, float* x12, int64_t ldx12,
                           float* x21, int64_t ldx21, float* x22,
                           int64_t ldx22, float* theta, float* phi,
                           float* taup1, float* taup2, float* tauq1,
                           float* tauq2 );
int64_t LAPACKE_sorbdb_work_64( int matrix_layout, char trans, char signs,
                                int64_t m, int64_t p, int64_t q,
                                float* x11, int64_t ldx11, float* x12,
                                int64_t ldx12, float* x21, int64_t ldx21,
                                float* x22, int64_t ldx22, float* theta,
                                float* phi, float* taup1, float* taup2,
                                float* tauq1, float* tauq2, float* work,
                                int64_t lwork );
int64_t LAPACKE_sorcsd_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, char jobv2t, char trans, char signs,
                           int64_t m, int64_t p, int64_t q, float* x11,
                           int64_t ldx11, float* x12, int64_t ldx12,
                           float* x21, int64_t ldx21, float* x22,
                           int64_t ldx22, float* theta, float* u1,
                           int64_t ldu1, float* u2, int64_t ldu2,
                           float* v1t, int64_t ldv1t, float* v2t,
                           int64_t ldv2t );
int64_t LAPACKE_sorcsd_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, char jobv2t, char trans,
                                char signs, int64_t m, int64_t p,
                                int64_t q, float* x11, int64_t ldx11,
                                float* x12, int64_t ldx12, float* x21,
                                int64_t ldx21, float* x22, int64_t ldx22,
                                float* theta, float* u1, int64_t ldu1,
                                float* u2, int64_t ldu2, float* v1t,
                                int64_t ldv1t, float* v2t, int64_t ldv2t,
                                float* work, int64_t lwork,
                                int64_t* iwork );
int64_t LAPACKE_sorcsd2by1_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, int64_t m, int64_t p, int64_t q,
                           float* x11, int64_t ldx11, float* x21, int64_t ldx21,
                           float* theta, float* u1, int64_t ldu1, float* u2,
                           int64_t ldu2, float* v1t, int64_t ldv1t);
int64_t LAPACKE_sorcsd2by1_work_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, int64_t m, int64_t p, int64_t q,
                           float* x11, int64_t ldx11, float* x21, int64_t ldx21,
                           float* theta, float* u1, int64_t ldu1, float* u2,
                           int64_t ldu2, float* v1t, int64_t ldv1t,
                           float* work, int64_t lwork, int64_t* iwork );
int64_t LAPACKE_ssyconv_64( int matrix_layout, char uplo, char way, int64_t n,
                            float* a, int64_t lda, const int64_t* ipiv, float* e );
int64_t LAPACKE_ssyconv_work_64( int matrix_layout, char uplo, char way,
                                 int64_t n, float* a, int64_t lda,
                                 const int64_t* ipiv, float* e );
int64_t LAPACKE_skyconv_64( int matrix_layout, char uplo, char way, int64_t n,
                            float* a, int64_t lda, const int64_t* ipiv, float* e );
int64_t LAPACKE_skyconv_work_64( int matrix_layout, char uplo, char way,
                                 int64_t n, float* a, int64_t lda,
                                 const int64_t* ipiv, float* e );
int64_t LAPACKE_ssyswapr_64( int matrix_layout, char uplo, int64_t n,
                             float* a, int64_t lda, int64_t i1,
                             int64_t i2 );
int64_t LAPACKE_ssyswapr_work_64( int matrix_layout, char uplo, int64_t n,
                                  float* a, int64_t lda, int64_t i1,
                                  int64_t i2 );
int64_t LAPACKE_skyswapr_64( int matrix_layout, char uplo, int64_t n,
                             float* a, int64_t lda, int64_t i1,
                             int64_t i2 );
int64_t LAPACKE_skyswapr_work_64( int matrix_layout, char uplo, int64_t n,
                                  float* a, int64_t lda, int64_t i1,
                                  int64_t i2 );
int64_t LAPACKE_ssytri2_64( int matrix_layout, char uplo, int64_t n, float* a,
                            int64_t lda, const int64_t* ipiv );
int64_t LAPACKE_ssytri2_work_64( int matrix_layout, char uplo, int64_t n,
                                 float* a, int64_t lda,
                                 const int64_t* ipiv,
                                 float* work, int64_t lwork );
int64_t LAPACKE_skytri2_64( int matrix_layout, char uplo, int64_t n, float* a,
                            int64_t lda, const int64_t* ipiv );
int64_t LAPACKE_skytri2_work_64( int matrix_layout, char uplo, int64_t n,
                                 float* a, int64_t lda,
                                 const int64_t* ipiv,
                                 float* work, int64_t lwork );
int64_t LAPACKE_ssytri2x_64( int matrix_layout, char uplo, int64_t n,
                             float* a, int64_t lda, const int64_t* ipiv,
                             int64_t nb );
int64_t LAPACKE_ssytri2x_work_64( int matrix_layout, char uplo, int64_t n,
                                  float* a, int64_t lda,
                                  const int64_t* ipiv, float* work,
                                  int64_t nb );
int64_t LAPACKE_skytri2x_64( int matrix_layout, char uplo, int64_t n,
                             float* a, int64_t lda, const int64_t* ipiv,
                             int64_t nb );
int64_t LAPACKE_skytri2x_work_64( int matrix_layout, char uplo, int64_t n,
                                  float* a, int64_t lda,
                                  const int64_t* ipiv, float* work,
                                  int64_t nb );
int64_t LAPACKE_ssytrs2_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const float* a, int64_t lda,
                            const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_ssytrs2_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const float* a,
                                 int64_t lda, const int64_t* ipiv,
                                 float* b, int64_t ldb, float* work );
int64_t LAPACKE_skytrs2_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const float* a, int64_t lda,
                            const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_skytrs2_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const float* a,
                                 int64_t lda, const int64_t* ipiv,
                                 float* b, int64_t ldb, float* work );
int64_t LAPACKE_zbbcsd_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, char jobv2t, char trans, int64_t m,
                           int64_t p, int64_t q, double* theta,
                           double* phi, lapack_complex_double* u1,
                           int64_t ldu1, lapack_complex_double* u2,
                           int64_t ldu2, lapack_complex_double* v1t,
                           int64_t ldv1t, lapack_complex_double* v2t,
                           int64_t ldv2t, double* b11d, double* b11e,
                           double* b12d, double* b12e, double* b21d,
                           double* b21e, double* b22d, double* b22e );
int64_t LAPACKE_zbbcsd_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, char jobv2t, char trans,
                                int64_t m, int64_t p, int64_t q,
                                double* theta, double* phi,
                                lapack_complex_double* u1, int64_t ldu1,
                                lapack_complex_double* u2, int64_t ldu2,
                                lapack_complex_double* v1t, int64_t ldv1t,
                                lapack_complex_double* v2t, int64_t ldv2t,
                                double* b11d, double* b11e, double* b12d,
                                double* b12e, double* b21d, double* b21e,
                                double* b22d, double* b22e, double* rwork,
                                int64_t lrwork );
int64_t LAPACKE_zheswapr_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_double* a, int64_t lda,
                             int64_t i1, int64_t i2 );
int64_t LAPACKE_zheswapr_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_double* a, int64_t lda,
                                  int64_t i1, int64_t i2 );
int64_t LAPACKE_zhetri2_64( int matrix_layout, char uplo, int64_t n,
                            lapack_complex_double* a, int64_t lda,
                            const int64_t* ipiv );
int64_t LAPACKE_zhetri2_work_64( int matrix_layout, char uplo, int64_t n,
                                 lapack_complex_double* a, int64_t lda,
                                 const int64_t* ipiv,
                                 lapack_complex_double* work, int64_t lwork );
int64_t LAPACKE_zhetri2x_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_double* a, int64_t lda,
                             const int64_t* ipiv, int64_t nb );
int64_t LAPACKE_zhetri2x_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_double* a, int64_t lda,
                                  const int64_t* ipiv,
                                  lapack_complex_double* work, int64_t nb );
int64_t LAPACKE_zhetrs2_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_double* a,
                            int64_t lda, const int64_t* ipiv,
                            lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zhetrs2_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_double* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* work );
int64_t LAPACKE_zsyconv_64( int matrix_layout, char uplo, char way, int64_t n,
                            lapack_complex_double* a, int64_t lda,
                            const int64_t* ipiv, lapack_complex_double* e );
int64_t LAPACKE_zsyconv_work_64( int matrix_layout, char uplo, char way,
                                 int64_t n, lapack_complex_double* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_double* e );
int64_t LAPACKE_zsyswapr_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_double* a, int64_t lda,
                             int64_t i1, int64_t i2 );
int64_t LAPACKE_zsyswapr_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_double* a, int64_t lda,
                                  int64_t i1, int64_t i2 );
int64_t LAPACKE_zsytri2_64( int matrix_layout, char uplo, int64_t n,
                            lapack_complex_double* a, int64_t lda,
                            const int64_t* ipiv );
int64_t LAPACKE_zsytri2_work_64( int matrix_layout, char uplo, int64_t n,
                                 lapack_complex_double* a, int64_t lda,
                                 const int64_t* ipiv,
                                 lapack_complex_double* work, int64_t lwork );
int64_t LAPACKE_zsytri2x_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_double* a, int64_t lda,
                             const int64_t* ipiv, int64_t nb );
int64_t LAPACKE_zsytri2x_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_double* a, int64_t lda,
                                  const int64_t* ipiv,
                                  lapack_complex_double* work, int64_t nb );
int64_t LAPACKE_zsytrs2_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_double* a,
                            int64_t lda, const int64_t* ipiv,
                            lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zsytrs2_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_double* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* work );
int64_t LAPACKE_zunbdb_64( int matrix_layout, char trans, char signs,
                           int64_t m, int64_t p, int64_t q,
                           lapack_complex_double* x11, int64_t ldx11,
                           lapack_complex_double* x12, int64_t ldx12,
                           lapack_complex_double* x21, int64_t ldx21,
                           lapack_complex_double* x22, int64_t ldx22,
                           double* theta, double* phi,
                           lapack_complex_double* taup1,
                           lapack_complex_double* taup2,
                           lapack_complex_double* tauq1,
                           lapack_complex_double* tauq2 );
int64_t LAPACKE_zunbdb_work_64( int matrix_layout, char trans, char signs,
                                int64_t m, int64_t p, int64_t q,
                                lapack_complex_double* x11, int64_t ldx11,
                                lapack_complex_double* x12, int64_t ldx12,
                                lapack_complex_double* x21, int64_t ldx21,
                                lapack_complex_double* x22, int64_t ldx22,
                                double* theta, double* phi,
                                lapack_complex_double* taup1,
                                lapack_complex_double* taup2,
                                lapack_complex_double* tauq1,
                                lapack_complex_double* tauq2,
                                lapack_complex_double* work, int64_t lwork );
int64_t LAPACKE_zuncsd_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, char jobv2t, char trans, char signs,
                           int64_t m, int64_t p, int64_t q,
                           lapack_complex_double* x11, int64_t ldx11,
                           lapack_complex_double* x12, int64_t ldx12,
                           lapack_complex_double* x21, int64_t ldx21,
                           lapack_complex_double* x22, int64_t ldx22,
                           double* theta, lapack_complex_double* u1,
                           int64_t ldu1, lapack_complex_double* u2,
                           int64_t ldu2, lapack_complex_double* v1t,
                           int64_t ldv1t, lapack_complex_double* v2t,
                           int64_t ldv2t );
int64_t LAPACKE_zuncsd_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, char jobv2t, char trans,
                                char signs, int64_t m, int64_t p,
                                int64_t q, lapack_complex_double* x11,
                                int64_t ldx11, lapack_complex_double* x12,
                                int64_t ldx12, lapack_complex_double* x21,
                                int64_t ldx21, lapack_complex_double* x22,
                                int64_t ldx22, double* theta,
                                lapack_complex_double* u1, int64_t ldu1,
                                lapack_complex_double* u2, int64_t ldu2,
                                lapack_complex_double* v1t, int64_t ldv1t,
                                lapack_complex_double* v2t, int64_t ldv2t,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, int64_t lrwork,
                                int64_t* iwork );
int64_t LAPACKE_zuncsd2by1_64( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, int64_t m, int64_t p, int64_t q,
                           lapack_complex_double* x11, int64_t ldx11,
                           lapack_complex_double* x21, int64_t ldx21,
                           double* theta, lapack_complex_double* u1,
                           int64_t ldu1, lapack_complex_double* u2,
                           int64_t ldu2, lapack_complex_double* v1t, int64_t ldv1t );
int64_t LAPACKE_zuncsd2by1_work_64( int matrix_layout, char jobu1, char jobu2,
                                char jobv1t, int64_t m, int64_t p,
                                int64_t q, lapack_complex_double* x11, int64_t ldx11,
                                lapack_complex_double* x21, int64_t ldx21,
                                double* theta, lapack_complex_double* u1,
                                int64_t ldu1, lapack_complex_double* u2,
                                int64_t ldu2, lapack_complex_double* v1t,
                                int64_t ldv1t, lapack_complex_double* work,
                                int64_t lwork, double* rwork, int64_t lrwork,
                                int64_t* iwork );

//LAPACK 3.4.0
int64_t LAPACKE_sgemqrt_64( int matrix_layout, char side, char trans,
                            int64_t m, int64_t n, int64_t k,
                            int64_t nb, const float* v, int64_t ldv,
                            const float* t, int64_t ldt, float* c,
                            int64_t ldc );
int64_t LAPACKE_dgemqrt_64( int matrix_layout, char side, char trans,
                            int64_t m, int64_t n, int64_t k,
                            int64_t nb, const double* v, int64_t ldv,
                            const double* t, int64_t ldt, double* c,
                            int64_t ldc );
int64_t LAPACKE_cgemqrt_64( int matrix_layout, char side, char trans,
                            int64_t m, int64_t n, int64_t k,
                            int64_t nb, const lapack_complex_float* v,
                            int64_t ldv, const lapack_complex_float* t,
                            int64_t ldt, lapack_complex_float* c,
                            int64_t ldc );
int64_t LAPACKE_zgemqrt_64( int matrix_layout, char side, char trans,
                            int64_t m, int64_t n, int64_t k,
                            int64_t nb, const lapack_complex_double* v,
                            int64_t ldv, const lapack_complex_double* t,
                            int64_t ldt, lapack_complex_double* c,
                            int64_t ldc );

int64_t LAPACKE_sgeqrt_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nb, float* a, int64_t lda, float* t,
                           int64_t ldt );
int64_t LAPACKE_dgeqrt_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nb, double* a, int64_t lda, double* t,
                           int64_t ldt );
int64_t LAPACKE_cgeqrt_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nb, lapack_complex_float* a,
                           int64_t lda, lapack_complex_float* t,
                           int64_t ldt );
int64_t LAPACKE_zgeqrt_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t nb, lapack_complex_double* a,
                           int64_t lda, lapack_complex_double* t,
                           int64_t ldt );

int64_t LAPACKE_sgeqrt2_64( int matrix_layout, int64_t m, int64_t n,
                            float* a, int64_t lda, float* t,
                            int64_t ldt );
int64_t LAPACKE_dgeqrt2_64( int matrix_layout, int64_t m, int64_t n,
                            double* a, int64_t lda, double* t,
                            int64_t ldt );
int64_t LAPACKE_cgeqrt2_64( int matrix_layout, int64_t m, int64_t n,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_zgeqrt2_64( int matrix_layout, int64_t m, int64_t n,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_sgeqrt3_64( int matrix_layout, int64_t m, int64_t n,
                            float* a, int64_t lda, float* t,
                            int64_t ldt );
int64_t LAPACKE_dgeqrt3_64( int matrix_layout, int64_t m, int64_t n,
                            double* a, int64_t lda, double* t,
                            int64_t ldt );
int64_t LAPACKE_cgeqrt3_64( int matrix_layout, int64_t m, int64_t n,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_zgeqrt3_64( int matrix_layout, int64_t m, int64_t n,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_stpmqrt_64( int matrix_layout, char side, char trans,
                            int64_t m, int64_t n, int64_t k,
                            int64_t l, int64_t nb, const float* v,
                            int64_t ldv, const float* t, int64_t ldt,
                            float* a, int64_t lda, float* b,
                            int64_t ldb );
int64_t LAPACKE_dtpmqrt_64( int matrix_layout, char side, char trans,
                            int64_t m, int64_t n, int64_t k,
                            int64_t l, int64_t nb, const double* v,
                            int64_t ldv, const double* t, int64_t ldt,
                            double* a, int64_t lda, double* b,
                            int64_t ldb );
int64_t LAPACKE_ctpmqrt_64( int matrix_layout, char side, char trans,
                            int64_t m, int64_t n, int64_t k,
                            int64_t l, int64_t nb,
                            const lapack_complex_float* v, int64_t ldv,
                            const lapack_complex_float* t, int64_t ldt,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_ztpmqrt_64( int matrix_layout, char side, char trans,
                            int64_t m, int64_t n, int64_t k,
                            int64_t l, int64_t nb,
                            const lapack_complex_double* v, int64_t ldv,
                            const lapack_complex_double* t, int64_t ldt,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_stpqrt_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t l, int64_t nb, float* a,
                           int64_t lda, float* b, int64_t ldb, float* t,
                           int64_t ldt );

int64_t LAPACKE_dtpqrt_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t l, int64_t nb, double* a,
                           int64_t lda, double* b, int64_t ldb, double* t,
                           int64_t ldt );
int64_t LAPACKE_ctpqrt_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t l, int64_t nb,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb,
                           lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_ztpqrt_64( int matrix_layout, int64_t m, int64_t n,
                           int64_t l, int64_t nb,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb,
                           lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_stpqrt2_64( int matrix_layout,
                            int64_t m, int64_t n, int64_t l,
                            float* a, int64_t lda,
                            float* b, int64_t ldb,
                            float* t, int64_t ldt );
int64_t LAPACKE_dtpqrt2_64( int matrix_layout,
                            int64_t m, int64_t n, int64_t l,
                            double* a, int64_t lda,
                            double* b, int64_t ldb,
                            double* t, int64_t ldt );
int64_t LAPACKE_ctpqrt2_64( int matrix_layout,
                            int64_t m, int64_t n, int64_t l,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* b, int64_t ldb,
                            lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_ztpqrt2_64( int matrix_layout,
                            int64_t m, int64_t n, int64_t l,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* b, int64_t ldb,
                            lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_stprfb_64( int matrix_layout, char side, char trans, char direct,
                           char storev, int64_t m, int64_t n,
                           int64_t k, int64_t l, const float* v,
                           int64_t ldv, const float* t, int64_t ldt,
                           float* a, int64_t lda, float* b, int64_t ldb );
int64_t LAPACKE_dtprfb_64( int matrix_layout, char side, char trans, char direct,
                           char storev, int64_t m, int64_t n,
                           int64_t k, int64_t l, const double* v,
                           int64_t ldv, const double* t, int64_t ldt,
                           double* a, int64_t lda, double* b, int64_t ldb );
int64_t LAPACKE_ctprfb_64( int matrix_layout, char side, char trans, char direct,
                           char storev, int64_t m, int64_t n,
                           int64_t k, int64_t l,
                           const lapack_complex_float* v, int64_t ldv,
                           const lapack_complex_float* t, int64_t ldt,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_ztprfb_64( int matrix_layout, char side, char trans, char direct,
                           char storev, int64_t m, int64_t n,
                           int64_t k, int64_t l,
                           const lapack_complex_double* v, int64_t ldv,
                           const lapack_complex_double* t, int64_t ldt,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sgemqrt_work_64( int matrix_layout, char side, char trans,
                                 int64_t m, int64_t n, int64_t k,
                                 int64_t nb, const float* v, int64_t ldv,
                                 const float* t, int64_t ldt, float* c,
                                 int64_t ldc, float* work );
int64_t LAPACKE_dgemqrt_work_64( int matrix_layout, char side, char trans,
                                 int64_t m, int64_t n, int64_t k,
                                 int64_t nb, const double* v, int64_t ldv,
                                 const double* t, int64_t ldt, double* c,
                                 int64_t ldc, double* work );
int64_t LAPACKE_cgemqrt_work_64( int matrix_layout, char side, char trans,
                                 int64_t m, int64_t n, int64_t k,
                                 int64_t nb, const lapack_complex_float* v,
                                 int64_t ldv, const lapack_complex_float* t,
                                 int64_t ldt, lapack_complex_float* c,
                                 int64_t ldc, lapack_complex_float* work );
int64_t LAPACKE_zgemqrt_work_64( int matrix_layout, char side, char trans,
                                 int64_t m, int64_t n, int64_t k,
                                 int64_t nb, const lapack_complex_double* v,
                                 int64_t ldv, const lapack_complex_double* t,
                                 int64_t ldt, lapack_complex_double* c,
                                 int64_t ldc, lapack_complex_double* work );

int64_t LAPACKE_sgeqrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nb, float* a, int64_t lda,
                                float* t, int64_t ldt, float* work );
int64_t LAPACKE_dgeqrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nb, double* a, int64_t lda,
                                double* t, int64_t ldt, double* work );
int64_t LAPACKE_cgeqrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nb, lapack_complex_float* a,
                                int64_t lda, lapack_complex_float* t,
                                int64_t ldt, lapack_complex_float* work );
int64_t LAPACKE_zgeqrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t nb, lapack_complex_double* a,
                                int64_t lda, lapack_complex_double* t,
                                int64_t ldt, lapack_complex_double* work );

int64_t LAPACKE_sgeqrt2_work_64( int matrix_layout, int64_t m, int64_t n,
                                 float* a, int64_t lda, float* t,
                                 int64_t ldt );
int64_t LAPACKE_dgeqrt2_work_64( int matrix_layout, int64_t m, int64_t n,
                                 double* a, int64_t lda, double* t,
                                 int64_t ldt );
int64_t LAPACKE_cgeqrt2_work_64( int matrix_layout, int64_t m, int64_t n,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_zgeqrt2_work_64( int matrix_layout, int64_t m, int64_t n,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_sgeqrt3_work_64( int matrix_layout, int64_t m, int64_t n,
                                 float* a, int64_t lda, float* t,
                                 int64_t ldt );
int64_t LAPACKE_dgeqrt3_work_64( int matrix_layout, int64_t m, int64_t n,
                                 double* a, int64_t lda, double* t,
                                 int64_t ldt );
int64_t LAPACKE_cgeqrt3_work_64( int matrix_layout, int64_t m, int64_t n,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_zgeqrt3_work_64( int matrix_layout, int64_t m, int64_t n,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_stpmqrt_work_64( int matrix_layout, char side, char trans,
                                 int64_t m, int64_t n, int64_t k,
                                 int64_t l, int64_t nb, const float* v,
                                 int64_t ldv, const float* t, int64_t ldt,
                                 float* a, int64_t lda, float* b,
                                 int64_t ldb, float* work );
int64_t LAPACKE_dtpmqrt_work_64( int matrix_layout, char side, char trans,
                                 int64_t m, int64_t n, int64_t k,
                                 int64_t l, int64_t nb, const double* v,
                                 int64_t ldv, const double* t,
                                 int64_t ldt, double* a, int64_t lda,
                                 double* b, int64_t ldb, double* work );
int64_t LAPACKE_ctpmqrt_work_64( int matrix_layout, char side, char trans,
                                 int64_t m, int64_t n, int64_t k,
                                 int64_t l, int64_t nb,
                                 const lapack_complex_float* v, int64_t ldv,
                                 const lapack_complex_float* t, int64_t ldt,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* work );
int64_t LAPACKE_ztpmqrt_work_64( int matrix_layout, char side, char trans,
                                 int64_t m, int64_t n, int64_t k,
                                 int64_t l, int64_t nb,
                                 const lapack_complex_double* v, int64_t ldv,
                                 const lapack_complex_double* t, int64_t ldt,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* work );

int64_t LAPACKE_stpqrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t l, int64_t nb, float* a,
                                int64_t lda, float* b, int64_t ldb,
                                float* t, int64_t ldt, float* work );
int64_t LAPACKE_dtpqrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t l, int64_t nb, double* a,
                                int64_t lda, double* b, int64_t ldb,
                                double* t, int64_t ldt, double* work );
int64_t LAPACKE_ctpqrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t l, int64_t nb,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* t, int64_t ldt,
                                lapack_complex_float* work );
int64_t LAPACKE_ztpqrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                int64_t l, int64_t nb,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* t, int64_t ldt,
                                lapack_complex_double* work );

int64_t LAPACKE_stpqrt2_work_64( int matrix_layout,
                                 int64_t m, int64_t n, int64_t l,
                                 float* a, int64_t lda,
                                 float* b, int64_t ldb,
                                 float* t, int64_t ldt );
int64_t LAPACKE_dtpqrt2_work_64( int matrix_layout,
                                 int64_t m, int64_t n, int64_t l,
                                 double* a, int64_t lda,
                                 double* b, int64_t ldb,
                                 double* t, int64_t ldt );
int64_t LAPACKE_ctpqrt2_work_64( int matrix_layout,
                                 int64_t m, int64_t n, int64_t l,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_ztpqrt2_work_64( int matrix_layout,
                                 int64_t m, int64_t n, int64_t l,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_stprfb_work_64( int matrix_layout, char side, char trans,
                                char direct, char storev, int64_t m,
                                int64_t n, int64_t k, int64_t l,
                                const float* v, int64_t ldv, const float* t,
                                int64_t ldt, float* a, int64_t lda,
                                float* b, int64_t ldb, float* work,
                                int64_t ldwork );
int64_t LAPACKE_dtprfb_work_64( int matrix_layout, char side, char trans,
                                char direct, char storev, int64_t m,
                                int64_t n, int64_t k, int64_t l,
                                const double* v, int64_t ldv,
                                const double* t, int64_t ldt, double* a,
                                int64_t lda, double* b, int64_t ldb,
                                double* work, int64_t ldwork );
int64_t LAPACKE_ctprfb_work_64( int matrix_layout, char side, char trans,
                                char direct, char storev, int64_t m,
                                int64_t n, int64_t k, int64_t l,
                                const lapack_complex_float* v, int64_t ldv,
                                const lapack_complex_float* t, int64_t ldt,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* b, int64_t ldb,
                                lapack_complex_float* work, int64_t ldwork );
int64_t LAPACKE_ztprfb_work_64( int matrix_layout, char side, char trans,
                                char direct, char storev, int64_t m,
                                int64_t n, int64_t k, int64_t l,
                                const lapack_complex_double* v, int64_t ldv,
                                const lapack_complex_double* t, int64_t ldt,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* b, int64_t ldb,
                                lapack_complex_double* work, int64_t ldwork );
//LAPACK 3.X.X
int64_t LAPACKE_ssysv_rook_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* a, int64_t lda,
                               int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dsysv_rook_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* a, int64_t lda,
                               int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_csysv_rook_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zsysv_rook_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_ssytrf_rook_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dsytrf_rook_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, int64_t* ipiv );
int64_t LAPACKE_csytrf_rook_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_zsytrf_rook_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv );

int64_t LAPACKE_ssytrs_rook_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_dsytrs_rook_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const double* a, int64_t lda,
                           const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_csytrs_rook_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zsytrs_rook_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_chetrf_rook_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_zhetrf_rook_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv );

int64_t LAPACKE_chetrs_rook_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_float* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zhetrs_rook_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const lapack_complex_double* a,
                           int64_t lda, const int64_t* ipiv,
                           lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_csyr_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_float alpha,
                             const lapack_complex_float* x, int64_t incx,
                             lapack_complex_float* a, int64_t lda );
int64_t LAPACKE_zsyr_64( int matrix_layout, char uplo, int64_t n,
                             lapack_complex_double alpha,
                             const lapack_complex_double* x, int64_t incx,
                             lapack_complex_double* a, int64_t lda );

int64_t LAPACKE_ssysv_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                    int64_t nrhs, float* a, int64_t lda,
                                    int64_t* ipiv, float* b, int64_t ldb,
                                    float* work, int64_t lwork );
int64_t LAPACKE_dsysv_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                    int64_t nrhs, double* a, int64_t lda,
                                    int64_t* ipiv, double* b, int64_t ldb,
                                    double* work, int64_t lwork );
int64_t LAPACKE_csysv_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                    int64_t nrhs, lapack_complex_float* a,
                                    int64_t lda, int64_t* ipiv,
                                    lapack_complex_float* b, int64_t ldb,
                                    lapack_complex_float* work,
                                    int64_t lwork );
int64_t LAPACKE_zsysv_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                    int64_t nrhs, lapack_complex_double* a,
                                    int64_t lda, int64_t* ipiv,
                                    lapack_complex_double* b, int64_t ldb,
                                    lapack_complex_double* work,
                                    int64_t lwork );

int64_t LAPACKE_ssytrf_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda, int64_t* ipiv,
                                float* work, int64_t lwork );
int64_t LAPACKE_dsytrf_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda, int64_t* ipiv,
                                double* work, int64_t lwork );
int64_t LAPACKE_csytrf_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_float* work,
                                int64_t lwork );
int64_t LAPACKE_zsytrf_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_double* work,
                                int64_t lwork );

int64_t LAPACKE_ssytrs_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const int64_t* ipiv, float* b,
                                int64_t ldb );
int64_t LAPACKE_dsytrs_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const double* a,
                                int64_t lda, const int64_t* ipiv,
                                double* b, int64_t ldb );
int64_t LAPACKE_csytrs_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zsytrs_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_chetrf_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_float* work,
                                int64_t lwork );
int64_t LAPACKE_zhetrf_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_double* work,
                                int64_t lwork );

int64_t LAPACKE_chetrs_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_float* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zhetrs_rook_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const lapack_complex_double* a,
                                int64_t lda, const int64_t* ipiv,
                                lapack_complex_double* b, int64_t ldb );


int64_t LAPACKE_csyr_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_float alpha,
                                  const lapack_complex_float* x,
                                  int64_t incx, lapack_complex_float* a,
                                  int64_t lda );
int64_t LAPACKE_zsyr_work_64( int matrix_layout, char uplo, int64_t n,
                                  lapack_complex_double alpha,
                                  const lapack_complex_double* x,
                                  int64_t incx, lapack_complex_double* a,
                                  int64_t lda );
void LAPACKE_ilaver_64( int64_t* vers_major,
                     int64_t* vers_minor,
                     int64_t* vers_patch );
// LAPACK 3.7.0
int64_t LAPACKE_ssysv_aa_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, float* a, int64_t lda,
                          int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_ssysv_aa_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* a, int64_t lda,
                               int64_t* ipiv, float* b, int64_t ldb,
                               float* work, int64_t lwork );
int64_t LAPACKE_dsysv_aa_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, double* a, int64_t lda,
                          int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_dsysv_aa_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* a, int64_t lda,
                               int64_t* ipiv, double* b, int64_t ldb,
                               double* work, int64_t lwork );
int64_t LAPACKE_csysv_aa_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, int64_t* ipiv,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_csysv_aa_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_float* b, int64_t ldb,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zsysv_aa_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, int64_t* ipiv,
                          lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zsysv_aa_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_double* b, int64_t ldb,
                               lapack_complex_double* work, int64_t lwork );
int64_t LAPACKE_chesv_aa_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, int64_t* ipiv,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_chesv_aa_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_float* b, int64_t ldb,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zhesv_aa_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, int64_t* ipiv,
                          lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zhesv_aa_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, int64_t* ipiv,
                               lapack_complex_double* b, int64_t ldb,
                               lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_ssytrf_aa_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, int64_t* ipiv );
int64_t LAPACKE_dsytrf_aa_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, int64_t* ipiv );
int64_t LAPACKE_csytrf_aa_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_zsytrf_aa_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_chetrf_aa_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           int64_t* ipiv );
int64_t LAPACKE_zhetrf_aa_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           int64_t* ipiv );

int64_t LAPACKE_ssytrf_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda, int64_t* ipiv,
                                float* work, int64_t lwork );
int64_t LAPACKE_dsytrf_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda, int64_t* ipiv,
                                double* work, int64_t lwork );
int64_t LAPACKE_csytrf_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_float* work,
                                int64_t lwork );
int64_t LAPACKE_zsytrf_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_double* work,
                                int64_t lwork );
int64_t LAPACKE_chetrf_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_float* work,
                                int64_t lwork );
int64_t LAPACKE_zhetrf_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                int64_t* ipiv, lapack_complex_double* work,
                                int64_t lwork );


int64_t LAPACKE_csytrs_aa_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_float* a,
                            int64_t lda, const int64_t* ipiv,
                            lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_csytrs_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_float* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_chetrs_aa_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_float* a,
                            int64_t lda, const int64_t* ipiv,
                            lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_chetrs_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_float* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_dsytrs_aa_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const double* a, int64_t lda,
                            const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_dsytrs_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const double* a,
                                 int64_t lda, const int64_t* ipiv,
                                 double* b, int64_t ldb, double* work, int64_t lwork );
int64_t LAPACKE_ssytrs_aa_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_ssytrs_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const int64_t* ipiv, float* b,
                                int64_t ldb, float* work, int64_t lwork );
int64_t LAPACKE_zsytrs_aa_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_double* a,
                            int64_t lda, const int64_t* ipiv,
                            lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zsytrs_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_double* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* work,  int64_t lwork);
int64_t LAPACKE_zhetrs_aa_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_double* a,
                            int64_t lda, const int64_t* ipiv,
                            lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zhetrs_aa_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_double* a,
                                 int64_t lda, const int64_t* ipiv,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* work,  int64_t lwork);


int64_t LAPACKE_ssysv_rk_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, float* a, int64_t lda,
                          float* e, int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_ssysv_rk_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* a, int64_t lda,
                               float* e, int64_t* ipiv, float* b, int64_t ldb,
                               float* work, int64_t lwork );
int64_t LAPACKE_dsysv_rk_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, double* a, int64_t lda,
                          double* e, int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_dsysv_rk_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* a, int64_t lda,
                               double* e, int64_t* ipiv, double* b, int64_t ldb,
                               double* work, int64_t lwork );
int64_t LAPACKE_csysv_rk_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* e, int64_t* ipiv,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_csysv_rk_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* e, int64_t* ipiv,
                               lapack_complex_float* b, int64_t ldb,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zsysv_rk_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* e, int64_t* ipiv,
                          lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zsysv_rk_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* e, int64_t* ipiv,
                               lapack_complex_double* b, int64_t ldb,
                               lapack_complex_double* work, int64_t lwork );
int64_t LAPACKE_chesv_rk_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* e, int64_t* ipiv,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_chesv_rk_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* e, int64_t* ipiv,
                               lapack_complex_float* b, int64_t ldb,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zhesv_rk_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* e, int64_t* ipiv,
                          lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zhesv_rk_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* e, int64_t* ipiv,
                               lapack_complex_double* b, int64_t ldb,
                               lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_ssytrf_rk_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, float* e, int64_t* ipiv );
int64_t LAPACKE_dsytrf_rk_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, double* e, int64_t* ipiv );
int64_t LAPACKE_csytrf_rk_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* e, int64_t* ipiv );
int64_t LAPACKE_zsytrf_rk_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* e, int64_t* ipiv );
int64_t LAPACKE_chetrf_rk_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           lapack_complex_float* e, int64_t* ipiv );
int64_t LAPACKE_zhetrf_rk_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           lapack_complex_double* e, int64_t* ipiv );
int64_t LAPACKE_ssytrf_rk_work_64( int matrix_layout, char uplo, int64_t n,
                                float* a, int64_t lda, float* e, int64_t* ipiv,
                                float* work, int64_t lwork );
int64_t LAPACKE_dsytrf_rk_work_64( int matrix_layout, char uplo, int64_t n,
                                double* a, int64_t lda, double* e, int64_t* ipiv,
                                double* work, int64_t lwork );
int64_t LAPACKE_csytrf_rk_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* e,
                                int64_t* ipiv, lapack_complex_float* work,
                                int64_t lwork );
int64_t LAPACKE_zsytrf_rk_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* e,
                                int64_t* ipiv, lapack_complex_double* work,
                                int64_t lwork );
int64_t LAPACKE_chetrf_rk_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                lapack_complex_float* e,
                                int64_t* ipiv, lapack_complex_float* work,
                                int64_t lwork );
int64_t LAPACKE_zhetrf_rk_work_64( int matrix_layout, char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                lapack_complex_double* e,
                                int64_t* ipiv, lapack_complex_double* work,
                                int64_t lwork );

int64_t LAPACKE_csytrs_3_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_float* a,
                            int64_t lda, const lapack_complex_float* e,
                            const int64_t* ipiv,
                            lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_csytrs_3_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_float* a,
                                 int64_t lda, const lapack_complex_float* e,
                                 const int64_t* ipiv,
                                 lapack_complex_float* b, int64_t ldb);
int64_t LAPACKE_chetrs_3_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_float* a,
                            int64_t lda, const lapack_complex_float* e,
                            const int64_t* ipiv,
                            lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_chetrs_3_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_float* a,
                                 int64_t lda, const lapack_complex_float* e,
                                 const int64_t* ipiv,
                                 lapack_complex_float* b, int64_t ldb);
int64_t LAPACKE_dsytrs_3_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const double* a, int64_t lda,
                            const double* e,
                            const int64_t* ipiv, double* b, int64_t ldb );
int64_t LAPACKE_dsytrs_3_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const double* a,
                                 int64_t lda, const double* e,
                                 const int64_t* ipiv,
                                 double* b, int64_t ldb);
int64_t LAPACKE_ssytrs_3_64( int matrix_layout, char uplo, int64_t n,
                           int64_t nrhs, const float* a, int64_t lda,
                           const float* e,
                           const int64_t* ipiv, float* b, int64_t ldb );
int64_t LAPACKE_ssytrs_3_work_64( int matrix_layout, char uplo, int64_t n,
                                int64_t nrhs, const float* a, int64_t lda,
                                const float* e, const int64_t* ipiv, float* b,
                                int64_t ldb);
int64_t LAPACKE_zsytrs_3_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_double* a,
                            int64_t lda, const lapack_complex_double* e,
                            const int64_t* ipiv,
                            lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zsytrs_3_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_double* a,
                                 int64_t lda, const lapack_complex_double* e,
                                 const int64_t* ipiv,
                                 lapack_complex_double* b, int64_t ldb);
int64_t LAPACKE_zhetrs_3_64( int matrix_layout, char uplo, int64_t n,
                            int64_t nrhs, const lapack_complex_double* a,
                            int64_t lda, const lapack_complex_double* e,
                            const int64_t* ipiv,
                            lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zhetrs_3_work_64( int matrix_layout, char uplo, int64_t n,
                                 int64_t nrhs, const lapack_complex_double* a,
                                 int64_t lda, const lapack_complex_double* e,
                                 const int64_t* ipiv,
                                 lapack_complex_double* b, int64_t ldb);

int64_t LAPACKE_ssytri_3_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, const float* e, const int64_t* ipiv );
int64_t LAPACKE_dsytri_3_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, const double* e, const int64_t* ipiv );
int64_t LAPACKE_csytri_3_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* e, const int64_t* ipiv );
int64_t LAPACKE_zsytri_3_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* e, const int64_t* ipiv );
int64_t LAPACKE_chetri_3_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* e, const int64_t* ipiv );
int64_t LAPACKE_zhetri_3_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* e, const int64_t* ipiv );
int64_t LAPACKE_ssytri_3_work_64( int matrix_layout, char uplo, int64_t n, float* a,
                           int64_t lda, const float* e, const int64_t* ipiv,
                           float* work, int64_t lwork  );
int64_t LAPACKE_dsytri_3_work_64( int matrix_layout, char uplo, int64_t n, double* a,
                           int64_t lda, const double* e, const int64_t* ipiv,
                           double* work, int64_t lwork  );
int64_t LAPACKE_csytri_3_work_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* e, const int64_t* ipiv,
                           lapack_complex_float* work, int64_t lwork  );
int64_t LAPACKE_zsytri_3_work_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* e, const int64_t* ipiv,
                           lapack_complex_double* work, int64_t lwork  );
int64_t LAPACKE_chetri_3_work_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* e, const int64_t* ipiv,
                           lapack_complex_float* work, int64_t lwork  );
int64_t LAPACKE_zhetri_3_work_64( int matrix_layout, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* e, const int64_t* ipiv,
                           lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_ssycon_3_64( int matrix_layout, char uplo, int64_t n,
                           const float* a, int64_t lda, const float* e,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_dsycon_3_64( int matrix_layout, char uplo, int64_t n,
                           const double* a, int64_t lda, const double* e,
                           const int64_t* ipiv, double anorm,
                           double* rcond );
int64_t LAPACKE_csycon_3_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* e,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_zsycon_3_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* e,
                           const int64_t* ipiv, double anorm,
                           double* rcond );
int64_t LAPACKE_checon_3_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* e,
                           const int64_t* ipiv, float anorm, float* rcond );
int64_t LAPACKE_zhecon_3_64( int matrix_layout, char uplo, int64_t n,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* e,
                           const int64_t* ipiv, double anorm,
                           double* rcond );
int64_t LAPACKE_ssycon_3_work_64( int matrix_layout, char uplo, int64_t n,
                                const float* a, int64_t lda, const float* e,
                                const int64_t* ipiv, float anorm,
                                float* rcond, float* work, int64_t* iwork );
int64_t LAPACKE_dsycon_3_work_64( int matrix_layout, char uplo, int64_t n,
                                const double* a, int64_t lda, const double* e,
                                const int64_t* ipiv, double anorm,
                                double* rcond, double* work,
                                int64_t* iwork );
int64_t LAPACKE_csycon_3_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* e,
                                const int64_t* ipiv, float anorm,
                                float* rcond, lapack_complex_float* work );
int64_t LAPACKE_zsycon_3_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* e,
                                const int64_t* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );
int64_t LAPACKE_checon_3_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* e,
                                const int64_t* ipiv, float anorm,
                                float* rcond, lapack_complex_float* work );
int64_t LAPACKE_zhecon_3_work_64( int matrix_layout, char uplo, int64_t n,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* e,
                                const int64_t* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

int64_t LAPACKE_sgelq_64( int matrix_layout, int64_t m, int64_t n,
                          float* a, int64_t lda,
                          float* t, int64_t tsize );
int64_t LAPACKE_dgelq_64( int matrix_layout, int64_t m, int64_t n,
                          double* a, int64_t lda,
                          double* t, int64_t tsize );
int64_t LAPACKE_cgelq_64( int matrix_layout, int64_t m, int64_t n,
                          lapack_complex_float* a, int64_t lda,
                          lapack_complex_float* t, int64_t tsize );
int64_t LAPACKE_zgelq_64( int matrix_layout, int64_t m, int64_t n,
                          lapack_complex_double* a, int64_t lda,
                          lapack_complex_double* t, int64_t tsize );

int64_t LAPACKE_sgelq_work_64( int matrix_layout, int64_t m, int64_t n,
                               float* a, int64_t lda,
                               float* t, int64_t tsize,
                               float* work, int64_t lwork );
int64_t LAPACKE_dgelq_work_64( int matrix_layout, int64_t m, int64_t n,
                               double* a, int64_t lda,
                               double* t, int64_t tsize,
                               double* work, int64_t lwork );
int64_t LAPACKE_cgelq_work_64( int matrix_layout, int64_t m, int64_t n,
                               lapack_complex_float* a, int64_t lda,
                               lapack_complex_float* t, int64_t tsize,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgelq_work_64( int matrix_layout, int64_t m, int64_t n,
                               lapack_complex_double* a, int64_t lda,
                               lapack_complex_double* t, int64_t tsize,
                               lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgemlq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const float* a, int64_t lda,
                           const float* t, int64_t tsize,
                           float* c, int64_t ldc );
int64_t LAPACKE_dgemlq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const double* a, int64_t lda,
                           const double* t, int64_t tsize,
                           double* c, int64_t ldc );
int64_t LAPACKE_cgemlq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* t, int64_t tsize,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zgemlq_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* t, int64_t tsize,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_sgemlq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const float* a, int64_t lda,
                                const float* t, int64_t tsize,
                                float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dgemlq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const double* a, int64_t lda,
                                const double* t, int64_t tsize,
                                double* c, int64_t ldc,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgemlq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* t, int64_t tsize,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgemlq_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* t, int64_t tsize,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgeqr_64( int matrix_layout, int64_t m, int64_t n,
                          float* a, int64_t lda,
                          float* t, int64_t tsize );
int64_t LAPACKE_dgeqr_64( int matrix_layout, int64_t m, int64_t n,
                          double* a, int64_t lda,
                          double* t, int64_t tsize );
int64_t LAPACKE_cgeqr_64( int matrix_layout, int64_t m, int64_t n,
                          lapack_complex_float* a, int64_t lda,
                          lapack_complex_float* t, int64_t tsize );
int64_t LAPACKE_zgeqr_64( int matrix_layout, int64_t m, int64_t n,
                          lapack_complex_double* a, int64_t lda,
                          lapack_complex_double* t, int64_t tsize );

int64_t LAPACKE_sgeqr_work_64( int matrix_layout, int64_t m, int64_t n,
                               float* a, int64_t lda,
                               float* t, int64_t tsize,
                               float* work, int64_t lwork );
int64_t LAPACKE_dgeqr_work_64( int matrix_layout, int64_t m, int64_t n,
                               double* a, int64_t lda,
                               double* t, int64_t tsize,
                               double* work, int64_t lwork );
int64_t LAPACKE_cgeqr_work_64( int matrix_layout, int64_t m, int64_t n,
                               lapack_complex_float* a, int64_t lda,
                               lapack_complex_float* t, int64_t tsize,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgeqr_work_64( int matrix_layout, int64_t m, int64_t n,
                               lapack_complex_double* a, int64_t lda,
                               lapack_complex_double* t, int64_t tsize,
                               lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgemqr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const float* a, int64_t lda,
                           const float* t, int64_t tsize,
                           float* c, int64_t ldc );
int64_t LAPACKE_dgemqr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const double* a, int64_t lda,
                           const double* t, int64_t tsize,
                           double* c, int64_t ldc );
int64_t LAPACKE_cgemqr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_float* a, int64_t lda,
                           const lapack_complex_float* t, int64_t tsize,
                           lapack_complex_float* c, int64_t ldc );
int64_t LAPACKE_zgemqr_64( int matrix_layout, char side, char trans,
                           int64_t m, int64_t n, int64_t k,
                           const lapack_complex_double* a, int64_t lda,
                           const lapack_complex_double* t, int64_t tsize,
                           lapack_complex_double* c, int64_t ldc );

int64_t LAPACKE_sgemqr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const float* a, int64_t lda,
                                const float* t, int64_t tsize,
                                float* c, int64_t ldc,
                                float* work, int64_t lwork );
int64_t LAPACKE_dgemqr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const double* a, int64_t lda,
                                const double* t, int64_t tsize,
                                double* c, int64_t ldc,
                                double* work, int64_t lwork );
int64_t LAPACKE_cgemqr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_float* a, int64_t lda,
                                const lapack_complex_float* t, int64_t tsize,
                                lapack_complex_float* c, int64_t ldc,
                                lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgemqr_work_64( int matrix_layout, char side, char trans,
                                int64_t m, int64_t n, int64_t k,
                                const lapack_complex_double* a, int64_t lda,
                                const lapack_complex_double* t, int64_t tsize,
                                lapack_complex_double* c, int64_t ldc,
                                lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgetsls_64( int matrix_layout, char trans, int64_t m,
                            int64_t n, int64_t nrhs, float* a,
                            int64_t lda, float* b, int64_t ldb );
int64_t LAPACKE_dgetsls_64( int matrix_layout, char trans, int64_t m,
                            int64_t n, int64_t nrhs, double* a,
                            int64_t lda, double* b, int64_t ldb );
int64_t LAPACKE_cgetsls_64( int matrix_layout, char trans, int64_t m,
                            int64_t n, int64_t nrhs,
                            lapack_complex_float* a, int64_t lda,
                            lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zgetsls_64( int matrix_layout, char trans, int64_t m,
                            int64_t n, int64_t nrhs,
                            lapack_complex_double* a, int64_t lda,
                            lapack_complex_double* b, int64_t ldb );

int64_t LAPACKE_sgetsls_work_64( int matrix_layout, char trans, int64_t m,
                                 int64_t n, int64_t nrhs, float* a,
                                 int64_t lda, float* b, int64_t ldb,
                                 float* work, int64_t lwork );
int64_t LAPACKE_dgetsls_work_64( int matrix_layout, char trans, int64_t m,
                                 int64_t n, int64_t nrhs, double* a,
                                 int64_t lda, double* b, int64_t ldb,
                                 double* work, int64_t lwork );
int64_t LAPACKE_cgetsls_work_64( int matrix_layout, char trans, int64_t m,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_float* a, int64_t lda,
                                 lapack_complex_float* b, int64_t ldb,
                                 lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgetsls_work_64( int matrix_layout, char trans, int64_t m,
                                 int64_t n, int64_t nrhs,
                                 lapack_complex_double* a, int64_t lda,
                                 lapack_complex_double* b, int64_t ldb,
                                 lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_sgetsqrhrt_64( int matrix_layout, int64_t m, int64_t n,
                               int64_t mb1, int64_t nb1, int64_t nb2,
                               float* a, int64_t lda,
                               float* t, int64_t ldt );
int64_t LAPACKE_dgetsqrhrt_64( int matrix_layout, int64_t m, int64_t n,
                               int64_t mb1, int64_t nb1, int64_t nb2,
                               double* a, int64_t lda,
                               double* t, int64_t ldt );
int64_t LAPACKE_cgetsqrhrt_64( int matrix_layout, int64_t m, int64_t n,
                               int64_t mb1, int64_t nb1, int64_t nb2,
                               lapack_complex_float* a, int64_t lda,
                               lapack_complex_float* t, int64_t ldt );
int64_t LAPACKE_zgetsqrhrt_64( int matrix_layout, int64_t m, int64_t n,
                               int64_t mb1, int64_t nb1, int64_t nb2,
                               lapack_complex_double* a, int64_t lda,
                               lapack_complex_double* t, int64_t ldt );

int64_t LAPACKE_sgetsqrhrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                    int64_t mb1, int64_t nb1, int64_t nb2,
                                    float* a, int64_t lda,
                                    float* t, int64_t ldt,
                                    float* work, int64_t lwork );
int64_t LAPACKE_dgetsqrhrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                    int64_t mb1, int64_t nb1, int64_t nb2,
                                    double* a, int64_t lda,
                                    double* t, int64_t ldt,
                                    double* work, int64_t lwork );
int64_t LAPACKE_cgetsqrhrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                    int64_t mb1, int64_t nb1, int64_t nb2,
                                    lapack_complex_float* a, int64_t lda,
                                    lapack_complex_float* t, int64_t ldt,
                                    lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zgetsqrhrt_work_64( int matrix_layout, int64_t m, int64_t n,
                                    int64_t mb1, int64_t nb1, int64_t nb2,
                                    lapack_complex_double* a, int64_t lda,
                                    lapack_complex_double* t, int64_t ldt,
                                    lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_ssyev_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          float* a, int64_t lda, float* w );
int64_t LAPACKE_dsyev_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          double* a, int64_t lda, double* w );

int64_t LAPACKE_ssyevd_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           float* a, int64_t lda, float* w );
int64_t LAPACKE_dsyevd_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           double* a, int64_t lda, double* w );

int64_t LAPACKE_ssyevr_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, float* a, int64_t lda, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* isuppz );
int64_t LAPACKE_dsyevr_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, double* a, int64_t lda, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* isuppz );

int64_t LAPACKE_ssyevx_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, float* a, int64_t lda, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_dsyevx_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, double* a, int64_t lda, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* ifail );

int64_t LAPACKE_ssyev_2stage_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, float* a, int64_t lda, float* w,
                               float* work, int64_t lwork );
int64_t LAPACKE_dsyev_2stage_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, double* a, int64_t lda,
                               double* w, double* work, int64_t lwork );

int64_t LAPACKE_ssyevd_2stage_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, float* a, int64_t lda,
                                float* w, float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dsyevd_2stage_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, double* a, int64_t lda,
                                double* w, double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_ssyevr_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, float* a,
                                int64_t lda, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w, float* z,
                                int64_t ldz, int64_t* isuppz, float* work,
                                int64_t lwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_dsyevr_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, double* a,
                                int64_t lda, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w, double* z,
                                int64_t ldz, int64_t* isuppz,
                                double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_ssyevx_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, float* a,
                                int64_t lda, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w, float* z,
                                int64_t ldz, float* work, int64_t lwork,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_dsyevx_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, double* a,
                                int64_t lda, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w, double* z,
                                int64_t ldz, double* work, int64_t lwork,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_cheev_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          lapack_complex_float* a, int64_t lda, float* w );
int64_t LAPACKE_zheev_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          lapack_complex_double* a, int64_t lda, double* w );

int64_t LAPACKE_cheevd_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           lapack_complex_float* a, int64_t lda, float* w );
int64_t LAPACKE_zheevd_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           lapack_complex_double* a, int64_t lda,
                           double* w );

int64_t LAPACKE_cheevr_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_float* a,
                           int64_t lda, float vl, float vu, int64_t il,
                           int64_t iu, float abstol, int64_t* m, float* w,
                           lapack_complex_float* z, int64_t ldz,
                           int64_t* isuppz );
int64_t LAPACKE_zheevr_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_double* a,
                           int64_t lda, double vl, double vu, int64_t il,
                           int64_t iu, double abstol, int64_t* m,
                           double* w, lapack_complex_double* z, int64_t ldz,
                           int64_t* isuppz );

int64_t LAPACKE_cheevx_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_float* a,
                           int64_t lda, float vl, float vu, int64_t il,
                           int64_t iu, float abstol, int64_t* m, float* w,
                           lapack_complex_float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_zheevx_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, lapack_complex_double* a,
                           int64_t lda, double vl, double vu, int64_t il,
                           int64_t iu, double abstol, int64_t* m,
                           double* w, lapack_complex_double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_cheev_2stage_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, lapack_complex_float* a,
                               int64_t lda, float* w,
                               lapack_complex_float* work, int64_t lwork,
                               float* rwork );
int64_t LAPACKE_zheev_2stage_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, lapack_complex_double* a,
                               int64_t lda, double* w,
                               lapack_complex_double* work, int64_t lwork,
                               double* rwork );

int64_t LAPACKE_cheevd_2stage_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, lapack_complex_float* a,
                                int64_t lda, float* w,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_zheevd_2stage_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, lapack_complex_double* a,
                                int64_t lda, double* w,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_cheevr_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                float vl, float vu, int64_t il,
                                int64_t iu, float abstol, int64_t* m,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, int64_t* isuppz,
                                lapack_complex_float* work, int64_t lwork,
                                float* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_zheevr_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                double vl, double vu, int64_t il,
                                int64_t iu, double abstol, int64_t* m,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, int64_t* isuppz,
                                lapack_complex_double* work, int64_t lwork,
                                double* rwork, int64_t lrwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_cheevx_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_float* a, int64_t lda,
                                float vl, float vu, int64_t il,
                                int64_t iu, float abstol, int64_t* m,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t* iwork, int64_t* ifail );
int64_t LAPACKE_zheevx_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n,
                                lapack_complex_double* a, int64_t lda,
                                double vl, double vu, int64_t il,
                                int64_t iu, double abstol, int64_t* m,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t* iwork, int64_t* ifail );

int64_t LAPACKE_ssbev_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t kd, float* ab, int64_t ldab, float* w,
                          float* z, int64_t ldz );
int64_t LAPACKE_dsbev_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t kd, double* ab, int64_t ldab, double* w,
                          double* z, int64_t ldz );

int64_t LAPACKE_ssbevd_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t kd, float* ab, int64_t ldab, float* w,
                           float* z, int64_t ldz );
int64_t LAPACKE_dsbevd_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t kd, double* ab, int64_t ldab,
                           double* w, double* z, int64_t ldz );

int64_t LAPACKE_ssbevx_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t kd, float* ab,
                           int64_t ldab, float* q, int64_t ldq, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, float* z, int64_t ldz,
                           int64_t* ifail );
int64_t LAPACKE_dsbevx_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t kd, double* ab,
                           int64_t ldab, double* q, int64_t ldq,
                           double vl, double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w, double* z,
                           int64_t ldz, int64_t* ifail );

int64_t LAPACKE_ssbev_2stage_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t kd, float* ab,
                               int64_t ldab, float* w, float* z,
                               int64_t ldz, float* work, int64_t lwork );
int64_t LAPACKE_dsbev_2stage_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t kd, double* ab,
                               int64_t ldab, double* w, double* z,
                               int64_t ldz, double* work, int64_t lwork );

int64_t LAPACKE_ssbevd_2stage_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t kd, float* ab,
                                int64_t ldab, float* w, float* z,
                                int64_t ldz, float* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );
int64_t LAPACKE_dsbevd_2stage_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t kd, double* ab,
                                int64_t ldab, double* w, double* z,
                                int64_t ldz, double* work, int64_t lwork,
                                int64_t* iwork, int64_t liwork );

int64_t LAPACKE_ssbevx_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t kd,
                                float* ab, int64_t ldab, float* q,
                                int64_t ldq, float vl, float vu,
                                int64_t il, int64_t iu, float abstol,
                                int64_t* m, float* w, float* z,
                                int64_t ldz, float* work, int64_t lwork, int64_t* iwork,
                                int64_t* ifail );
int64_t LAPACKE_dsbevx_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t kd,
                                double* ab, int64_t ldab, double* q,
                                int64_t ldq, double vl, double vu,
                                int64_t il, int64_t iu, double abstol,
                                int64_t* m, double* w, double* z,
                                int64_t ldz, double* work, int64_t lwork, int64_t* iwork,
                                int64_t* ifail );

int64_t LAPACKE_chbev_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t kd, lapack_complex_float* ab,
                          int64_t ldab, float* w, lapack_complex_float* z,
                          int64_t ldz );
int64_t LAPACKE_zhbev_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                          int64_t kd, lapack_complex_double* ab,
                          int64_t ldab, double* w, lapack_complex_double* z,
                          int64_t ldz );

int64_t LAPACKE_chbevd_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t kd, lapack_complex_float* ab,
                           int64_t ldab, float* w, lapack_complex_float* z,
                           int64_t ldz );
int64_t LAPACKE_zhbevd_2stage_64( int matrix_layout, char jobz, char uplo, int64_t n,
                           int64_t kd, lapack_complex_double* ab,
                           int64_t ldab, double* w, lapack_complex_double* z,
                           int64_t ldz );

int64_t LAPACKE_chbevx_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t kd,
                           lapack_complex_float* ab, int64_t ldab,
                           lapack_complex_float* q, int64_t ldq, float vl,
                           float vu, int64_t il, int64_t iu, float abstol,
                           int64_t* m, float* w, lapack_complex_float* z,
                           int64_t ldz, int64_t* ifail );
int64_t LAPACKE_zhbevx_2stage_64( int matrix_layout, char jobz, char range, char uplo,
                           int64_t n, int64_t kd,
                           lapack_complex_double* ab, int64_t ldab,
                           lapack_complex_double* q, int64_t ldq, double vl,
                           double vu, int64_t il, int64_t iu,
                           double abstol, int64_t* m, double* w,
                           lapack_complex_double* z, int64_t ldz,
                           int64_t* ifail );

int64_t LAPACKE_chbev_2stage_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t kd,
                               lapack_complex_float* ab, int64_t ldab,
                               float* w, lapack_complex_float* z,
                               int64_t ldz, lapack_complex_float* work,
                               int64_t lwork, float* rwork );
int64_t LAPACKE_zhbev_2stage_work_64( int matrix_layout, char jobz, char uplo,
                               int64_t n, int64_t kd,
                               lapack_complex_double* ab, int64_t ldab,
                               double* w, lapack_complex_double* z,
                               int64_t ldz, lapack_complex_double* work,
                               int64_t lwork, double* rwork );

int64_t LAPACKE_chbevd_2stage_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t kd,
                                lapack_complex_float* ab, int64_t ldab,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                int64_t lwork, float* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );
int64_t LAPACKE_zhbevd_2stage_work_64( int matrix_layout, char jobz, char uplo,
                                int64_t n, int64_t kd,
                                lapack_complex_double* ab, int64_t ldab,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                int64_t lwork, double* rwork,
                                int64_t lrwork, int64_t* iwork,
                                int64_t liwork );

int64_t LAPACKE_chbevx_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t kd,
                                lapack_complex_float* ab, int64_t ldab,
                                lapack_complex_float* q, int64_t ldq,
                                float vl, float vu, int64_t il,
                                int64_t iu, float abstol, int64_t* m,
                                float* w, lapack_complex_float* z,
                                int64_t ldz, lapack_complex_float* work,
                                int64_t lwork, float* rwork, int64_t* iwork,
                                int64_t* ifail );
int64_t LAPACKE_zhbevx_2stage_work_64( int matrix_layout, char jobz, char range,
                                char uplo, int64_t n, int64_t kd,
                                lapack_complex_double* ab, int64_t ldab,
                                lapack_complex_double* q, int64_t ldq,
                                double vl, double vu, int64_t il,
                                int64_t iu, double abstol, int64_t* m,
                                double* w, lapack_complex_double* z,
                                int64_t ldz, lapack_complex_double* work,
                                int64_t lwork, double* rwork, int64_t* iwork,
                                int64_t* ifail );

int64_t LAPACKE_ssygv_2stage_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, float* a, int64_t lda,
                          float* b, int64_t ldb, float* w );
int64_t LAPACKE_dsygv_2stage_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, double* a, int64_t lda,
                          double* b, int64_t ldb, double* w );
int64_t LAPACKE_ssygv_2stage_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, float* a,
                               int64_t lda, float* b, int64_t ldb,
                               float* w, float* work, int64_t lwork );
int64_t LAPACKE_dsygv_2stage_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, double* a,
                               int64_t lda, double* b, int64_t ldb,
                               double* w, double* work, int64_t lwork );

int64_t LAPACKE_chegv_2stage_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* b,
                          int64_t ldb, float* w );
int64_t LAPACKE_zhegv_2stage_64( int matrix_layout, int64_t itype, char jobz,
                          char uplo, int64_t n, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* b,
                          int64_t ldb, double* w );
int64_t LAPACKE_chegv_2stage_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* b,
                               int64_t ldb, float* w,
                               lapack_complex_float* work, int64_t lwork,
                               float* rwork );
int64_t LAPACKE_zhegv_2stage_work_64( int matrix_layout, int64_t itype, char jobz,
                               char uplo, int64_t n,
                               lapack_complex_double* a, int64_t lda,
                               lapack_complex_double* b, int64_t ldb,
                               double* w, lapack_complex_double* work,
                               int64_t lwork, double* rwork );

//LAPACK 3.8.0
int64_t LAPACKE_ssysv_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, float* a, int64_t lda,
                          float* tb, int64_t ltb, int64_t* ipiv,
                          int64_t* ipiv2, float* b, int64_t ldb );
int64_t LAPACKE_ssysv_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* a, int64_t lda,
                               float* tb, int64_t ltb, int64_t* ipiv,
                               int64_t* ipiv2, float* b, int64_t ldb,
                               float* work, int64_t lwork );
int64_t LAPACKE_dsysv_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, double* a, int64_t lda,
                          double* tb, int64_t ltb,
                          int64_t* ipiv, int64_t* ipiv2,
                          double* b, int64_t ldb );
int64_t LAPACKE_dsysv_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* a, int64_t lda,
                               double* tb, int64_t ltb,
                               int64_t* ipiv, int64_t* ipiv2,
                               double* b, int64_t ldb,
                               double* work, int64_t lwork );
int64_t LAPACKE_csysv_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_csysv_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_float* b, int64_t ldb,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zsysv_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                          lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zsysv_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_double* b, int64_t ldb,
                               lapack_complex_double* work, int64_t lwork );
int64_t LAPACKE_chesv_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_chesv_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_float* b, int64_t ldb,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zhesv_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                          lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zhesv_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_double* b, int64_t ldb,
                               lapack_complex_double* work, int64_t lwork );

int64_t LAPACKE_ssytrf_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          float* a, int64_t lda,
                          float* tb, int64_t ltb, int64_t* ipiv,
                          int64_t* ipiv2 );
int64_t LAPACKE_ssytrf_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               float* a, int64_t lda,
                               float* tb, int64_t ltb, int64_t* ipiv,
                               int64_t* ipiv2,
                               float* work, int64_t lwork );
int64_t LAPACKE_dsytrf_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          double* a, int64_t lda,
                          double* tb, int64_t ltb,
                          int64_t* ipiv, int64_t* ipiv2 );
int64_t LAPACKE_dsytrf_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               double* a, int64_t lda,
                               double* tb, int64_t ltb,
                               int64_t* ipiv, int64_t* ipiv2,
                               double* work, int64_t lwork );
int64_t LAPACKE_csytrf_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2 );
int64_t LAPACKE_csytrf_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zsytrf_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2 );
int64_t LAPACKE_zsytrf_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_double* work, int64_t lwork );
int64_t LAPACKE_chetrf_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2 );
int64_t LAPACKE_chetrf_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_float* work, int64_t lwork );
int64_t LAPACKE_zhetrf_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2 );
int64_t LAPACKE_zhetrf_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_double* work, int64_t lwork );


int64_t LAPACKE_ssytrs_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, float* a, int64_t lda,
                          float* tb, int64_t ltb, int64_t* ipiv,
                          int64_t* ipiv2, float* b, int64_t ldb );
int64_t LAPACKE_ssytrs_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, float* a, int64_t lda,
                               float* tb, int64_t ltb, int64_t* ipiv,
                               int64_t* ipiv2, float* b, int64_t ldb );
int64_t LAPACKE_dsytrs_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, double* a, int64_t lda,
                          double* tb, int64_t ltb,
                          int64_t* ipiv, int64_t* ipiv2,
                          double* b, int64_t ldb );
int64_t LAPACKE_dsytrs_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, double* a, int64_t lda,
                               double* tb, int64_t ltb,
                               int64_t* ipiv, int64_t* ipiv2,
                               double* b, int64_t ldb );
int64_t LAPACKE_csytrs_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_csytrs_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zsytrs_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                          lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zsytrs_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_chetrs_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_float* a,
                          int64_t lda, lapack_complex_float* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                          lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_chetrs_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_float* a,
                               int64_t lda, lapack_complex_float* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_float* b, int64_t ldb );
int64_t LAPACKE_zhetrs_aa_2stage_64( int matrix_layout, char uplo, int64_t n,
                          int64_t nrhs, lapack_complex_double* a,
                          int64_t lda, lapack_complex_double* tb,
                          int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                          lapack_complex_double* b, int64_t ldb );
int64_t LAPACKE_zhetrs_aa_2stage_work_64( int matrix_layout, char uplo, int64_t n,
                               int64_t nrhs, lapack_complex_double* a,
                               int64_t lda, lapack_complex_double* tb,
                               int64_t ltb, int64_t* ipiv, int64_t* ipiv2,
                               lapack_complex_double* b, int64_t ldb );
//LAPACK 3.10.0
int64_t LAPACKE_sorhr_col_64( int matrix_layout, int64_t m, int64_t n,
                              int64_t nb, float* a,
                              int64_t lda, float* t,
                              int64_t ldt, float* d );
int64_t LAPACKE_sorhr_col_work_64( int matrix_layout, int64_t m, int64_t n,
                                   int64_t nb, float* a,
                                   int64_t lda, float* t,
                                   int64_t ldt, float* d );
int64_t LAPACKE_dorhr_col_64( int matrix_layout, int64_t m, int64_t n,
                              int64_t nb, double* a,
                              int64_t lda, double* t,
                              int64_t ldt, double* d );
int64_t LAPACKE_dorhr_col_work_64( int matrix_layout, int64_t m, int64_t n,
                                   int64_t nb, double* a,
                                   int64_t lda, double* t,
                                   int64_t ldt, double* d );
int64_t LAPACKE_cunhr_col_64( int matrix_layout, int64_t m, int64_t n,
                              int64_t nb, lapack_complex_float* a,
                              int64_t lda, lapack_complex_float* t,
                              int64_t ldt, lapack_complex_float* d );
int64_t LAPACKE_cunhr_col_work_64( int matrix_layout, int64_t m, int64_t n,
                                   int64_t nb, lapack_complex_float* a,
                                   int64_t lda, lapack_complex_float* t,
                                   int64_t ldt, lapack_complex_float* d );
int64_t LAPACKE_zunhr_col_64( int matrix_layout, int64_t m, int64_t n,
                              int64_t nb, lapack_complex_double* a,
                              int64_t lda, lapack_complex_double* t,
                              int64_t ldt, lapack_complex_double* d );
int64_t LAPACKE_zunhr_col_work_64( int matrix_layout, int64_t m, int64_t n,
                                   int64_t nb, lapack_complex_double* a,
                                   int64_t lda, lapack_complex_double* t,
                                   int64_t ldt, lapack_complex_double* d );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LAPACKE_64_H_ */
