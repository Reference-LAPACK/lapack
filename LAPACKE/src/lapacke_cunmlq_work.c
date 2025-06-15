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
*****************************************************************************
* Contents: Native middle-level C interface to LAPACK function cunmlq
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_cunmlq_work)( int matrix_layout, char side, char trans,
                                lapack_int m, lapack_int n, lapack_int k,
                                const lapack_complex_float* a, lapack_int lda,
                                const lapack_complex_float* tau,
                                lapack_complex_float* c, lapack_int ldc,
                                lapack_complex_float* work, lapack_int lwork )
{
    lapack_int info = 0;
    lapack_int r;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_cunmlq( &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work,
                       &lwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        r = API_SUFFIX(LAPACKE_lsame)( side, 'l' ) ? m : n;
        lapack_int lda_t = MAX(1,k);
        lapack_int ldc_t = MAX(1,m);
        lapack_complex_float* a_t = NULL;
        lapack_complex_float* c_t = NULL;
        /* Check leading dimension(s) */
        if( lda < r ) {
            info = -8;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_cunmlq_work", info );
            return info;
        }
        if( ldc < n ) {
            info = -11;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_cunmlq_work", info );
            return info;
        }
        /* Query optimal working array(s) size if requested */
        if( lwork == -1 ) {
            LAPACK_cunmlq( &side, &trans, &m, &n, &k, a, &lda_t, tau, c, &ldc_t,
                           work, &lwork, &info );
            return (info < 0) ? (info - 1) : info;
        }
        /* Allocate memory for temporary array(s) */
        if( API_SUFFIX(LAPACKE_lsame)( side, 'l' ) ) {
            a_t = (lapack_complex_float*)
                LAPACKE_malloc( sizeof(lapack_complex_float) * lda_t * MAX(1,m) );
        } else {
            a_t = (lapack_complex_float*)
                LAPACKE_malloc( sizeof(lapack_complex_float) * lda_t * MAX(1,n) );
        }
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        c_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) * ldc_t * MAX(1,n) );
        if( c_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_1;
        }
        /* Transpose input matrices */
        API_SUFFIX(LAPACKE_cge_trans)( matrix_layout, k, r, a, lda, a_t, lda_t );
        API_SUFFIX(LAPACKE_cge_trans)( matrix_layout, m, n, c, ldc, c_t, ldc_t );
        /* Call LAPACK function and adjust info */
        LAPACK_cunmlq( &side, &trans, &m, &n, &k, a_t, &lda_t, tau, c_t, &ldc_t,
                       work, &lwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        API_SUFFIX(LAPACKE_cge_trans)( LAPACK_COL_MAJOR, m, n, c_t, ldc_t, c, ldc );
        /* Release memory and exit */
        LAPACKE_free( c_t );
exit_level_1:
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_cunmlq_work", info );
        }
    } else {
        info = -1;
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_cunmlq_work", info );
    }
    return info;
}
