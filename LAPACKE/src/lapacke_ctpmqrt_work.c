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
* Contents: Native middle-level C interface to LAPACK function ctpmqrt
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_ctpmqrt_work)( int matrix_layout, char side, char trans,
                                 lapack_int m, lapack_int n, lapack_int k,
                                 lapack_int l, lapack_int nb,
                                 const lapack_complex_float* v, lapack_int ldv,
                                 const lapack_complex_float* t, lapack_int ldt,
                                 lapack_complex_float* a, lapack_int lda,
                                 lapack_complex_float* b, lapack_int ldb,
                                 lapack_complex_float* work )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_ctpmqrt( &side, &trans, &m, &n, &k, &l, &nb, v, &ldv, t, &ldt, a,
                        &lda, b, &ldb, work, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int nrowsA, ncolsA, nrowsV;
        if      ( API_SUFFIX(LAPACKE_lsame)(side, 'l') ) { nrowsA = k; ncolsA = n; nrowsV = m; }
        else if ( API_SUFFIX(LAPACKE_lsame)(side, 'r') ) { nrowsA = m; ncolsA = k; nrowsV = n; }
        else {
            info = -2;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_ctpmqrt_work", info );
            return info;
        }
        lapack_int lda_t = MAX(1,nrowsA);
        lapack_int ldb_t = MAX(1,m);
        lapack_int ldt_t = MAX(1,nb);
        lapack_int ldv_t = MAX(1,nrowsV);
        lapack_complex_float* v_t = NULL;
        lapack_complex_float* t_t = NULL;
        lapack_complex_float* a_t = NULL;
        lapack_complex_float* b_t = NULL;
        /* Check leading dimension(s) */
        if( lda < ncolsA ) {
            info = -14;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_ctpmqrt_work", info );
            return info;
        }
        if( ldb < n ) {
            info = -16;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_ctpmqrt_work", info );
            return info;
        }
        if( ldt < k ) {
            info = -12;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_ctpmqrt_work", info );
            return info;
        }
        if( ldv < k ) {
            info = -10;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_ctpmqrt_work", info );
            return info;
        }
        /* Allocate memory for temporary array(s) */
        v_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) * ldv_t * MAX(1,k) );
        if( v_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        t_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) * ldt_t * MAX(1,k) );
        if( t_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_1;
        }
        a_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) * lda_t * MAX(1,ncolsA) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_2;
        }
        b_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) * ldb_t * MAX(1,n) );
        if( b_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_3;
        }
        /* Transpose input matrices */
        API_SUFFIX(LAPACKE_cge_trans)( LAPACK_ROW_MAJOR, nrowsV, k, v, ldv, v_t, ldv_t );
        API_SUFFIX(LAPACKE_cge_trans)( LAPACK_ROW_MAJOR, nb, k, t, ldt, t_t, ldt_t );
        API_SUFFIX(LAPACKE_cge_trans)( LAPACK_ROW_MAJOR, nrowsA, ncolsA, a, lda, a_t, lda_t );
        API_SUFFIX(LAPACKE_cge_trans)( LAPACK_ROW_MAJOR, m, n, b, ldb, b_t, ldb_t );
        /* Call LAPACK function and adjust info */
        LAPACK_ctpmqrt( &side, &trans, &m, &n, &k, &l, &nb, v_t, &ldv_t, t_t,
                        &ldt_t, a_t, &lda_t, b_t, &ldb_t, work, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        API_SUFFIX(LAPACKE_cge_trans)( LAPACK_COL_MAJOR, nrowsA, ncolsA, a_t, lda_t, a, lda );
        API_SUFFIX(LAPACKE_cge_trans)( LAPACK_COL_MAJOR, m, n, b_t, ldb_t, b, ldb );
        /* Release memory and exit */
        LAPACKE_free( b_t );
exit_level_3:
        LAPACKE_free( a_t );
exit_level_2:
        LAPACKE_free( t_t );
exit_level_1:
        LAPACKE_free( v_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_ctpmqrt_work", info );
        }
    } else {
        info = -1;
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_ctpmqrt_work", info );
    }
    return info;
}
