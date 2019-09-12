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
* Contents: Native middle-level C interface to LAPACK function dorhr
* Author: Julien Langou
* Generated July 2019
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_dorhr_work( int matrix_layout, 
	lapack_int m, lapack_int n, lapack_int mb1, lapack_int nb1, 
	double* a, lapack_int lda, double* t1, lapack_int ldt1, lapack_int nb2, 
	double* t2, lapack_int ldt2, double* d, double* work, lapack_int* lwork )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
	LAPACK_dorhr( &m, &n, &mb1, &nb1, a, &lda, t1, &ldt1, &nb2, 
		t2, &ldt2, d, work, &lwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int lda_t = MAX(1,n);
        lapack_int ldt1_t = MAX(1,MIN(nb1,n));
        lapack_int ldt2_t = MAX(1,MIN(nb2,n));
        lapack_int n_t1 = n * ceil((m-n)/(mb1-n));
        double* a_t = NULL;
        double* t1_t = NULL;
        double* t2_t = NULL;
        /* Check leading dimension(s) */
        if( lda < n ) {
            info = -7;
            LAPACKE_xerbla( "LAPACKE_dorhr_work", info );
            return info;
        }
        if( ldt1 < MIN(nb1,n) ) {
            info = -9;
            LAPACKE_xerbla( "LAPACKE_dorhr_work", info );
            return info;
        }
        if( ldt2 < MIN(nb2,n) ) {
            info = -12;
            LAPACKE_xerbla( "LAPACKE_dorhr_work", info );
            return info;
        }
        /* Allocate memory for temporary array(s) */
        a_t = (double*)LAPACKE_malloc( sizeof(double) * lda_t * MAX(1,n) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        t1_t = (double*)LAPACKE_malloc( sizeof(double) * ldt1_t * MAX(1,n*n_t1) );
        if( t1_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        t2_t = (double*)LAPACKE_malloc( sizeof(double) * ldt2_t * MAX(1,n) );
        if( t2_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        /* Transpose input matrices */
        LAPACKE_dge_trans( matrix_layout, m, n, a, lda, a_t, lda_t );
        LAPACKE_dge_trans( matrix_layout, n, n_t, t1, ldt1, t1_t, ldt1_t );
        /* Call LAPACK function and adjust info */
	LAPACK_dorhr( &m, &n, &mb1, &nb1, a_t, &lda_t, t1_t, &ldt1_t, &nb2, 
		t2_t, &ldt2_t, d, work, &lwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_dge_trans( matrix_layout, m, n, a_t, lda_t, a, lda );
        LAPACKE_dge_trans( matrix_layout, n, n_t1, t1_t, ldt1_t, t1, ldt1 );
        LAPACKE_dge_trans( matrix_layout, n, n, t2_t, ldt2_t, t2, ldt2 );
        /* Release memory and exit */
        LAPACKE_free( t2_t );
        LAPACKE_free( t1_t );
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_dorhr_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_dorhr_work", info );
    }
    return info;
}
