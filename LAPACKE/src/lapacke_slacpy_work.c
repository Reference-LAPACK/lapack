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
* Contents: Native middle-level C interface to LAPACK function slacpy
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_slacpy_work)( int matrix_layout, char uplo, lapack_int m,
                                lapack_int n, const float* a, lapack_int lda,
                                float* b, lapack_int ldb )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_slacpy( &uplo, &m, &n, a, &lda, b, &ldb );
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        char uplo_t = uplo;
        /* Check leading dimension(s) */
        if( lda < n ) {
            info = -6;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_slacpy_work", info );
            return info;
        }
        if( ldb < n ) {
            info = -8;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_slacpy_work", info );
            return info;
        }
        /* Fix (issue #729): the previous code transposed the FULL m-by-n
         * matrix through temporary buffers, which (a) read the uninitialized
         * complementary region of the destination temporary and (b) wrote it
         * back over the part of B that `uplo` requires to stay untouched,
         * corrupting the caller's data and reading uninitialized memory. A
         * row-major m-by-n matrix with leading dimension lda is bit-identical
         * to a column-major n-by-m matrix (its transpose). Copying triangle
         * `uplo` of the row-major matrix therefore equals copying the OPPOSITE
         * triangle of that transpose, so swap U<->L and the m/n extents and
         * call the Fortran kernel directly on a and b. This copies exactly the
         * requested triangle and never touches the complementary region -- no
         * full-matrix transpose, no temporary allocation, no uninitialized
         * reads. uplo values other than 'U'/'L' (e.g. 'A'/'a', a full copy)
         * are layout-agnostic and pass through unchanged. */
        if( uplo == 'U' || uplo == 'u' ) {
            uplo_t = 'L';
        } else if( uplo == 'L' || uplo == 'l' ) {
            uplo_t = 'U';
        }
        LAPACK_slacpy( &uplo_t, &n, &m, a, &lda, b, &ldb );
        info = 0;  /* LAPACK call is ok! */
    } else {
        info = -1;
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_slacpy_work", info );
    }
    return info;
}
