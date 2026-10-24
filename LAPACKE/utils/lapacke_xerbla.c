/*****************************************************************************
  Copyright (c) 2010, Intel Corp.
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
* Contents: Native C interface to LAPACK xerbla
* Author: Intel Corporation
*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "lapacke_utils.h"

static lapack_xerbla_callback_t active_callback = NULL;

void API_SUFFIX(LAPACKE_xerbla) LAPACKE_XERBLA_ARGS
{
    if (active_callback) {
        active_callback(name, info);
        return;
    }
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        printf( "Not enough memory to allocate work array in %s\n", name );
    } else if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
        printf( "Not enough memory to transpose matrix in %s\n", name );
    } else if( info < 0 ) {
        printf( "Wrong parameter %d in %s\n", -(int) info, name );
    }
}

/* from CBLAS */
static inline size_t LAPACKE_xerbla_trimmed_length(const char *name,
                                                   const size_t name_len)
{
    if (name == NULL) return 0;
    size_t actual_len = 0;
    while (actual_len < name_len && name[actual_len] != '\0') {
       actual_len++;
    }
    while (actual_len > 0 && name[actual_len - 1] == ' ') {
       actual_len--;
    }
    return actual_len;
}

void API_SUFFIX(LAPACKE_F77_xerbla) LAPACK_FORTRAN_XERBLA_ARGS
{
#ifdef LAPACK_FORTRAN_STRLEN_END
    size_t srname_len = len > 0 ? (size_t)len : 0;
#else
    size_t srname_len = 6;
#endif
    srname_len = MIN(31, LAPACKE_xerbla_trimmed_length(F77_srname,srname_len));
    char srname[32]; /* would use VLA but MSVC no like */
    strncpy(srname, F77_srname, srname_len);
    srname[srname_len] = '\0';
    API_SUFFIX(LAPACKE_xerbla)(srname, *info);
}

lapack_xerbla_callback_t API_SUFFIX(LAPACKE_get_xerbla)(void)
{
    return active_callback;
}

lapack_xerbla_callback_t API_SUFFIX(LAPACKE_set_xerbla)(lapack_xerbla_callback_t cb, lapack_fortran_xerbla_callback_t F77_cb)
{
    lapack_xerbla_callback_t old_cb = active_callback;
    active_callback = cb;
    LAPACK_set_lapack_xerbla(
        cb ? (F77_cb ? F77_cb : API_SUFFIX(LAPACKE_F77_xerbla)) /* set up new handler */
        : NULL /* restore default */
    );
    return old_cb;
}
