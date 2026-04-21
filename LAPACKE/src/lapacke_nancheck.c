/*****************************************************************************
  Copyright (c) 2017, Intel Corp.
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
* Contents: Native C interface to control NaN checking
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

#include <stdlib.h>

static int nancheck_flag = -1;

void LAPACKE_set_nancheck( int flag )
{
    nancheck_flag = ( flag ) ? 1 : 0;
}

typedef struct {
  int found;
  char *value;
} _lapacke_env_var;

static inline const _lapacke_env_var LAPACKE_getenv(const char *var_name)
{
    size_t var_length = 0;

    /* Get the length of the environment variable value */
#if defined(_WIN32)
    errno_t result = getenv_s( &var_length, NULL, 0, var_name );
    if ( result != 0 || var_length == 0 ) {
        return (_lapacke_env_var){ 0, NULL };
    }
#else
    const char *env = getenv( var_name );
    if ( env == NULL ) {
        return (_lapacke_env_var){ 0, NULL };
    }
    var_length = strlen( env ) + 1;
#endif

    /* Allocate memory for the environment variable value */
    char *value = (char *)LAPACKE_malloc( var_length );
    if ( value == NULL ) {
        return (_lapacke_env_var){ 0, NULL };
    }

    /* Get the value of the environment variable */
#if defined(_WIN32)
    result = getenv_s( &var_length, value, var_length, var_name );
    if ( result != 0 ) {
        LAPACKE_free( value );
        return (_lapacke_env_var){ 0, NULL };
    }
#else
    memcpy( value, env, var_length );
#endif

    return (_lapacke_env_var){ 1, value };
}

static inline void LAPACKE_freenv(_lapacke_env_var *var)
{
    if ( var->found && var->value != NULL ) {
        LAPACKE_free( var->value );
    }
    var->found = 0;
    var->value = NULL;
}

int LAPACKE_get_nancheck( )
{
    if ( nancheck_flag != -1 ) {
        return nancheck_flag;
    }

    /* Check environment variable, once and only once */
    _lapacke_env_var lapacke_nancheck = LAPACKE_getenv( "LAPACKE_NANCHECK" );
    if ( !lapacke_nancheck.found ) {
        /* By default, NaN checking is enabled */
        nancheck_flag = 1;
    } else {
        nancheck_flag = atoi( lapacke_nancheck.value ) ? 1 : 0;
    }

    LAPACKE_freenv( &lapacke_nancheck );
    return nancheck_flag;
}
