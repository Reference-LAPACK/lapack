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
* Contents: Native C interface to LAPACK
* Author: Intel Corporation
* Generated May, 2011
*****************************************************************************/

#ifndef _LAPACKE_CONFIG_H_
#define _LAPACKE_CONFIG_H_

#ifdef __cplusplus
#if defined(LAPACK_COMPLEX_CPP)
#include <complex>
#endif
extern "C" {
#endif /* __cplusplus */

#include <stdlib.h>

#ifndef lapack_int
#if defined(LAPACK_ILP64)
#define lapack_int              long
#else
#define lapack_int              int
#endif
#endif

#ifndef lapack_logical
#define lapack_logical          lapack_int
#endif

#ifndef LAPACK_COMPLEX_CUSTOM

#if defined(LAPACK_COMPLEX_STRUCTURE)

typedef struct { float real, imag; } _lapack_complex_float;
typedef struct { double real, imag; } _lapack_complex_double;
#define lapack_complex_float  _lapack_complex_float
#define lapack_complex_double _lapack_complex_double
#define lapack_complex_float_real(z)  ((z).real)
#define lapack_complex_float_imag(z)  ((z).imag)
#define lapack_complex_double_real(z)  ((z).real)
#define lapack_complex_double_imag(z)  ((z).imag)

#elif defined(LAPACK_COMPLEX_C99)

#include <complex.h>
#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex
#define lapack_complex_float_real(z)       (creal(z))
#define lapack_complex_float_imag(z)       (cimag(z))
#define lapack_complex_double_real(z)       (creal(z))
#define lapack_complex_double_imag(z)       (cimag(z))

#elif defined(LAPACK_COMPLEX_CPP)

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#define lapack_complex_float_real(z)       ((z).real())
#define lapack_complex_float_imag(z)       ((z).imag())
#define lapack_complex_double_real(z)       ((z).real())
#define lapack_complex_double_imag(z)       ((z).imag())

#else

#include <complex.h>
#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex
#define lapack_complex_float_real(z)       (creal(z))
#define lapack_complex_float_imag(z)       (cimag(z))
#define lapack_complex_double_real(z)       (creal(z))
#define lapack_complex_double_imag(z)       (cimag(z))

#endif

lapack_complex_float lapack_make_complex_float( float re, float im );
lapack_complex_double lapack_make_complex_double( double re, double im );

#endif

#ifndef LAPACK_malloc
#define LAPACK_malloc( size )   malloc( size )
#endif

#ifndef LAPACK_free
#define LAPACK_free( p )        free( p )
#endif

#ifndef LAPACKE_NAME
#if defined(LAPACK_NAME_PATTERN_LC)
#define LAPACKE_NAME(lcname,UCNAME)  lapacke_##lcname
#elif defined(LAPACK_NAME_PATTERN_UC)
#define LAPACKE_NAME(lcname,UCNAME)  LAPACKE_##UCNAME
#elif defined(LAPACK_NAME_PATTERN_MC)
#define LAPACKE_NAME(lcname,UCNAME)  LAPACKE_##lcname
#elif defined(LAPACK_NAME_PATTERN_LC_SHORT)
#define LAPACKE_NAME(lcname,UCNAME)  c##lcname
#elif defined(LAPACKE_NAME_PATTERN_UC_SHORT)
#define LAPACKE_NAME(lcname,UCNAME)  C##UCNAME
#elif defined(LAPACK_NAME_PATTERN_MC_SHORT)
#define LAPACKE_NAME(lcname,UCNAME)  C##lcname
#else
#define LAPACK_NAME_PATTERN_MC
#endif
#endif

#ifndef LAPACK_NAME_PATTERN_MC

#define LAPACKE_lsame    LAPACKE_NAME(lsame,LSAME)
#define LAPACKE_xerbla   LAPACKE_NAME(xerbla,XERBLA)

#define LAPACKE_sbdsdc   LAPACKE_NAME(sbdsdc,SBDSDC)
#define LAPACKE_dbdsdc   LAPACKE_NAME(dbdsdc,DBDSDC)

#define LAPACKE_sbdsqr   LAPACKE_NAME(sbdsqr,SBDSQR)
#define LAPACKE_dbdsqr   LAPACKE_NAME(dbdsqr,DBDSQR)
#define LAPACKE_cbdsqr   LAPACKE_NAME(cbdsqr,CBDSQR)
#define LAPACKE_zbdsqr   LAPACKE_NAME(zbdsqr,ZBDSQR)

#define LAPACKE_sdisna   LAPACKE_NAME(sdisna,SDISNA)
#define LAPACKE_ddisna   LAPACKE_NAME(ddisna,DDISNA)

#define LAPACKE_sgbbrd   LAPACKE_NAME(sgbbrd,SGBBRD)
#define LAPACKE_dgbbrd   LAPACKE_NAME(dgbbrd,DGBBRD)
#define LAPACKE_cgbbrd   LAPACKE_NAME(cgbbrd,CGBBRD)
#define LAPACKE_zgbbrd   LAPACKE_NAME(zgbbrd,ZGBBRD)

#define LAPACKE_sgbcon   LAPACKE_NAME(sgbcon,SGBCON)
#define LAPACKE_dgbcon   LAPACKE_NAME(dgbcon,DGBCON)
#define LAPACKE_cgbcon   LAPACKE_NAME(cgbcon,CGBCON)
#define LAPACKE_zgbcon   LAPACKE_NAME(zgbcon,ZGBCON)

#define LAPACKE_sgbequ   LAPACKE_NAME(sgbequ,SGBEQU)
#define LAPACKE_dgbequ   LAPACKE_NAME(dgbequ,DGBEQU)
#define LAPACKE_cgbequ   LAPACKE_NAME(cgbequ,CGBEQU)
#define LAPACKE_zgbequ   LAPACKE_NAME(zgbequ,ZGBEQU)

#define LAPACKE_sgbequb   LAPACKE_NAME(sgbequb,SGBEQUB)
#define LAPACKE_dgbequb   LAPACKE_NAME(dgbequb,DGBEQUB)
#define LAPACKE_cgbequb   LAPACKE_NAME(cgbequb,CGBEQUB)
#define LAPACKE_zgbequb   LAPACKE_NAME(zgbequb,ZGBEQUB)

#define LAPACKE_sgbrfs   LAPACKE_NAME(sgbrfs,SGBRFS)
#define LAPACKE_dgbrfs   LAPACKE_NAME(dgbrfs,DGBRFS)
#define LAPACKE_cgbrfs   LAPACKE_NAME(cgbrfs,CGBRFS)
#define LAPACKE_zgbrfs   LAPACKE_NAME(zgbrfs,ZGBRFS)

#define LAPACKE_sgbrfsx   LAPACKE_NAME(sgbrfsx,SGBRFSX)
#define LAPACKE_dgbrfsx   LAPACKE_NAME(dgbrfsx,DGBRFSX)
#define LAPACKE_cgbrfsx   LAPACKE_NAME(cgbrfsx,CGBRFSX)
#define LAPACKE_zgbrfsx   LAPACKE_NAME(zgbrfsx,ZGBRFSX)

#define LAPACKE_sgbsv   LAPACKE_NAME(sgbsv,SGBSV)
#define LAPACKE_dgbsv   LAPACKE_NAME(dgbsv,DGBSV)
#define LAPACKE_cgbsv   LAPACKE_NAME(cgbsv,CGBSV)
#define LAPACKE_zgbsv   LAPACKE_NAME(zgbsv,ZGBSV)

#define LAPACKE_sgbsvx   LAPACKE_NAME(sgbsvx,SGBSVX)
#define LAPACKE_dgbsvx   LAPACKE_NAME(dgbsvx,DGBSVX)
#define LAPACKE_cgbsvx   LAPACKE_NAME(cgbsvx,CGBSVX)
#define LAPACKE_zgbsvx   LAPACKE_NAME(zgbsvx,ZGBSVX)

#define LAPACKE_sgbsvxx   LAPACKE_NAME(sgbsvxx,SGBSVXX)
#define LAPACKE_dgbsvxx   LAPACKE_NAME(dgbsvxx,DGBSVXX)
#define LAPACKE_cgbsvxx   LAPACKE_NAME(cgbsvxx,CGBSVXX)
#define LAPACKE_zgbsvxx   LAPACKE_NAME(zgbsvxx,ZGBSVXX)

#define LAPACKE_sgbtrf   LAPACKE_NAME(sgbtrf,SGBTRF)
#define LAPACKE_dgbtrf   LAPACKE_NAME(dgbtrf,DGBTRF)
#define LAPACKE_cgbtrf   LAPACKE_NAME(cgbtrf,CGBTRF)
#define LAPACKE_zgbtrf   LAPACKE_NAME(zgbtrf,ZGBTRF)

#define LAPACKE_sgbtrs   LAPACKE_NAME(sgbtrs,SGBTRS)
#define LAPACKE_dgbtrs   LAPACKE_NAME(dgbtrs,DGBTRS)
#define LAPACKE_cgbtrs   LAPACKE_NAME(cgbtrs,CGBTRS)
#define LAPACKE_zgbtrs   LAPACKE_NAME(zgbtrs,ZGBTRS)

#define LAPACKE_sgebak   LAPACKE_NAME(sgebak,SGEBAK)
#define LAPACKE_dgebak   LAPACKE_NAME(dgebak,DGEBAK)
#define LAPACKE_cgebak   LAPACKE_NAME(cgebak,CGEBAK)
#define LAPACKE_zgebak   LAPACKE_NAME(zgebak,ZGEBAK)

#define LAPACKE_sgebal   LAPACKE_NAME(sgebal,SGEBAL)
#define LAPACKE_dgebal   LAPACKE_NAME(dgebal,DGEBAL)
#define LAPACKE_cgebal   LAPACKE_NAME(cgebal,CGEBAL)
#define LAPACKE_zgebal   LAPACKE_NAME(zgebal,ZGEBAL)

#define LAPACKE_sgebrd   LAPACKE_NAME(sgebrd,SGEBRD)
#define LAPACKE_dgebrd   LAPACKE_NAME(dgebrd,DGEBRD)
#define LAPACKE_cgebrd   LAPACKE_NAME(cgebrd,CGEBRD)
#define LAPACKE_zgebrd   LAPACKE_NAME(zgebrd,ZGEBRD)

#define LAPACKE_sgecon   LAPACKE_NAME(sgecon,SGECON)
#define LAPACKE_dgecon   LAPACKE_NAME(dgecon,DGECON)
#define LAPACKE_cgecon   LAPACKE_NAME(cgecon,CGECON)
#define LAPACKE_zgecon   LAPACKE_NAME(zgecon,ZGECON)

#define LAPACKE_sgeequ   LAPACKE_NAME(sgeequ,SGEEQU)
#define LAPACKE_dgeequ   LAPACKE_NAME(dgeequ,DGEEQU)
#define LAPACKE_cgeequ   LAPACKE_NAME(cgeequ,CGEEQU)
#define LAPACKE_zgeequ   LAPACKE_NAME(zgeequ,ZGEEQU)

#define LAPACKE_sgeequb   LAPACKE_NAME(sgeequb,SGEEQUB)
#define LAPACKE_dgeequb   LAPACKE_NAME(dgeequb,DGEEQUB)
#define LAPACKE_cgeequb   LAPACKE_NAME(cgeequb,CGEEQUB)
#define LAPACKE_zgeequb   LAPACKE_NAME(zgeequb,ZGEEQUB)

#define LAPACKE_sgees   LAPACKE_NAME(sgees,SGEES)
#define LAPACKE_dgees   LAPACKE_NAME(dgees,DGEES)
#define LAPACKE_cgees   LAPACKE_NAME(cgees,CGEES)
#define LAPACKE_zgees   LAPACKE_NAME(zgees,ZGEES)

#define LAPACKE_sgeesx   LAPACKE_NAME(sgeesx,SGEESX)
#define LAPACKE_dgeesx   LAPACKE_NAME(dgeesx,DGEESX)
#define LAPACKE_cgeesx   LAPACKE_NAME(cgeesx,CGEESX)
#define LAPACKE_zgeesx   LAPACKE_NAME(zgeesx,ZGEESX)

#define LAPACKE_sgeev   LAPACKE_NAME(sgeev,SGEEV)
#define LAPACKE_dgeev   LAPACKE_NAME(dgeev,DGEEV)
#define LAPACKE_cgeev   LAPACKE_NAME(cgeev,CGEEV)
#define LAPACKE_zgeev   LAPACKE_NAME(zgeev,ZGEEV)

#define LAPACKE_sgeevx   LAPACKE_NAME(sgeevx,SGEEVX)
#define LAPACKE_dgeevx   LAPACKE_NAME(dgeevx,DGEEVX)
#define LAPACKE_cgeevx   LAPACKE_NAME(cgeevx,CGEEVX)
#define LAPACKE_zgeevx   LAPACKE_NAME(zgeevx,ZGEEVX)

#define LAPACKE_sgehrd   LAPACKE_NAME(sgehrd,SGEHRD)
#define LAPACKE_dgehrd   LAPACKE_NAME(dgehrd,DGEHRD)
#define LAPACKE_cgehrd   LAPACKE_NAME(cgehrd,CGEHRD)
#define LAPACKE_zgehrd   LAPACKE_NAME(zgehrd,ZGEHRD)

#define LAPACKE_sgejsv   LAPACKE_NAME(sgejsv,SGEJSV)
#define LAPACKE_dgejsv   LAPACKE_NAME(dgejsv,DGEJSV)

#define LAPACKE_sgelq2   LAPACKE_NAME(sgelq2,SGELQ2)
#define LAPACKE_dgelq2   LAPACKE_NAME(dgelq2,DGELQ2)
#define LAPACKE_cgelq2   LAPACKE_NAME(cgelq2,CGELQ2)
#define LAPACKE_zgelq2   LAPACKE_NAME(zgelq2,ZGELQ2)

#define LAPACKE_sgelqf   LAPACKE_NAME(sgelqf,SGELQF)
#define LAPACKE_dgelqf   LAPACKE_NAME(dgelqf,DGELQF)
#define LAPACKE_cgelqf   LAPACKE_NAME(cgelqf,CGELQF)
#define LAPACKE_zgelqf   LAPACKE_NAME(zgelqf,ZGELQF)

#define LAPACKE_sgels   LAPACKE_NAME(sgels,SGELS)
#define LAPACKE_dgels   LAPACKE_NAME(dgels,DGELS)
#define LAPACKE_cgels   LAPACKE_NAME(cgels,CGELS)
#define LAPACKE_zgels   LAPACKE_NAME(zgels,ZGELS)

#define LAPACKE_sgelsd   LAPACKE_NAME(sgelsd,SGELSD)
#define LAPACKE_dgelsd   LAPACKE_NAME(dgelsd,DGELSD)
#define LAPACKE_cgelsd   LAPACKE_NAME(cgelsd,CGELSD)
#define LAPACKE_zgelsd   LAPACKE_NAME(zgelsd,ZGELSD)

#define LAPACKE_sgelss   LAPACKE_NAME(sgelss,SGELSS)
#define LAPACKE_dgelss   LAPACKE_NAME(dgelss,DGELSS)
#define LAPACKE_cgelss   LAPACKE_NAME(cgelss,CGELSS)
#define LAPACKE_zgelss   LAPACKE_NAME(zgelss,ZGELSS)

#define LAPACKE_sgelsy   LAPACKE_NAME(sgelsy,SGELSY)
#define LAPACKE_dgelsy   LAPACKE_NAME(dgelsy,DGELSY)
#define LAPACKE_cgelsy   LAPACKE_NAME(cgelsy,CGELSY)
#define LAPACKE_zgelsy   LAPACKE_NAME(zgelsy,ZGELSY)

#define LAPACKE_sgeqlf   LAPACKE_NAME(sgeqlf,SGEQLF)
#define LAPACKE_dgeqlf   LAPACKE_NAME(dgeqlf,DGEQLF)
#define LAPACKE_cgeqlf   LAPACKE_NAME(cgeqlf,CGEQLF)
#define LAPACKE_zgeqlf   LAPACKE_NAME(zgeqlf,ZGEQLF)

#define LAPACKE_sgeqp3   LAPACKE_NAME(sgeqp3,SGEQP3)
#define LAPACKE_dgeqp3   LAPACKE_NAME(dgeqp3,DGEQP3)
#define LAPACKE_cgeqp3   LAPACKE_NAME(cgeqp3,CGEQP3)
#define LAPACKE_zgeqp3   LAPACKE_NAME(zgeqp3,ZGEQP3)

#define LAPACKE_sgeqpf   LAPACKE_NAME(sgeqpf,SGEQPF)
#define LAPACKE_dgeqpf   LAPACKE_NAME(dgeqpf,DGEQPF)
#define LAPACKE_cgeqpf   LAPACKE_NAME(cgeqpf,CGEQPF)
#define LAPACKE_zgeqpf   LAPACKE_NAME(zgeqpf,ZGEQPF)

#define LAPACKE_sgeqr2   LAPACKE_NAME(sgeqr2,SGEQR2)
#define LAPACKE_dgeqr2   LAPACKE_NAME(dgeqr2,DGEQR2)
#define LAPACKE_cgeqr2   LAPACKE_NAME(cgeqr2,CGEQR2)
#define LAPACKE_zgeqr2   LAPACKE_NAME(zgeqr2,ZGEQR2)

#define LAPACKE_sgeqrf   LAPACKE_NAME(sgeqrf,SGEQRF)
#define LAPACKE_dgeqrf   LAPACKE_NAME(dgeqrf,DGEQRF)
#define LAPACKE_cgeqrf   LAPACKE_NAME(cgeqrf,CGEQRF)
#define LAPACKE_zgeqrf   LAPACKE_NAME(zgeqrf,ZGEQRF)

#define LAPACKE_sgeqrfp   LAPACKE_NAME(sgeqrfp,SGEQRFP)
#define LAPACKE_dgeqrfp   LAPACKE_NAME(dgeqrfp,DGEQRFP)
#define LAPACKE_cgeqrfp   LAPACKE_NAME(cgeqrfp,CGEQRFP)
#define LAPACKE_zgeqrfp   LAPACKE_NAME(zgeqrfp,ZGEQRFP)

#define LAPACKE_sgerfs   LAPACKE_NAME(sgerfs,SGERFS)
#define LAPACKE_dgerfs   LAPACKE_NAME(dgerfs,DGERFS)
#define LAPACKE_cgerfs   LAPACKE_NAME(cgerfs,CGERFS)
#define LAPACKE_zgerfs   LAPACKE_NAME(zgerfs,ZGERFS)

#define LAPACKE_sgerfsx   LAPACKE_NAME(sgerfsx,SGERFSX)
#define LAPACKE_dgerfsx   LAPACKE_NAME(dgerfsx,DGERFSX)
#define LAPACKE_cgerfsx   LAPACKE_NAME(cgerfsx,CGERFSX)
#define LAPACKE_zgerfsx   LAPACKE_NAME(zgerfsx,ZGERFSX)

#define LAPACKE_sgerqf   LAPACKE_NAME(sgerqf,SGERQF)
#define LAPACKE_dgerqf   LAPACKE_NAME(dgerqf,DGERQF)
#define LAPACKE_cgerqf   LAPACKE_NAME(cgerqf,CGERQF)
#define LAPACKE_zgerqf   LAPACKE_NAME(zgerqf,ZGERQF)

#define LAPACKE_sgesdd   LAPACKE_NAME(sgesdd,SGESDD)
#define LAPACKE_dgesdd   LAPACKE_NAME(dgesdd,DGESDD)
#define LAPACKE_cgesdd   LAPACKE_NAME(cgesdd,CGESDD)
#define LAPACKE_zgesdd   LAPACKE_NAME(zgesdd,ZGESDD)

#define LAPACKE_sgesv   LAPACKE_NAME(sgesv,SGESV)
#define LAPACKE_dgesv   LAPACKE_NAME(dgesv,DGESV)
#define LAPACKE_cgesv   LAPACKE_NAME(cgesv,CGESV)
#define LAPACKE_zgesv   LAPACKE_NAME(zgesv,ZGESV)
#define LAPACKE_dsgesv   LAPACKE_NAME(dsgesv,DSGESV)
#define LAPACKE_zcgesv   LAPACKE_NAME(zcgesv,ZCGESV)

#define LAPACKE_sgesvd   LAPACKE_NAME(sgesvd,SGESVD)
#define LAPACKE_dgesvd   LAPACKE_NAME(dgesvd,DGESVD)
#define LAPACKE_cgesvd   LAPACKE_NAME(cgesvd,CGESVD)
#define LAPACKE_zgesvd   LAPACKE_NAME(zgesvd,ZGESVD)

#define LAPACKE_sgesvj   LAPACKE_NAME(sgesvj,SGESVJ)
#define LAPACKE_dgesvj   LAPACKE_NAME(dgesvj,DGESVJ)

#define LAPACKE_sgesvx   LAPACKE_NAME(sgesvx,SGESVX)
#define LAPACKE_dgesvx   LAPACKE_NAME(dgesvx,DGESVX)
#define LAPACKE_cgesvx   LAPACKE_NAME(cgesvx,CGESVX)
#define LAPACKE_zgesvx   LAPACKE_NAME(zgesvx,ZGESVX)

#define LAPACKE_sgesvxx   LAPACKE_NAME(sgesvxx,SGESVXX)
#define LAPACKE_dgesvxx   LAPACKE_NAME(dgesvxx,DGESVXX)
#define LAPACKE_cgesvxx   LAPACKE_NAME(cgesvxx,CGESVXX)
#define LAPACKE_zgesvxx   LAPACKE_NAME(zgesvxx,ZGESVXX)

#define LAPACKE_sgetf2   LAPACKE_NAME(sgetf2,SGETF2)
#define LAPACKE_dgetf2   LAPACKE_NAME(dgetf2,DGETF2)
#define LAPACKE_cgetf2   LAPACKE_NAME(cgetf2,CGETF2)
#define LAPACKE_zgetf2   LAPACKE_NAME(zgetf2,ZGETF2)

#define LAPACKE_sgetrf   LAPACKE_NAME(sgetrf,SGETRF)
#define LAPACKE_dgetrf   LAPACKE_NAME(dgetrf,DGETRF)
#define LAPACKE_cgetrf   LAPACKE_NAME(cgetrf,CGETRF)
#define LAPACKE_zgetrf   LAPACKE_NAME(zgetrf,ZGETRF)

#define LAPACKE_sgetri   LAPACKE_NAME(sgetri,SGETRI)
#define LAPACKE_dgetri   LAPACKE_NAME(dgetri,DGETRI)
#define LAPACKE_cgetri   LAPACKE_NAME(cgetri,CGETRI)
#define LAPACKE_zgetri   LAPACKE_NAME(zgetri,ZGETRI)

#define LAPACKE_sgetrs   LAPACKE_NAME(sgetrs,SGETRS)
#define LAPACKE_dgetrs   LAPACKE_NAME(dgetrs,DGETRS)
#define LAPACKE_cgetrs   LAPACKE_NAME(cgetrs,CGETRS)
#define LAPACKE_zgetrs   LAPACKE_NAME(zgetrs,ZGETRS)

#define LAPACKE_sggbak   LAPACKE_NAME(sggbak,SGGBAK)
#define LAPACKE_dggbak   LAPACKE_NAME(dggbak,DGGBAK)
#define LAPACKE_cggbak   LAPACKE_NAME(cggbak,CGGBAK)
#define LAPACKE_zggbak   LAPACKE_NAME(zggbak,ZGGBAK)

#define LAPACKE_sggbal   LAPACKE_NAME(sggbal,SGGBAL)
#define LAPACKE_dggbal   LAPACKE_NAME(dggbal,DGGBAL)
#define LAPACKE_cggbal   LAPACKE_NAME(cggbal,CGGBAL)
#define LAPACKE_zggbal   LAPACKE_NAME(zggbal,ZGGBAL)

#define LAPACKE_sgges   LAPACKE_NAME(sgges,SGGES)
#define LAPACKE_dgges   LAPACKE_NAME(dgges,DGGES)
#define LAPACKE_cgges   LAPACKE_NAME(cgges,CGGES)
#define LAPACKE_zgges   LAPACKE_NAME(zgges,ZGGES)

#define LAPACKE_sggesx   LAPACKE_NAME(sggesx,SGGESX)
#define LAPACKE_dggesx   LAPACKE_NAME(dggesx,DGGESX)
#define LAPACKE_cggesx   LAPACKE_NAME(cggesx,CGGESX)
#define LAPACKE_zggesx   LAPACKE_NAME(zggesx,ZGGESX)

#define LAPACKE_sggev   LAPACKE_NAME(sggev,SGGEV)
#define LAPACKE_dggev   LAPACKE_NAME(dggev,DGGEV)
#define LAPACKE_cggev   LAPACKE_NAME(cggev,CGGEV)
#define LAPACKE_zggev   LAPACKE_NAME(zggev,ZGGEV)

#define LAPACKE_sggevx   LAPACKE_NAME(sggevx,SGGEVX)
#define LAPACKE_dggevx   LAPACKE_NAME(dggevx,DGGEVX)
#define LAPACKE_cggevx   LAPACKE_NAME(cggevx,CGGEVX)
#define LAPACKE_zggevx   LAPACKE_NAME(zggevx,ZGGEVX)

#define LAPACKE_sggglm   LAPACKE_NAME(sggglm,SGGGLM)
#define LAPACKE_dggglm   LAPACKE_NAME(dggglm,DGGGLM)
#define LAPACKE_cggglm   LAPACKE_NAME(cggglm,CGGGLM)
#define LAPACKE_zggglm   LAPACKE_NAME(zggglm,ZGGGLM)

#define LAPACKE_sgghrd   LAPACKE_NAME(sgghrd,SGGHRD)
#define LAPACKE_dgghrd   LAPACKE_NAME(dgghrd,DGGHRD)
#define LAPACKE_cgghrd   LAPACKE_NAME(cgghrd,CGGHRD)
#define LAPACKE_zgghrd   LAPACKE_NAME(zgghrd,ZGGHRD)

#define LAPACKE_sgglse   LAPACKE_NAME(sgglse,SGGLSE)
#define LAPACKE_dgglse   LAPACKE_NAME(dgglse,DGGLSE)
#define LAPACKE_cgglse   LAPACKE_NAME(cgglse,CGGLSE)
#define LAPACKE_zgglse   LAPACKE_NAME(zgglse,ZGGLSE)

#define LAPACKE_sggqrf   LAPACKE_NAME(sggqrf,SGGQRF)
#define LAPACKE_dggqrf   LAPACKE_NAME(dggqrf,DGGQRF)
#define LAPACKE_cggqrf   LAPACKE_NAME(cggqrf,CGGQRF)
#define LAPACKE_zggqrf   LAPACKE_NAME(zggqrf,ZGGQRF)

#define LAPACKE_sggrqf   LAPACKE_NAME(sggrqf,SGGRQF)
#define LAPACKE_dggrqf   LAPACKE_NAME(dggrqf,DGGRQF)
#define LAPACKE_cggrqf   LAPACKE_NAME(cggrqf,CGGRQF)
#define LAPACKE_zggrqf   LAPACKE_NAME(zggrqf,ZGGRQF)

#define LAPACKE_sggsvd   LAPACKE_NAME(sggsvd,SGGSVD)
#define LAPACKE_dggsvd   LAPACKE_NAME(dggsvd,DGGSVD)
#define LAPACKE_cggsvd   LAPACKE_NAME(cggsvd,CGGSVD)
#define LAPACKE_zggsvd   LAPACKE_NAME(zggsvd,ZGGSVD)

#define LAPACKE_sggsvp   LAPACKE_NAME(sggsvp,SGGSVP)
#define LAPACKE_dggsvp   LAPACKE_NAME(dggsvp,DGGSVP)
#define LAPACKE_cggsvp   LAPACKE_NAME(cggsvp,CGGSVP)
#define LAPACKE_zggsvp   LAPACKE_NAME(zggsvp,ZGGSVP)

#define LAPACKE_sgtcon   LAPACKE_NAME(sgtcon,SGTCON)
#define LAPACKE_dgtcon   LAPACKE_NAME(dgtcon,DGTCON)
#define LAPACKE_cgtcon   LAPACKE_NAME(cgtcon,CGTCON)
#define LAPACKE_zgtcon   LAPACKE_NAME(zgtcon,ZGTCON)

#define LAPACKE_sgtrfs   LAPACKE_NAME(sgtrfs,SGTRFS)
#define LAPACKE_dgtrfs   LAPACKE_NAME(dgtrfs,DGTRFS)
#define LAPACKE_cgtrfs   LAPACKE_NAME(cgtrfs,CGTRFS)
#define LAPACKE_zgtrfs   LAPACKE_NAME(zgtrfs,ZGTRFS)

#define LAPACKE_sgtsv   LAPACKE_NAME(sgtsv,SGTSV)
#define LAPACKE_dgtsv   LAPACKE_NAME(dgtsv,DGTSV)
#define LAPACKE_cgtsv   LAPACKE_NAME(cgtsv,CGTSV)
#define LAPACKE_zgtsv   LAPACKE_NAME(zgtsv,ZGTSV)

#define LAPACKE_sgtsvx   LAPACKE_NAME(sgtsvx,SGTSVX)
#define LAPACKE_dgtsvx   LAPACKE_NAME(dgtsvx,DGTSVX)
#define LAPACKE_cgtsvx   LAPACKE_NAME(cgtsvx,CGTSVX)
#define LAPACKE_zgtsvx   LAPACKE_NAME(zgtsvx,ZGTSVX)

#define LAPACKE_sgttrf   LAPACKE_NAME(sgttrf,SGTTRF)
#define LAPACKE_dgttrf   LAPACKE_NAME(dgttrf,DGTTRF)
#define LAPACKE_cgttrf   LAPACKE_NAME(cgttrf,CGTTRF)
#define LAPACKE_zgttrf   LAPACKE_NAME(zgttrf,ZGTTRF)

#define LAPACKE_sgttrs   LAPACKE_NAME(sgttrs,SGTTRS)
#define LAPACKE_dgttrs   LAPACKE_NAME(dgttrs,DGTTRS)
#define LAPACKE_cgttrs   LAPACKE_NAME(cgttrs,CGTTRS)
#define LAPACKE_zgttrs   LAPACKE_NAME(zgttrs,ZGTTRS)

#define LAPACKE_chbev   LAPACKE_NAME(chbev,CHBEV)
#define LAPACKE_zhbev   LAPACKE_NAME(zhbev,ZHBEV)

#define LAPACKE_chbevd   LAPACKE_NAME(chbevd,CHBEVD)
#define LAPACKE_zhbevd   LAPACKE_NAME(zhbevd,ZHBEVD)

#define LAPACKE_chbevx   LAPACKE_NAME(chbevx,CHBEVX)
#define LAPACKE_zhbevx   LAPACKE_NAME(zhbevx,ZHBEVX)

#define LAPACKE_chbgst   LAPACKE_NAME(chbgst,CHBGST)
#define LAPACKE_zhbgst   LAPACKE_NAME(zhbgst,ZHBGST)

#define LAPACKE_chbgv   LAPACKE_NAME(chbgv,CHBGV)
#define LAPACKE_zhbgv   LAPACKE_NAME(zhbgv,ZHBGV)

#define LAPACKE_chbgvd   LAPACKE_NAME(chbgvd,CHBGVD)
#define LAPACKE_zhbgvd   LAPACKE_NAME(zhbgvd,ZHBGVD)

#define LAPACKE_chbgvx   LAPACKE_NAME(chbgvx,CHBGVX)
#define LAPACKE_zhbgvx   LAPACKE_NAME(zhbgvx,ZHBGVX)

#define LAPACKE_chbtrd   LAPACKE_NAME(chbtrd,CHBTRD)
#define LAPACKE_zhbtrd   LAPACKE_NAME(zhbtrd,ZHBTRD)

#define LAPACKE_checon   LAPACKE_NAME(checon,CHECON)
#define LAPACKE_zhecon   LAPACKE_NAME(zhecon,ZHECON)

#define LAPACKE_cheequb   LAPACKE_NAME(cheequb,CHEEQUB)
#define LAPACKE_zheequb   LAPACKE_NAME(zheequb,ZHEEQUB)

#define LAPACKE_cheev   LAPACKE_NAME(cheev,CHEEV)
#define LAPACKE_zheev   LAPACKE_NAME(zheev,ZHEEV)

#define LAPACKE_cheevd   LAPACKE_NAME(cheevd,CHEEVD)
#define LAPACKE_zheevd   LAPACKE_NAME(zheevd,ZHEEVD)

#define LAPACKE_cheevr   LAPACKE_NAME(cheevr,CHEEVR)
#define LAPACKE_zheevr   LAPACKE_NAME(zheevr,ZHEEVR)

#define LAPACKE_cheevx   LAPACKE_NAME(cheevx,CHEEVX)
#define LAPACKE_zheevx   LAPACKE_NAME(zheevx,ZHEEVX)

#define LAPACKE_chegst   LAPACKE_NAME(chegst,CHEGST)
#define LAPACKE_zhegst   LAPACKE_NAME(zhegst,ZHEGST)

#define LAPACKE_chegv   LAPACKE_NAME(chegv,CHEGV)
#define LAPACKE_zhegv   LAPACKE_NAME(zhegv,ZHEGV)

#define LAPACKE_chegvd   LAPACKE_NAME(chegvd,CHEGVD)
#define LAPACKE_zhegvd   LAPACKE_NAME(zhegvd,ZHEGVD)

#define LAPACKE_chegvx   LAPACKE_NAME(chegvx,CHEGVX)
#define LAPACKE_zhegvx   LAPACKE_NAME(zhegvx,ZHEGVX)

#define LAPACKE_cherfs   LAPACKE_NAME(cherfs,CHERFS)
#define LAPACKE_zherfs   LAPACKE_NAME(zherfs,ZHERFS)

#define LAPACKE_cherfsx   LAPACKE_NAME(cherfsx,CHERFSX)
#define LAPACKE_zherfsx   LAPACKE_NAME(zherfsx,ZHERFSX)

#define LAPACKE_chesv   LAPACKE_NAME(chesv,CHESV)
#define LAPACKE_zhesv   LAPACKE_NAME(zhesv,ZHESV)

#define LAPACKE_chesvx   LAPACKE_NAME(chesvx,CHESVX)
#define LAPACKE_zhesvx   LAPACKE_NAME(zhesvx,ZHESVX)

#define LAPACKE_chesvxx   LAPACKE_NAME(chesvxx,CHESVXX)
#define LAPACKE_zhesvxx   LAPACKE_NAME(zhesvxx,ZHESVXX)

#define LAPACKE_chetrd   LAPACKE_NAME(chetrd,CHETRD)
#define LAPACKE_zhetrd   LAPACKE_NAME(zhetrd,ZHETRD)

#define LAPACKE_chetrf   LAPACKE_NAME(chetrf,CHETRF)
#define LAPACKE_zhetrf   LAPACKE_NAME(zhetrf,ZHETRF)

#define LAPACKE_chetri   LAPACKE_NAME(chetri,CHETRI)
#define LAPACKE_zhetri   LAPACKE_NAME(zhetri,ZHETRI)

#define LAPACKE_chetrs   LAPACKE_NAME(chetrs,CHETRS)
#define LAPACKE_zhetrs   LAPACKE_NAME(zhetrs,ZHETRS)

#define LAPACKE_chfrk   LAPACKE_NAME(chfrk,CHFRK)
#define LAPACKE_zhfrk   LAPACKE_NAME(zhfrk,ZHFRK)

#define LAPACKE_shgeqz   LAPACKE_NAME(shgeqz,SHGEQZ)
#define LAPACKE_dhgeqz   LAPACKE_NAME(dhgeqz,DHGEQZ)
#define LAPACKE_chgeqz   LAPACKE_NAME(chgeqz,CHGEQZ)
#define LAPACKE_zhgeqz   LAPACKE_NAME(zhgeqz,ZHGEQZ)

#define LAPACKE_chpcon   LAPACKE_NAME(chpcon,CHPCON)
#define LAPACKE_zhpcon   LAPACKE_NAME(zhpcon,ZHPCON)

#define LAPACKE_chpev   LAPACKE_NAME(chpev,CHPEV)
#define LAPACKE_zhpev   LAPACKE_NAME(zhpev,ZHPEV)

#define LAPACKE_chpevd   LAPACKE_NAME(chpevd,CHPEVD)
#define LAPACKE_zhpevd   LAPACKE_NAME(zhpevd,ZHPEVD)

#define LAPACKE_chpevx   LAPACKE_NAME(chpevx,CHPEVX)
#define LAPACKE_zhpevx   LAPACKE_NAME(zhpevx,ZHPEVX)

#define LAPACKE_chpgst   LAPACKE_NAME(chpgst,CHPGST)
#define LAPACKE_zhpgst   LAPACKE_NAME(zhpgst,ZHPGST)

#define LAPACKE_chpgv   LAPACKE_NAME(chpgv,CHPGV)
#define LAPACKE_zhpgv   LAPACKE_NAME(zhpgv,ZHPGV)

#define LAPACKE_chpgvd   LAPACKE_NAME(chpgvd,CHPGVD)
#define LAPACKE_zhpgvd   LAPACKE_NAME(zhpgvd,ZHPGVD)

#define LAPACKE_chpgvx   LAPACKE_NAME(chpgvx,CHPGVX)
#define LAPACKE_zhpgvx   LAPACKE_NAME(zhpgvx,ZHPGVX)

#define LAPACKE_chprfs   LAPACKE_NAME(chprfs,CHPRFS)
#define LAPACKE_zhprfs   LAPACKE_NAME(zhprfs,ZHPRFS)

#define LAPACKE_chpsv   LAPACKE_NAME(chpsv,CHPSV)
#define LAPACKE_zhpsv   LAPACKE_NAME(zhpsv,ZHPSV)

#define LAPACKE_chpsvx   LAPACKE_NAME(chpsvx,CHPSVX)
#define LAPACKE_zhpsvx   LAPACKE_NAME(zhpsvx,ZHPSVX)

#define LAPACKE_chptrd   LAPACKE_NAME(chptrd,CHPTRD)
#define LAPACKE_zhptrd   LAPACKE_NAME(zhptrd,ZHPTRD)

#define LAPACKE_chptrf   LAPACKE_NAME(chptrf,CHPTRF)
#define LAPACKE_zhptrf   LAPACKE_NAME(zhptrf,ZHPTRF)

#define LAPACKE_chptri   LAPACKE_NAME(chptri,CHPTRI)
#define LAPACKE_zhptri   LAPACKE_NAME(zhptri,ZHPTRI)

#define LAPACKE_chptrs   LAPACKE_NAME(chptrs,CHPTRS)
#define LAPACKE_zhptrs   LAPACKE_NAME(zhptrs,ZHPTRS)

#define LAPACKE_shsein   LAPACKE_NAME(shsein,SHSEIN)
#define LAPACKE_dhsein   LAPACKE_NAME(dhsein,DHSEIN)
#define LAPACKE_chsein   LAPACKE_NAME(chsein,CHSEIN)
#define LAPACKE_zhsein   LAPACKE_NAME(zhsein,ZHSEIN)

#define LAPACKE_shseqr   LAPACKE_NAME(shseqr,SHSEQR)
#define LAPACKE_dhseqr   LAPACKE_NAME(dhseqr,DHSEQR)
#define LAPACKE_chseqr   LAPACKE_NAME(chseqr,CHSEQR)
#define LAPACKE_zhseqr   LAPACKE_NAME(zhseqr,ZHSEQR)

#define LAPACKE_clacgv   LAPACKE_NAME(clacgv,CLACGV)
#define LAPACKE_zlacgv   LAPACKE_NAME(zlacgv,ZLACGV)

#define LAPACKE_slacpy   LAPACKE_NAME(slacpy,SLACPY)
#define LAPACKE_dlacpy   LAPACKE_NAME(dlacpy,DLACPY)
#define LAPACKE_clacpy   LAPACKE_NAME(clacpy,CLACPY)
#define LAPACKE_zlacpy   LAPACKE_NAME(zlacpy,ZLACPY)

#define LAPACKE_zlag2c   LAPACKE_NAME(zlag2c,ZLAG2C)

#define LAPACKE_slag2d   LAPACKE_NAME(slag2d,SLAG2D)

#define LAPACKE_dlag2s   LAPACKE_NAME(dlag2s,DLAG2S)

#define LAPACKE_clag2z   LAPACKE_NAME(clag2z,CLAG2Z)

#define LAPACKE_slagge   LAPACKE_NAME(slagge,SLAGGE)
#define LAPACKE_dlagge   LAPACKE_NAME(dlagge,DLAGGE)
#define LAPACKE_clagge   LAPACKE_NAME(clagge,CLAGGE)
#define LAPACKE_zlagge   LAPACKE_NAME(zlagge,ZLAGGE)

#define LAPACKE_slamch   LAPACKE_NAME(slamch,SLAMCH)
#define LAPACKE_dlamch   LAPACKE_NAME(dlamch,DLAMCH)

#define LAPACKE_slange   LAPACKE_NAME(slange,SLANGE)
#define LAPACKE_dlange   LAPACKE_NAME(dlange,DLANGE)
#define LAPACKE_clange   LAPACKE_NAME(clange,CLANGE)
#define LAPACKE_zlange   LAPACKE_NAME(zlange,ZLANGE)

#define LAPACKE_clanhe   LAPACKE_NAME(clanhe,CLANHE)
#define LAPACKE_zlanhe   LAPACKE_NAME(zlanhe,ZLANHE)

#define LAPACKE_slansy   LAPACKE_NAME(slansy,SLANSY)
#define LAPACKE_dlansy   LAPACKE_NAME(dlansy,DLANSY)
#define LAPACKE_clansy   LAPACKE_NAME(clansy,CLANSY)
#define LAPACKE_zlansy   LAPACKE_NAME(zlansy,ZLANSY)

#define LAPACKE_slantr   LAPACKE_NAME(slantr,SLANTR)
#define LAPACKE_dlantr   LAPACKE_NAME(dlantr,DLANTR)
#define LAPACKE_clantr   LAPACKE_NAME(clantr,CLANTR)
#define LAPACKE_zlantr   LAPACKE_NAME(zlantr,ZLANTR)

#define LAPACKE_slarfb   LAPACKE_NAME(slarfb,SLARFB)
#define LAPACKE_dlarfb   LAPACKE_NAME(dlarfb,DLARFB)
#define LAPACKE_clarfb   LAPACKE_NAME(clarfb,CLARFB)
#define LAPACKE_zlarfb   LAPACKE_NAME(zlarfb,ZLARFB)

#define LAPACKE_slarfg   LAPACKE_NAME(slarfg,SLARFG)
#define LAPACKE_dlarfg   LAPACKE_NAME(dlarfg,DLARFG)
#define LAPACKE_clarfg   LAPACKE_NAME(clarfg,CLARFG)
#define LAPACKE_zlarfg   LAPACKE_NAME(zlarfg,ZLARFG)

#define LAPACKE_slarft   LAPACKE_NAME(slarft,SLARFT)
#define LAPACKE_dlarft   LAPACKE_NAME(dlarft,DLARFT)
#define LAPACKE_clarft   LAPACKE_NAME(clarft,CLARFT)
#define LAPACKE_zlarft   LAPACKE_NAME(zlarft,ZLARFT)

#define LAPACKE_slarfx   LAPACKE_NAME(slarfx,SLARFX)
#define LAPACKE_dlarfx   LAPACKE_NAME(dlarfx,DLARFX)
#define LAPACKE_clarfx   LAPACKE_NAME(clarfx,CLARFX)
#define LAPACKE_zlarfx   LAPACKE_NAME(zlarfx,ZLARFX)

#define LAPACKE_slarnv   LAPACKE_NAME(slarnv,SLARNV)
#define LAPACKE_dlarnv   LAPACKE_NAME(dlarnv,DLARNV)
#define LAPACKE_clarnv   LAPACKE_NAME(clarnv,CLARNV)
#define LAPACKE_zlarnv   LAPACKE_NAME(zlarnv,ZLARNV)

#define LAPACKE_slaset   LAPACKE_NAME(slaset,SLASET)
#define LAPACKE_dlaset   LAPACKE_NAME(dlaset,DLASET)
#define LAPACKE_claset   LAPACKE_NAME(claset,CLASET)
#define LAPACKE_zlaset   LAPACKE_NAME(zlaset,ZLASET)

#define LAPACKE_slasrt   LAPACKE_NAME(slasrt,SLASRT)
#define LAPACKE_dlasrt   LAPACKE_NAME(dlasrt,DLASRT)

#define LAPACKE_slaswp   LAPACKE_NAME(slaswp,SLASWP)
#define LAPACKE_dlaswp   LAPACKE_NAME(dlaswp,DLASWP)
#define LAPACKE_claswp   LAPACKE_NAME(claswp,CLASWP)
#define LAPACKE_zlaswp   LAPACKE_NAME(zlaswp,ZLASWP)

#define LAPACKE_slatms   LAPACKE_NAME(slatms,SLATMS)
#define LAPACKE_dlatms   LAPACKE_NAME(dlatms,DLATMS)
#define LAPACKE_clatms   LAPACKE_NAME(clatms,CLATMS)
#define LAPACKE_zlatms   LAPACKE_NAME(zlatms,ZLATMS)

#define LAPACKE_slauum   LAPACKE_NAME(slauum,SLAUUM)
#define LAPACKE_dlauum   LAPACKE_NAME(dlauum,DLAUUM)
#define LAPACKE_clauum   LAPACKE_NAME(clauum,CLAUUM)
#define LAPACKE_zlauum   LAPACKE_NAME(zlauum,ZLAUUM)

#define LAPACKE_sopgtr   LAPACKE_NAME(sopgtr,SOPGTR)
#define LAPACKE_dopgtr   LAPACKE_NAME(dopgtr,DOPGTR)

#define LAPACKE_sopmtr   LAPACKE_NAME(sopmtr,SOPMTR)
#define LAPACKE_dopmtr   LAPACKE_NAME(dopmtr,DOPMTR)

#define LAPACKE_sorgbr   LAPACKE_NAME(sorgbr,SORGBR)
#define LAPACKE_dorgbr   LAPACKE_NAME(dorgbr,DORGBR)

#define LAPACKE_sorghr   LAPACKE_NAME(sorghr,SORGHR)
#define LAPACKE_dorghr   LAPACKE_NAME(dorghr,DORGHR)

#define LAPACKE_sorglq   LAPACKE_NAME(sorglq,SORGLQ)
#define LAPACKE_dorglq   LAPACKE_NAME(dorglq,DORGLQ)

#define LAPACKE_sorgql   LAPACKE_NAME(sorgql,SORGQL)
#define LAPACKE_dorgql   LAPACKE_NAME(dorgql,DORGQL)

#define LAPACKE_sorgqr   LAPACKE_NAME(sorgqr,SORGQR)
#define LAPACKE_dorgqr   LAPACKE_NAME(dorgqr,DORGQR)

#define LAPACKE_sorgrq   LAPACKE_NAME(sorgrq,SORGRQ)
#define LAPACKE_dorgrq   LAPACKE_NAME(dorgrq,DORGRQ)

#define LAPACKE_sorgtr   LAPACKE_NAME(sorgtr,SORGTR)
#define LAPACKE_dorgtr   LAPACKE_NAME(dorgtr,DORGTR)

#define LAPACKE_sormbr   LAPACKE_NAME(sormbr,SORMBR)
#define LAPACKE_dormbr   LAPACKE_NAME(dormbr,DORMBR)

#define LAPACKE_sormhr   LAPACKE_NAME(sormhr,SORMHR)
#define LAPACKE_dormhr   LAPACKE_NAME(dormhr,DORMHR)

#define LAPACKE_sormlq   LAPACKE_NAME(sormlq,SORMLQ)
#define LAPACKE_dormlq   LAPACKE_NAME(dormlq,DORMLQ)

#define LAPACKE_sormql   LAPACKE_NAME(sormql,SORMQL)
#define LAPACKE_dormql   LAPACKE_NAME(dormql,DORMQL)

#define LAPACKE_sormqr   LAPACKE_NAME(sormqr,SORMQR)
#define LAPACKE_dormqr   LAPACKE_NAME(dormqr,DORMQR)

#define LAPACKE_sormrq   LAPACKE_NAME(sormrq,SORMRQ)
#define LAPACKE_dormrq   LAPACKE_NAME(dormrq,DORMRQ)

#define LAPACKE_sormrz   LAPACKE_NAME(sormrz,SORMRZ)
#define LAPACKE_dormrz   LAPACKE_NAME(dormrz,DORMRZ)

#define LAPACKE_sormtr   LAPACKE_NAME(sormtr,SORMTR)
#define LAPACKE_dormtr   LAPACKE_NAME(dormtr,DORMTR)

#define LAPACKE_spbcon   LAPACKE_NAME(spbcon,SPBCON)
#define LAPACKE_dpbcon   LAPACKE_NAME(dpbcon,DPBCON)
#define LAPACKE_cpbcon   LAPACKE_NAME(cpbcon,CPBCON)
#define LAPACKE_zpbcon   LAPACKE_NAME(zpbcon,ZPBCON)

#define LAPACKE_spbequ   LAPACKE_NAME(spbequ,SPBEQU)
#define LAPACKE_dpbequ   LAPACKE_NAME(dpbequ,DPBEQU)
#define LAPACKE_cpbequ   LAPACKE_NAME(cpbequ,CPBEQU)
#define LAPACKE_zpbequ   LAPACKE_NAME(zpbequ,ZPBEQU)

#define LAPACKE_spbrfs   LAPACKE_NAME(spbrfs,SPBRFS)
#define LAPACKE_dpbrfs   LAPACKE_NAME(dpbrfs,DPBRFS)
#define LAPACKE_cpbrfs   LAPACKE_NAME(cpbrfs,CPBRFS)
#define LAPACKE_zpbrfs   LAPACKE_NAME(zpbrfs,ZPBRFS)

#define LAPACKE_spbstf   LAPACKE_NAME(spbstf,SPBSTF)
#define LAPACKE_dpbstf   LAPACKE_NAME(dpbstf,DPBSTF)
#define LAPACKE_cpbstf   LAPACKE_NAME(cpbstf,CPBSTF)
#define LAPACKE_zpbstf   LAPACKE_NAME(zpbstf,ZPBSTF)

#define LAPACKE_spbsv   LAPACKE_NAME(spbsv,SPBSV)
#define LAPACKE_dpbsv   LAPACKE_NAME(dpbsv,DPBSV)
#define LAPACKE_cpbsv   LAPACKE_NAME(cpbsv,CPBSV)
#define LAPACKE_zpbsv   LAPACKE_NAME(zpbsv,ZPBSV)

#define LAPACKE_spbsvx   LAPACKE_NAME(spbsvx,SPBSVX)
#define LAPACKE_dpbsvx   LAPACKE_NAME(dpbsvx,DPBSVX)
#define LAPACKE_cpbsvx   LAPACKE_NAME(cpbsvx,CPBSVX)
#define LAPACKE_zpbsvx   LAPACKE_NAME(zpbsvx,ZPBSVX)

#define LAPACKE_spbtrf   LAPACKE_NAME(spbtrf,SPBTRF)
#define LAPACKE_dpbtrf   LAPACKE_NAME(dpbtrf,DPBTRF)
#define LAPACKE_cpbtrf   LAPACKE_NAME(cpbtrf,CPBTRF)
#define LAPACKE_zpbtrf   LAPACKE_NAME(zpbtrf,ZPBTRF)

#define LAPACKE_spbtrs   LAPACKE_NAME(spbtrs,SPBTRS)
#define LAPACKE_dpbtrs   LAPACKE_NAME(dpbtrs,DPBTRS)
#define LAPACKE_cpbtrs   LAPACKE_NAME(cpbtrs,CPBTRS)
#define LAPACKE_zpbtrs   LAPACKE_NAME(zpbtrs,ZPBTRS)

#define LAPACKE_spftrf   LAPACKE_NAME(spftrf,SPFTRF)
#define LAPACKE_dpftrf   LAPACKE_NAME(dpftrf,DPFTRF)
#define LAPACKE_cpftrf   LAPACKE_NAME(cpftrf,CPFTRF)
#define LAPACKE_zpftrf   LAPACKE_NAME(zpftrf,ZPFTRF)

#define LAPACKE_spftri   LAPACKE_NAME(spftri,SPFTRI)
#define LAPACKE_dpftri   LAPACKE_NAME(dpftri,DPFTRI)
#define LAPACKE_cpftri   LAPACKE_NAME(cpftri,CPFTRI)
#define LAPACKE_zpftri   LAPACKE_NAME(zpftri,ZPFTRI)

#define LAPACKE_spftrs   LAPACKE_NAME(spftrs,SPFTRS)
#define LAPACKE_dpftrs   LAPACKE_NAME(dpftrs,DPFTRS)
#define LAPACKE_cpftrs   LAPACKE_NAME(cpftrs,CPFTRS)
#define LAPACKE_zpftrs   LAPACKE_NAME(zpftrs,ZPFTRS)

#define LAPACKE_spocon   LAPACKE_NAME(spocon,SPOCON)
#define LAPACKE_dpocon   LAPACKE_NAME(dpocon,DPOCON)
#define LAPACKE_cpocon   LAPACKE_NAME(cpocon,CPOCON)
#define LAPACKE_zpocon   LAPACKE_NAME(zpocon,ZPOCON)

#define LAPACKE_spoequ   LAPACKE_NAME(spoequ,SPOEQU)
#define LAPACKE_dpoequ   LAPACKE_NAME(dpoequ,DPOEQU)
#define LAPACKE_cpoequ   LAPACKE_NAME(cpoequ,CPOEQU)
#define LAPACKE_zpoequ   LAPACKE_NAME(zpoequ,ZPOEQU)

#define LAPACKE_spoequb   LAPACKE_NAME(spoequb,SPOEQUB)
#define LAPACKE_dpoequb   LAPACKE_NAME(dpoequb,DPOEQUB)
#define LAPACKE_cpoequb   LAPACKE_NAME(cpoequb,CPOEQUB)
#define LAPACKE_zpoequb   LAPACKE_NAME(zpoequb,ZPOEQUB)

#define LAPACKE_sporfs   LAPACKE_NAME(sporfs,SPORFS)
#define LAPACKE_dporfs   LAPACKE_NAME(dporfs,DPORFS)
#define LAPACKE_cporfs   LAPACKE_NAME(cporfs,CPORFS)
#define LAPACKE_zporfs   LAPACKE_NAME(zporfs,ZPORFS)

#define LAPACKE_sporfsx   LAPACKE_NAME(sporfsx,SPORFSX)
#define LAPACKE_dporfsx   LAPACKE_NAME(dporfsx,DPORFSX)
#define LAPACKE_cporfsx   LAPACKE_NAME(cporfsx,CPORFSX)
#define LAPACKE_zporfsx   LAPACKE_NAME(zporfsx,ZPORFSX)

#define LAPACKE_sposv   LAPACKE_NAME(sposv,SPOSV)
#define LAPACKE_dposv   LAPACKE_NAME(dposv,DPOSV)
#define LAPACKE_cposv   LAPACKE_NAME(cposv,CPOSV)
#define LAPACKE_zposv   LAPACKE_NAME(zposv,ZPOSV)
#define LAPACKE_dsposv   LAPACKE_NAME(dsposv,DSPOSV)
#define LAPACKE_zcposv   LAPACKE_NAME(zcposv,ZCPOSV)

#define LAPACKE_sposvx   LAPACKE_NAME(sposvx,SPOSVX)
#define LAPACKE_dposvx   LAPACKE_NAME(dposvx,DPOSVX)
#define LAPACKE_cposvx   LAPACKE_NAME(cposvx,CPOSVX)
#define LAPACKE_zposvx   LAPACKE_NAME(zposvx,ZPOSVX)

#define LAPACKE_sposvxx   LAPACKE_NAME(sposvxx,SPOSVXX)
#define LAPACKE_dposvxx   LAPACKE_NAME(dposvxx,DPOSVXX)
#define LAPACKE_cposvxx   LAPACKE_NAME(cposvxx,CPOSVXX)
#define LAPACKE_zposvxx   LAPACKE_NAME(zposvxx,ZPOSVXX)

#define LAPACKE_spotrf   LAPACKE_NAME(spotrf,SPOTRF)
#define LAPACKE_dpotrf   LAPACKE_NAME(dpotrf,DPOTRF)
#define LAPACKE_cpotrf   LAPACKE_NAME(cpotrf,CPOTRF)
#define LAPACKE_zpotrf   LAPACKE_NAME(zpotrf,ZPOTRF)

#define LAPACKE_spotri   LAPACKE_NAME(spotri,SPOTRI)
#define LAPACKE_dpotri   LAPACKE_NAME(dpotri,DPOTRI)
#define LAPACKE_cpotri   LAPACKE_NAME(cpotri,CPOTRI)
#define LAPACKE_zpotri   LAPACKE_NAME(zpotri,ZPOTRI)

#define LAPACKE_spotrs   LAPACKE_NAME(spotrs,SPOTRS)
#define LAPACKE_dpotrs   LAPACKE_NAME(dpotrs,DPOTRS)
#define LAPACKE_cpotrs   LAPACKE_NAME(cpotrs,CPOTRS)
#define LAPACKE_zpotrs   LAPACKE_NAME(zpotrs,ZPOTRS)

#define LAPACKE_sppcon   LAPACKE_NAME(sppcon,SPPCON)
#define LAPACKE_dppcon   LAPACKE_NAME(dppcon,DPPCON)
#define LAPACKE_cppcon   LAPACKE_NAME(cppcon,CPPCON)
#define LAPACKE_zppcon   LAPACKE_NAME(zppcon,ZPPCON)

#define LAPACKE_sppequ   LAPACKE_NAME(sppequ,SPPEQU)
#define LAPACKE_dppequ   LAPACKE_NAME(dppequ,DPPEQU)
#define LAPACKE_cppequ   LAPACKE_NAME(cppequ,CPPEQU)
#define LAPACKE_zppequ   LAPACKE_NAME(zppequ,ZPPEQU)

#define LAPACKE_spprfs   LAPACKE_NAME(spprfs,SPPRFS)
#define LAPACKE_dpprfs   LAPACKE_NAME(dpprfs,DPPRFS)
#define LAPACKE_cpprfs   LAPACKE_NAME(cpprfs,CPPRFS)
#define LAPACKE_zpprfs   LAPACKE_NAME(zpprfs,ZPPRFS)

#define LAPACKE_sppsv   LAPACKE_NAME(sppsv,SPPSV)
#define LAPACKE_dppsv   LAPACKE_NAME(dppsv,DPPSV)
#define LAPACKE_cppsv   LAPACKE_NAME(cppsv,CPPSV)
#define LAPACKE_zppsv   LAPACKE_NAME(zppsv,ZPPSV)

#define LAPACKE_sppsvx   LAPACKE_NAME(sppsvx,SPPSVX)
#define LAPACKE_dppsvx   LAPACKE_NAME(dppsvx,DPPSVX)
#define LAPACKE_cppsvx   LAPACKE_NAME(cppsvx,CPPSVX)
#define LAPACKE_zppsvx   LAPACKE_NAME(zppsvx,ZPPSVX)

#define LAPACKE_spptrf   LAPACKE_NAME(spptrf,SPPTRF)
#define LAPACKE_dpptrf   LAPACKE_NAME(dpptrf,DPPTRF)
#define LAPACKE_cpptrf   LAPACKE_NAME(cpptrf,CPPTRF)
#define LAPACKE_zpptrf   LAPACKE_NAME(zpptrf,ZPPTRF)

#define LAPACKE_spptri   LAPACKE_NAME(spptri,SPPTRI)
#define LAPACKE_dpptri   LAPACKE_NAME(dpptri,DPPTRI)
#define LAPACKE_cpptri   LAPACKE_NAME(cpptri,CPPTRI)
#define LAPACKE_zpptri   LAPACKE_NAME(zpptri,ZPPTRI)

#define LAPACKE_spptrs   LAPACKE_NAME(spptrs,SPPTRS)
#define LAPACKE_dpptrs   LAPACKE_NAME(dpptrs,DPPTRS)
#define LAPACKE_cpptrs   LAPACKE_NAME(cpptrs,CPPTRS)
#define LAPACKE_zpptrs   LAPACKE_NAME(zpptrs,ZPPTRS)

#define LAPACKE_spstrf   LAPACKE_NAME(spstrf,SPSTRF)
#define LAPACKE_dpstrf   LAPACKE_NAME(dpstrf,DPSTRF)
#define LAPACKE_cpstrf   LAPACKE_NAME(cpstrf,CPSTRF)
#define LAPACKE_zpstrf   LAPACKE_NAME(zpstrf,ZPSTRF)

#define LAPACKE_sptcon   LAPACKE_NAME(sptcon,SPTCON)
#define LAPACKE_dptcon   LAPACKE_NAME(dptcon,DPTCON)
#define LAPACKE_cptcon   LAPACKE_NAME(cptcon,CPTCON)
#define LAPACKE_zptcon   LAPACKE_NAME(zptcon,ZPTCON)

#define LAPACKE_spteqr   LAPACKE_NAME(spteqr,SPTEQR)
#define LAPACKE_dpteqr   LAPACKE_NAME(dpteqr,DPTEQR)
#define LAPACKE_cpteqr   LAPACKE_NAME(cpteqr,CPTEQR)
#define LAPACKE_zpteqr   LAPACKE_NAME(zpteqr,ZPTEQR)

#define LAPACKE_sptrfs   LAPACKE_NAME(sptrfs,SPTRFS)
#define LAPACKE_dptrfs   LAPACKE_NAME(dptrfs,DPTRFS)
#define LAPACKE_cptrfs   LAPACKE_NAME(cptrfs,CPTRFS)
#define LAPACKE_zptrfs   LAPACKE_NAME(zptrfs,ZPTRFS)

#define LAPACKE_sptsv   LAPACKE_NAME(sptsv,SPTSV)
#define LAPACKE_dptsv   LAPACKE_NAME(dptsv,DPTSV)
#define LAPACKE_cptsv   LAPACKE_NAME(cptsv,CPTSV)
#define LAPACKE_zptsv   LAPACKE_NAME(zptsv,ZPTSV)

#define LAPACKE_sptsvx   LAPACKE_NAME(sptsvx,SPTSVX)
#define LAPACKE_dptsvx   LAPACKE_NAME(dptsvx,DPTSVX)
#define LAPACKE_cptsvx   LAPACKE_NAME(cptsvx,CPTSVX)
#define LAPACKE_zptsvx   LAPACKE_NAME(zptsvx,ZPTSVX)

#define LAPACKE_spttrf   LAPACKE_NAME(spttrf,SPTTRF)
#define LAPACKE_dpttrf   LAPACKE_NAME(dpttrf,DPTTRF)
#define LAPACKE_cpttrf   LAPACKE_NAME(cpttrf,CPTTRF)
#define LAPACKE_zpttrf   LAPACKE_NAME(zpttrf,ZPTTRF)

#define LAPACKE_spttrs   LAPACKE_NAME(spttrs,SPTTRS)
#define LAPACKE_dpttrs   LAPACKE_NAME(dpttrs,DPTTRS)
#define LAPACKE_cpttrs   LAPACKE_NAME(cpttrs,CPTTRS)
#define LAPACKE_zpttrs   LAPACKE_NAME(zpttrs,ZPTTRS)

#define LAPACKE_ssbev   LAPACKE_NAME(ssbev,SSBEV)
#define LAPACKE_dsbev   LAPACKE_NAME(dsbev,DSBEV)

#define LAPACKE_ssbevd   LAPACKE_NAME(ssbevd,SSBEVD)
#define LAPACKE_dsbevd   LAPACKE_NAME(dsbevd,DSBEVD)

#define LAPACKE_ssbevx   LAPACKE_NAME(ssbevx,SSBEVX)
#define LAPACKE_dsbevx   LAPACKE_NAME(dsbevx,DSBEVX)

#define LAPACKE_ssbgst   LAPACKE_NAME(ssbgst,SSBGST)
#define LAPACKE_dsbgst   LAPACKE_NAME(dsbgst,DSBGST)

#define LAPACKE_ssbgv   LAPACKE_NAME(ssbgv,SSBGV)
#define LAPACKE_dsbgv   LAPACKE_NAME(dsbgv,DSBGV)

#define LAPACKE_ssbgvd   LAPACKE_NAME(ssbgvd,SSBGVD)
#define LAPACKE_dsbgvd   LAPACKE_NAME(dsbgvd,DSBGVD)

#define LAPACKE_ssbgvx   LAPACKE_NAME(ssbgvx,SSBGVX)
#define LAPACKE_dsbgvx   LAPACKE_NAME(dsbgvx,DSBGVX)

#define LAPACKE_ssbtrd   LAPACKE_NAME(ssbtrd,SSBTRD)
#define LAPACKE_dsbtrd   LAPACKE_NAME(dsbtrd,DSBTRD)

#define LAPACKE_ssfrk   LAPACKE_NAME(ssfrk,SSFRK)
#define LAPACKE_dsfrk   LAPACKE_NAME(dsfrk,DSFRK)

#define LAPACKE_sspcon   LAPACKE_NAME(sspcon,SSPCON)
#define LAPACKE_dspcon   LAPACKE_NAME(dspcon,DSPCON)
#define LAPACKE_cspcon   LAPACKE_NAME(cspcon,CSPCON)
#define LAPACKE_zspcon   LAPACKE_NAME(zspcon,ZSPCON)

#define LAPACKE_sspev   LAPACKE_NAME(sspev,SSPEV)
#define LAPACKE_dspev   LAPACKE_NAME(dspev,DSPEV)

#define LAPACKE_sspevd   LAPACKE_NAME(sspevd,SSPEVD)
#define LAPACKE_dspevd   LAPACKE_NAME(dspevd,DSPEVD)

#define LAPACKE_sspevx   LAPACKE_NAME(sspevx,SSPEVX)
#define LAPACKE_dspevx   LAPACKE_NAME(dspevx,DSPEVX)

#define LAPACKE_sspgst   LAPACKE_NAME(sspgst,SSPGST)
#define LAPACKE_dspgst   LAPACKE_NAME(dspgst,DSPGST)

#define LAPACKE_sspgv   LAPACKE_NAME(sspgv,SSPGV)
#define LAPACKE_dspgv   LAPACKE_NAME(dspgv,DSPGV)

#define LAPACKE_sspgvd   LAPACKE_NAME(sspgvd,SSPGVD)
#define LAPACKE_dspgvd   LAPACKE_NAME(dspgvd,DSPGVD)

#define LAPACKE_sspgvx   LAPACKE_NAME(sspgvx,SSPGVX)
#define LAPACKE_dspgvx   LAPACKE_NAME(dspgvx,DSPGVX)

#define LAPACKE_ssprfs   LAPACKE_NAME(ssprfs,SSPRFS)
#define LAPACKE_dsprfs   LAPACKE_NAME(dsprfs,DSPRFS)
#define LAPACKE_csprfs   LAPACKE_NAME(csprfs,CSPRFS)
#define LAPACKE_zsprfs   LAPACKE_NAME(zsprfs,ZSPRFS)

#define LAPACKE_sspsv   LAPACKE_NAME(sspsv,SSPSV)
#define LAPACKE_dspsv   LAPACKE_NAME(dspsv,DSPSV)
#define LAPACKE_cspsv   LAPACKE_NAME(cspsv,CSPSV)
#define LAPACKE_zspsv   LAPACKE_NAME(zspsv,ZSPSV)

#define LAPACKE_sspsvx   LAPACKE_NAME(sspsvx,SSPSVX)
#define LAPACKE_dspsvx   LAPACKE_NAME(dspsvx,DSPSVX)
#define LAPACKE_cspsvx   LAPACKE_NAME(cspsvx,CSPSVX)
#define LAPACKE_zspsvx   LAPACKE_NAME(zspsvx,ZSPSVX)

#define LAPACKE_ssptrd   LAPACKE_NAME(ssptrd,SSPTRD)
#define LAPACKE_dsptrd   LAPACKE_NAME(dsptrd,DSPTRD)

#define LAPACKE_ssptrf   LAPACKE_NAME(ssptrf,SSPTRF)
#define LAPACKE_dsptrf   LAPACKE_NAME(dsptrf,DSPTRF)
#define LAPACKE_csptrf   LAPACKE_NAME(csptrf,CSPTRF)
#define LAPACKE_zsptrf   LAPACKE_NAME(zsptrf,ZSPTRF)

#define LAPACKE_ssptri   LAPACKE_NAME(ssptri,SSPTRI)
#define LAPACKE_dsptri   LAPACKE_NAME(dsptri,DSPTRI)
#define LAPACKE_csptri   LAPACKE_NAME(csptri,CSPTRI)
#define LAPACKE_zsptri   LAPACKE_NAME(zsptri,ZSPTRI)

#define LAPACKE_ssptrs   LAPACKE_NAME(ssptrs,SSPTRS)
#define LAPACKE_dsptrs   LAPACKE_NAME(dsptrs,DSPTRS)
#define LAPACKE_csptrs   LAPACKE_NAME(csptrs,CSPTRS)
#define LAPACKE_zsptrs   LAPACKE_NAME(zsptrs,ZSPTRS)

#define LAPACKE_sstebz   LAPACKE_NAME(sstebz,SSTEBZ)
#define LAPACKE_dstebz   LAPACKE_NAME(dstebz,DSTEBZ)

#define LAPACKE_sstedc   LAPACKE_NAME(sstedc,SSTEDC)
#define LAPACKE_dstedc   LAPACKE_NAME(dstedc,DSTEDC)
#define LAPACKE_cstedc   LAPACKE_NAME(cstedc,CSTEDC)
#define LAPACKE_zstedc   LAPACKE_NAME(zstedc,ZSTEDC)

#define LAPACKE_sstegr   LAPACKE_NAME(sstegr,SSTEGR)
#define LAPACKE_dstegr   LAPACKE_NAME(dstegr,DSTEGR)
#define LAPACKE_cstegr   LAPACKE_NAME(cstegr,CSTEGR)
#define LAPACKE_zstegr   LAPACKE_NAME(zstegr,ZSTEGR)

#define LAPACKE_sstein   LAPACKE_NAME(sstein,SSTEIN)
#define LAPACKE_dstein   LAPACKE_NAME(dstein,DSTEIN)
#define LAPACKE_cstein   LAPACKE_NAME(cstein,CSTEIN)
#define LAPACKE_zstein   LAPACKE_NAME(zstein,ZSTEIN)

#define LAPACKE_sstemr   LAPACKE_NAME(sstemr,SSTEMR)
#define LAPACKE_dstemr   LAPACKE_NAME(dstemr,DSTEMR)
#define LAPACKE_cstemr   LAPACKE_NAME(cstemr,CSTEMR)
#define LAPACKE_zstemr   LAPACKE_NAME(zstemr,ZSTEMR)

#define LAPACKE_ssteqr   LAPACKE_NAME(ssteqr,SSTEQR)
#define LAPACKE_dsteqr   LAPACKE_NAME(dsteqr,DSTEQR)
#define LAPACKE_csteqr   LAPACKE_NAME(csteqr,CSTEQR)
#define LAPACKE_zsteqr   LAPACKE_NAME(zsteqr,ZSTEQR)

#define LAPACKE_ssterf   LAPACKE_NAME(ssterf,SSTERF)
#define LAPACKE_dsterf   LAPACKE_NAME(dsterf,DSTERF)

#define LAPACKE_sstev   LAPACKE_NAME(sstev,SSTEV)
#define LAPACKE_dstev   LAPACKE_NAME(dstev,DSTEV)

#define LAPACKE_sstevd   LAPACKE_NAME(sstevd,SSTEVD)
#define LAPACKE_dstevd   LAPACKE_NAME(dstevd,DSTEVD)

#define LAPACKE_sstevr   LAPACKE_NAME(sstevr,SSTEVR)
#define LAPACKE_dstevr   LAPACKE_NAME(dstevr,DSTEVR)

#define LAPACKE_sstevx   LAPACKE_NAME(sstevx,SSTEVX)
#define LAPACKE_dstevx   LAPACKE_NAME(dstevx,DSTEVX)

#define LAPACKE_ssycon   LAPACKE_NAME(ssycon,SSYCON)
#define LAPACKE_dsycon   LAPACKE_NAME(dsycon,DSYCON)
#define LAPACKE_csycon   LAPACKE_NAME(csycon,CSYCON)
#define LAPACKE_zsycon   LAPACKE_NAME(zsycon,ZSYCON)

#define LAPACKE_ssyequb   LAPACKE_NAME(ssyequb,SSYEQUB)
#define LAPACKE_dsyequb   LAPACKE_NAME(dsyequb,DSYEQUB)
#define LAPACKE_csyequb   LAPACKE_NAME(csyequb,CSYEQUB)
#define LAPACKE_zsyequb   LAPACKE_NAME(zsyequb,ZSYEQUB)

#define LAPACKE_ssyev   LAPACKE_NAME(ssyev,SSYEV)
#define LAPACKE_dsyev   LAPACKE_NAME(dsyev,DSYEV)

#define LAPACKE_ssyevd   LAPACKE_NAME(ssyevd,SSYEVD)
#define LAPACKE_dsyevd   LAPACKE_NAME(dsyevd,DSYEVD)

#define LAPACKE_ssyevr   LAPACKE_NAME(ssyevr,SSYEVR)
#define LAPACKE_dsyevr   LAPACKE_NAME(dsyevr,DSYEVR)

#define LAPACKE_ssyevx   LAPACKE_NAME(ssyevx,SSYEVX)
#define LAPACKE_dsyevx   LAPACKE_NAME(dsyevx,DSYEVX)

#define LAPACKE_ssygst   LAPACKE_NAME(ssygst,SSYGST)
#define LAPACKE_dsygst   LAPACKE_NAME(dsygst,DSYGST)

#define LAPACKE_ssygv   LAPACKE_NAME(ssygv,SSYGV)
#define LAPACKE_dsygv   LAPACKE_NAME(dsygv,DSYGV)

#define LAPACKE_ssygvd   LAPACKE_NAME(ssygvd,SSYGVD)
#define LAPACKE_dsygvd   LAPACKE_NAME(dsygvd,DSYGVD)

#define LAPACKE_ssygvx   LAPACKE_NAME(ssygvx,SSYGVX)
#define LAPACKE_dsygvx   LAPACKE_NAME(dsygvx,DSYGVX)

#define LAPACKE_ssyrfs   LAPACKE_NAME(ssyrfs,SSYRFS)
#define LAPACKE_dsyrfs   LAPACKE_NAME(dsyrfs,DSYRFS)
#define LAPACKE_csyrfs   LAPACKE_NAME(csyrfs,CSYRFS)
#define LAPACKE_zsyrfs   LAPACKE_NAME(zsyrfs,ZSYRFS)

#define LAPACKE_ssyrfsx   LAPACKE_NAME(ssyrfsx,SSYRFSX)
#define LAPACKE_dsyrfsx   LAPACKE_NAME(dsyrfsx,DSYRFSX)
#define LAPACKE_csyrfsx   LAPACKE_NAME(csyrfsx,CSYRFSX)
#define LAPACKE_zsyrfsx   LAPACKE_NAME(zsyrfsx,ZSYRFSX)

#define LAPACKE_ssysv   LAPACKE_NAME(ssysv,SSYSV)
#define LAPACKE_dsysv   LAPACKE_NAME(dsysv,DSYSV)
#define LAPACKE_csysv   LAPACKE_NAME(csysv,CSYSV)
#define LAPACKE_zsysv   LAPACKE_NAME(zsysv,ZSYSV)

#define LAPACKE_ssysvx   LAPACKE_NAME(ssysvx,SSYSVX)
#define LAPACKE_dsysvx   LAPACKE_NAME(dsysvx,DSYSVX)
#define LAPACKE_csysvx   LAPACKE_NAME(csysvx,CSYSVX)
#define LAPACKE_zsysvx   LAPACKE_NAME(zsysvx,ZSYSVX)

#define LAPACKE_ssysvxx   LAPACKE_NAME(ssysvxx,SSYSVXX)
#define LAPACKE_dsysvxx   LAPACKE_NAME(dsysvxx,DSYSVXX)
#define LAPACKE_csysvxx   LAPACKE_NAME(csysvxx,CSYSVXX)
#define LAPACKE_zsysvxx   LAPACKE_NAME(zsysvxx,ZSYSVXX)

#define LAPACKE_ssytrd   LAPACKE_NAME(ssytrd,SSYTRD)
#define LAPACKE_dsytrd   LAPACKE_NAME(dsytrd,DSYTRD)

#define LAPACKE_ssytrf   LAPACKE_NAME(ssytrf,SSYTRF)
#define LAPACKE_dsytrf   LAPACKE_NAME(dsytrf,DSYTRF)
#define LAPACKE_csytrf   LAPACKE_NAME(csytrf,CSYTRF)
#define LAPACKE_zsytrf   LAPACKE_NAME(zsytrf,ZSYTRF)

#define LAPACKE_ssytri   LAPACKE_NAME(ssytri,SSYTRI)
#define LAPACKE_dsytri   LAPACKE_NAME(dsytri,DSYTRI)
#define LAPACKE_csytri   LAPACKE_NAME(csytri,CSYTRI)
#define LAPACKE_zsytri   LAPACKE_NAME(zsytri,ZSYTRI)

#define LAPACKE_ssytrs   LAPACKE_NAME(ssytrs,SSYTRS)
#define LAPACKE_dsytrs   LAPACKE_NAME(dsytrs,DSYTRS)
#define LAPACKE_csytrs   LAPACKE_NAME(csytrs,CSYTRS)
#define LAPACKE_zsytrs   LAPACKE_NAME(zsytrs,ZSYTRS)

#define LAPACKE_stbcon   LAPACKE_NAME(stbcon,STBCON)
#define LAPACKE_dtbcon   LAPACKE_NAME(dtbcon,DTBCON)
#define LAPACKE_ctbcon   LAPACKE_NAME(ctbcon,CTBCON)
#define LAPACKE_ztbcon   LAPACKE_NAME(ztbcon,ZTBCON)

#define LAPACKE_stbrfs   LAPACKE_NAME(stbrfs,STBRFS)
#define LAPACKE_dtbrfs   LAPACKE_NAME(dtbrfs,DTBRFS)
#define LAPACKE_ctbrfs   LAPACKE_NAME(ctbrfs,CTBRFS)
#define LAPACKE_ztbrfs   LAPACKE_NAME(ztbrfs,ZTBRFS)

#define LAPACKE_stbtrs   LAPACKE_NAME(stbtrs,STBTRS)
#define LAPACKE_dtbtrs   LAPACKE_NAME(dtbtrs,DTBTRS)
#define LAPACKE_ctbtrs   LAPACKE_NAME(ctbtrs,CTBTRS)
#define LAPACKE_ztbtrs   LAPACKE_NAME(ztbtrs,ZTBTRS)

#define LAPACKE_stfsm   LAPACKE_NAME(stfsm,STFSM)
#define LAPACKE_dtfsm   LAPACKE_NAME(dtfsm,DTFSM)
#define LAPACKE_ctfsm   LAPACKE_NAME(ctfsm,CTFSM)
#define LAPACKE_ztfsm   LAPACKE_NAME(ztfsm,ZTFSM)

#define LAPACKE_stftri   LAPACKE_NAME(stftri,STFTRI)
#define LAPACKE_dtftri   LAPACKE_NAME(dtftri,DTFTRI)
#define LAPACKE_ctftri   LAPACKE_NAME(ctftri,CTFTRI)
#define LAPACKE_ztftri   LAPACKE_NAME(ztftri,ZTFTRI)

#define LAPACKE_stfttp   LAPACKE_NAME(stfttp,STFTTP)
#define LAPACKE_dtfttp   LAPACKE_NAME(dtfttp,DTFTTP)
#define LAPACKE_ctfttp   LAPACKE_NAME(ctfttp,CTFTTP)
#define LAPACKE_ztfttp   LAPACKE_NAME(ztfttp,ZTFTTP)

#define LAPACKE_stfttr   LAPACKE_NAME(stfttr,STFTTR)
#define LAPACKE_dtfttr   LAPACKE_NAME(dtfttr,DTFTTR)
#define LAPACKE_ctfttr   LAPACKE_NAME(ctfttr,CTFTTR)
#define LAPACKE_ztfttr   LAPACKE_NAME(ztfttr,ZTFTTR)

#define LAPACKE_stgevc   LAPACKE_NAME(stgevc,STGEVC)
#define LAPACKE_dtgevc   LAPACKE_NAME(dtgevc,DTGEVC)
#define LAPACKE_ctgevc   LAPACKE_NAME(ctgevc,CTGEVC)
#define LAPACKE_ztgevc   LAPACKE_NAME(ztgevc,ZTGEVC)

#define LAPACKE_stgexc   LAPACKE_NAME(stgexc,STGEXC)
#define LAPACKE_dtgexc   LAPACKE_NAME(dtgexc,DTGEXC)
#define LAPACKE_ctgexc   LAPACKE_NAME(ctgexc,CTGEXC)
#define LAPACKE_ztgexc   LAPACKE_NAME(ztgexc,ZTGEXC)

#define LAPACKE_stgsen   LAPACKE_NAME(stgsen,STGSEN)
#define LAPACKE_dtgsen   LAPACKE_NAME(dtgsen,DTGSEN)
#define LAPACKE_ctgsen   LAPACKE_NAME(ctgsen,CTGSEN)
#define LAPACKE_ztgsen   LAPACKE_NAME(ztgsen,ZTGSEN)

#define LAPACKE_stgsja   LAPACKE_NAME(stgsja,STGSJA)
#define LAPACKE_dtgsja   LAPACKE_NAME(dtgsja,DTGSJA)
#define LAPACKE_ctgsja   LAPACKE_NAME(ctgsja,CTGSJA)
#define LAPACKE_ztgsja   LAPACKE_NAME(ztgsja,ZTGSJA)

#define LAPACKE_stgsna   LAPACKE_NAME(stgsna,STGSNA)
#define LAPACKE_dtgsna   LAPACKE_NAME(dtgsna,DTGSNA)
#define LAPACKE_ctgsna   LAPACKE_NAME(ctgsna,CTGSNA)
#define LAPACKE_ztgsna   LAPACKE_NAME(ztgsna,ZTGSNA)

#define LAPACKE_stgsyl   LAPACKE_NAME(stgsyl,STGSYL)
#define LAPACKE_dtgsyl   LAPACKE_NAME(dtgsyl,DTGSYL)
#define LAPACKE_ctgsyl   LAPACKE_NAME(ctgsyl,CTGSYL)
#define LAPACKE_ztgsyl   LAPACKE_NAME(ztgsyl,ZTGSYL)

#define LAPACKE_stpcon   LAPACKE_NAME(stpcon,STPCON)
#define LAPACKE_dtpcon   LAPACKE_NAME(dtpcon,DTPCON)
#define LAPACKE_ctpcon   LAPACKE_NAME(ctpcon,CTPCON)
#define LAPACKE_ztpcon   LAPACKE_NAME(ztpcon,ZTPCON)

#define LAPACKE_stprfs   LAPACKE_NAME(stprfs,STPRFS)
#define LAPACKE_dtprfs   LAPACKE_NAME(dtprfs,DTPRFS)
#define LAPACKE_ctprfs   LAPACKE_NAME(ctprfs,CTPRFS)
#define LAPACKE_ztprfs   LAPACKE_NAME(ztprfs,ZTPRFS)

#define LAPACKE_stptri   LAPACKE_NAME(stptri,STPTRI)
#define LAPACKE_dtptri   LAPACKE_NAME(dtptri,DTPTRI)
#define LAPACKE_ctptri   LAPACKE_NAME(ctptri,CTPTRI)
#define LAPACKE_ztptri   LAPACKE_NAME(ztptri,ZTPTRI)

#define LAPACKE_stptrs   LAPACKE_NAME(stptrs,STPTRS)
#define LAPACKE_dtptrs   LAPACKE_NAME(dtptrs,DTPTRS)
#define LAPACKE_ctptrs   LAPACKE_NAME(ctptrs,CTPTRS)
#define LAPACKE_ztptrs   LAPACKE_NAME(ztptrs,ZTPTRS)

#define LAPACKE_stpttf   LAPACKE_NAME(stpttf,STPTTF)
#define LAPACKE_dtpttf   LAPACKE_NAME(dtpttf,DTPTTF)
#define LAPACKE_ctpttf   LAPACKE_NAME(ctpttf,CTPTTF)
#define LAPACKE_ztpttf   LAPACKE_NAME(ztpttf,ZTPTTF)

#define LAPACKE_stpttr   LAPACKE_NAME(stpttr,STPTTR)
#define LAPACKE_dtpttr   LAPACKE_NAME(dtpttr,DTPTTR)
#define LAPACKE_ctpttr   LAPACKE_NAME(ctpttr,CTPTTR)
#define LAPACKE_ztpttr   LAPACKE_NAME(ztpttr,ZTPTTR)

#define LAPACKE_strcon   LAPACKE_NAME(strcon,STRCON)
#define LAPACKE_dtrcon   LAPACKE_NAME(dtrcon,DTRCON)
#define LAPACKE_ctrcon   LAPACKE_NAME(ctrcon,CTRCON)
#define LAPACKE_ztrcon   LAPACKE_NAME(ztrcon,ZTRCON)

#define LAPACKE_strevc   LAPACKE_NAME(strevc,STREVC)
#define LAPACKE_dtrevc   LAPACKE_NAME(dtrevc,DTREVC)
#define LAPACKE_ctrevc   LAPACKE_NAME(ctrevc,CTREVC)
#define LAPACKE_ztrevc   LAPACKE_NAME(ztrevc,ZTREVC)

#define LAPACKE_strexc   LAPACKE_NAME(strexc,STREXC)
#define LAPACKE_dtrexc   LAPACKE_NAME(dtrexc,DTREXC)
#define LAPACKE_ctrexc   LAPACKE_NAME(ctrexc,CTREXC)
#define LAPACKE_ztrexc   LAPACKE_NAME(ztrexc,ZTREXC)

#define LAPACKE_strrfs   LAPACKE_NAME(strrfs,STRRFS)
#define LAPACKE_dtrrfs   LAPACKE_NAME(dtrrfs,DTRRFS)
#define LAPACKE_ctrrfs   LAPACKE_NAME(ctrrfs,CTRRFS)
#define LAPACKE_ztrrfs   LAPACKE_NAME(ztrrfs,ZTRRFS)

#define LAPACKE_strsen   LAPACKE_NAME(strsen,STRSEN)
#define LAPACKE_dtrsen   LAPACKE_NAME(dtrsen,DTRSEN)
#define LAPACKE_ctrsen   LAPACKE_NAME(ctrsen,CTRSEN)
#define LAPACKE_ztrsen   LAPACKE_NAME(ztrsen,ZTRSEN)

#define LAPACKE_strsna   LAPACKE_NAME(strsna,STRSNA)
#define LAPACKE_dtrsna   LAPACKE_NAME(dtrsna,DTRSNA)
#define LAPACKE_ctrsna   LAPACKE_NAME(ctrsna,CTRSNA)
#define LAPACKE_ztrsna   LAPACKE_NAME(ztrsna,ZTRSNA)

#define LAPACKE_strsyl   LAPACKE_NAME(strsyl,STRSYL)
#define LAPACKE_dtrsyl   LAPACKE_NAME(dtrsyl,DTRSYL)
#define LAPACKE_ctrsyl   LAPACKE_NAME(ctrsyl,CTRSYL)
#define LAPACKE_ztrsyl   LAPACKE_NAME(ztrsyl,ZTRSYL)

#define LAPACKE_strtri   LAPACKE_NAME(strtri,STRTRI)
#define LAPACKE_dtrtri   LAPACKE_NAME(dtrtri,DTRTRI)
#define LAPACKE_ctrtri   LAPACKE_NAME(ctrtri,CTRTRI)
#define LAPACKE_ztrtri   LAPACKE_NAME(ztrtri,ZTRTRI)

#define LAPACKE_strtrs   LAPACKE_NAME(strtrs,STRTRS)
#define LAPACKE_dtrtrs   LAPACKE_NAME(dtrtrs,DTRTRS)
#define LAPACKE_ctrtrs   LAPACKE_NAME(ctrtrs,CTRTRS)
#define LAPACKE_ztrtrs   LAPACKE_NAME(ztrtrs,ZTRTRS)

#define LAPACKE_strttf   LAPACKE_NAME(strttf,STRTTF)
#define LAPACKE_dtrttf   LAPACKE_NAME(dtrttf,DTRTTF)
#define LAPACKE_ctrttf   LAPACKE_NAME(ctrttf,CTRTTF)
#define LAPACKE_ztrttf   LAPACKE_NAME(ztrttf,ZTRTTF)

#define LAPACKE_strttp   LAPACKE_NAME(strttp,STRTTP)
#define LAPACKE_dtrttp   LAPACKE_NAME(dtrttp,DTRTTP)
#define LAPACKE_ctrttp   LAPACKE_NAME(ctrttp,CTRTTP)
#define LAPACKE_ztrttp   LAPACKE_NAME(ztrttp,ZTRTTP)

#define LAPACKE_stzrzf   LAPACKE_NAME(stzrzf,STZRZF)
#define LAPACKE_dtzrzf   LAPACKE_NAME(dtzrzf,DTZRZF)
#define LAPACKE_ctzrzf   LAPACKE_NAME(ctzrzf,CTZRZF)
#define LAPACKE_ztzrzf   LAPACKE_NAME(ztzrzf,ZTZRZF)

#define LAPACKE_cungbr   LAPACKE_NAME(cungbr,CUNGBR)
#define LAPACKE_zungbr   LAPACKE_NAME(zungbr,ZUNGBR)

#define LAPACKE_cunghr   LAPACKE_NAME(cunghr,CUNGHR)
#define LAPACKE_zunghr   LAPACKE_NAME(zunghr,ZUNGHR)

#define LAPACKE_cunglq   LAPACKE_NAME(cunglq,CUNGLQ)
#define LAPACKE_zunglq   LAPACKE_NAME(zunglq,ZUNGLQ)

#define LAPACKE_cungql   LAPACKE_NAME(cungql,CUNGQL)
#define LAPACKE_zungql   LAPACKE_NAME(zungql,ZUNGQL)

#define LAPACKE_cungqr   LAPACKE_NAME(cungqr,CUNGQR)
#define LAPACKE_zungqr   LAPACKE_NAME(zungqr,ZUNGQR)

#define LAPACKE_cungrq   LAPACKE_NAME(cungrq,CUNGRQ)
#define LAPACKE_zungrq   LAPACKE_NAME(zungrq,ZUNGRQ)

#define LAPACKE_cungtr   LAPACKE_NAME(cungtr,CUNGTR)
#define LAPACKE_zungtr   LAPACKE_NAME(zungtr,ZUNGTR)

#define LAPACKE_cunmbr   LAPACKE_NAME(cunmbr,CUNMBR)
#define LAPACKE_zunmbr   LAPACKE_NAME(zunmbr,ZUNMBR)

#define LAPACKE_cunmhr   LAPACKE_NAME(cunmhr,CUNMHR)
#define LAPACKE_zunmhr   LAPACKE_NAME(zunmhr,ZUNMHR)

#define LAPACKE_cunmlq   LAPACKE_NAME(cunmlq,CUNMLQ)
#define LAPACKE_zunmlq   LAPACKE_NAME(zunmlq,ZUNMLQ)

#define LAPACKE_cunmql   LAPACKE_NAME(cunmql,CUNMQL)
#define LAPACKE_zunmql   LAPACKE_NAME(zunmql,ZUNMQL)

#define LAPACKE_cunmqr   LAPACKE_NAME(cunmqr,CUNMQR)
#define LAPACKE_zunmqr   LAPACKE_NAME(zunmqr,ZUNMQR)

#define LAPACKE_cunmrq   LAPACKE_NAME(cunmrq,CUNMRQ)
#define LAPACKE_zunmrq   LAPACKE_NAME(zunmrq,ZUNMRQ)

#define LAPACKE_cunmrz   LAPACKE_NAME(cunmrz,CUNMRZ)
#define LAPACKE_zunmrz   LAPACKE_NAME(zunmrz,ZUNMRZ)

#define LAPACKE_cunmtr   LAPACKE_NAME(cunmtr,CUNMTR)
#define LAPACKE_zunmtr   LAPACKE_NAME(zunmtr,ZUNMTR)

#define LAPACKE_cupgtr   LAPACKE_NAME(cupgtr,CUPGTR)
#define LAPACKE_zupgtr   LAPACKE_NAME(zupgtr,ZUPGTR)

#define LAPACKE_cupmtr   LAPACKE_NAME(cupmtr,CUPMTR)
#define LAPACKE_zupmtr   LAPACKE_NAME(zupmtr,ZUPMTR)

#define LAPACKE_sbdsdc_work   LAPACKE_NAME(sbdsdc_work,SBDSDC_WORK)
#define LAPACKE_dbdsdc_work   LAPACKE_NAME(dbdsdc_work,DBDSDC_WORK)

#define LAPACKE_sbdsqr_work   LAPACKE_NAME(sbdsqr_work,SBDSQR_WORK)
#define LAPACKE_dbdsqr_work   LAPACKE_NAME(dbdsqr_work,DBDSQR_WORK)
#define LAPACKE_cbdsqr_work   LAPACKE_NAME(cbdsqr_work,CBDSQR_WORK)
#define LAPACKE_zbdsqr_work   LAPACKE_NAME(zbdsqr_work,ZBDSQR_WORK)

#define LAPACKE_sdisna_work   LAPACKE_NAME(sdisna_work,SDISNA_WORK)
#define LAPACKE_ddisna_work   LAPACKE_NAME(ddisna_work,DDISNA_WORK)

#define LAPACKE_sgbbrd_work   LAPACKE_NAME(sgbbrd_work,SGBBRD_WORK)
#define LAPACKE_dgbbrd_work   LAPACKE_NAME(dgbbrd_work,DGBBRD_WORK)
#define LAPACKE_cgbbrd_work   LAPACKE_NAME(cgbbrd_work,CGBBRD_WORK)
#define LAPACKE_zgbbrd_work   LAPACKE_NAME(zgbbrd_work,ZGBBRD_WORK)

#define LAPACKE_sgbcon_work   LAPACKE_NAME(sgbcon_work,SGBCON_WORK)
#define LAPACKE_dgbcon_work   LAPACKE_NAME(dgbcon_work,DGBCON_WORK)
#define LAPACKE_cgbcon_work   LAPACKE_NAME(cgbcon_work,CGBCON_WORK)
#define LAPACKE_zgbcon_work   LAPACKE_NAME(zgbcon_work,ZGBCON_WORK)

#define LAPACKE_sgbequ_work   LAPACKE_NAME(sgbequ_work,SGBEQU_WORK)
#define LAPACKE_dgbequ_work   LAPACKE_NAME(dgbequ_work,DGBEQU_WORK)
#define LAPACKE_cgbequ_work   LAPACKE_NAME(cgbequ_work,CGBEQU_WORK)
#define LAPACKE_zgbequ_work   LAPACKE_NAME(zgbequ_work,ZGBEQU_WORK)

#define LAPACKE_sgbequb_work   LAPACKE_NAME(sgbequb_work,SGBEQUB_WORK)
#define LAPACKE_dgbequb_work   LAPACKE_NAME(dgbequb_work,DGBEQUB_WORK)
#define LAPACKE_cgbequb_work   LAPACKE_NAME(cgbequb_work,CGBEQUB_WORK)
#define LAPACKE_zgbequb_work   LAPACKE_NAME(zgbequb_work,ZGBEQUB_WORK)

#define LAPACKE_sgbrfs_work   LAPACKE_NAME(sgbrfs_work,SGBRFS_WORK)
#define LAPACKE_dgbrfs_work   LAPACKE_NAME(dgbrfs_work,DGBRFS_WORK)
#define LAPACKE_cgbrfs_work   LAPACKE_NAME(cgbrfs_work,CGBRFS_WORK)
#define LAPACKE_zgbrfs_work   LAPACKE_NAME(zgbrfs_work,ZGBRFS_WORK)

#define LAPACKE_sgbrfsx_work   LAPACKE_NAME(sgbrfsx_work,SGBRFSX_WORK)
#define LAPACKE_dgbrfsx_work   LAPACKE_NAME(dgbrfsx_work,DGBRFSX_WORK)
#define LAPACKE_cgbrfsx_work   LAPACKE_NAME(cgbrfsx_work,CGBRFSX_WORK)
#define LAPACKE_zgbrfsx_work   LAPACKE_NAME(zgbrfsx_work,ZGBRFSX_WORK)

#define LAPACKE_sgbsv_work   LAPACKE_NAME(sgbsv_work,SGBSV_WORK)
#define LAPACKE_dgbsv_work   LAPACKE_NAME(dgbsv_work,DGBSV_WORK)
#define LAPACKE_cgbsv_work   LAPACKE_NAME(cgbsv_work,CGBSV_WORK)
#define LAPACKE_zgbsv_work   LAPACKE_NAME(zgbsv_work,ZGBSV_WORK)

#define LAPACKE_sgbsvx_work   LAPACKE_NAME(sgbsvx_work,SGBSVX_WORK)
#define LAPACKE_dgbsvx_work   LAPACKE_NAME(dgbsvx_work,DGBSVX_WORK)
#define LAPACKE_cgbsvx_work   LAPACKE_NAME(cgbsvx_work,CGBSVX_WORK)
#define LAPACKE_zgbsvx_work   LAPACKE_NAME(zgbsvx_work,ZGBSVX_WORK)

#define LAPACKE_sgbsvxx_work   LAPACKE_NAME(sgbsvxx_work,SGBSVXX_WORK)
#define LAPACKE_dgbsvxx_work   LAPACKE_NAME(dgbsvxx_work,DGBSVXX_WORK)
#define LAPACKE_cgbsvxx_work   LAPACKE_NAME(cgbsvxx_work,CGBSVXX_WORK)
#define LAPACKE_zgbsvxx_work   LAPACKE_NAME(zgbsvxx_work,ZGBSVXX_WORK)

#define LAPACKE_sgbtrf_work   LAPACKE_NAME(sgbtrf_work,SGBTRF_WORK)
#define LAPACKE_dgbtrf_work   LAPACKE_NAME(dgbtrf_work,DGBTRF_WORK)
#define LAPACKE_cgbtrf_work   LAPACKE_NAME(cgbtrf_work,CGBTRF_WORK)
#define LAPACKE_zgbtrf_work   LAPACKE_NAME(zgbtrf_work,ZGBTRF_WORK)

#define LAPACKE_sgbtrs_work   LAPACKE_NAME(sgbtrs_work,SGBTRS_WORK)
#define LAPACKE_dgbtrs_work   LAPACKE_NAME(dgbtrs_work,DGBTRS_WORK)
#define LAPACKE_cgbtrs_work   LAPACKE_NAME(cgbtrs_work,CGBTRS_WORK)
#define LAPACKE_zgbtrs_work   LAPACKE_NAME(zgbtrs_work,ZGBTRS_WORK)

#define LAPACKE_sgebak_work   LAPACKE_NAME(sgebak_work,SGEBAK_WORK)
#define LAPACKE_dgebak_work   LAPACKE_NAME(dgebak_work,DGEBAK_WORK)
#define LAPACKE_cgebak_work   LAPACKE_NAME(cgebak_work,CGEBAK_WORK)
#define LAPACKE_zgebak_work   LAPACKE_NAME(zgebak_work,ZGEBAK_WORK)

#define LAPACKE_sgebal_work   LAPACKE_NAME(sgebal_work,SGEBAL_WORK)
#define LAPACKE_dgebal_work   LAPACKE_NAME(dgebal_work,DGEBAL_WORK)
#define LAPACKE_cgebal_work   LAPACKE_NAME(cgebal_work,CGEBAL_WORK)
#define LAPACKE_zgebal_work   LAPACKE_NAME(zgebal_work,ZGEBAL_WORK)

#define LAPACKE_sgebrd_work   LAPACKE_NAME(sgebrd_work,SGEBRD_WORK)
#define LAPACKE_dgebrd_work   LAPACKE_NAME(dgebrd_work,DGEBRD_WORK)
#define LAPACKE_cgebrd_work   LAPACKE_NAME(cgebrd_work,CGEBRD_WORK)
#define LAPACKE_zgebrd_work   LAPACKE_NAME(zgebrd_work,ZGEBRD_WORK)

#define LAPACKE_sgecon_work   LAPACKE_NAME(sgecon_work,SGECON_WORK)
#define LAPACKE_dgecon_work   LAPACKE_NAME(dgecon_work,DGECON_WORK)
#define LAPACKE_cgecon_work   LAPACKE_NAME(cgecon_work,CGECON_WORK)
#define LAPACKE_zgecon_work   LAPACKE_NAME(zgecon_work,ZGECON_WORK)

#define LAPACKE_sgeequ_work   LAPACKE_NAME(sgeequ_work,SGEEQU_WORK)
#define LAPACKE_dgeequ_work   LAPACKE_NAME(dgeequ_work,DGEEQU_WORK)
#define LAPACKE_cgeequ_work   LAPACKE_NAME(cgeequ_work,CGEEQU_WORK)
#define LAPACKE_zgeequ_work   LAPACKE_NAME(zgeequ_work,ZGEEQU_WORK)

#define LAPACKE_sgeequb_work   LAPACKE_NAME(sgeequb_work,SGEEQUB_WORK)
#define LAPACKE_dgeequb_work   LAPACKE_NAME(dgeequb_work,DGEEQUB_WORK)
#define LAPACKE_cgeequb_work   LAPACKE_NAME(cgeequb_work,CGEEQUB_WORK)
#define LAPACKE_zgeequb_work   LAPACKE_NAME(zgeequb_work,ZGEEQUB_WORK)

#define LAPACKE_sgees_work   LAPACKE_NAME(sgees_work,SGEES_WORK)
#define LAPACKE_dgees_work   LAPACKE_NAME(dgees_work,DGEES_WORK)
#define LAPACKE_cgees_work   LAPACKE_NAME(cgees_work,CGEES_WORK)
#define LAPACKE_zgees_work   LAPACKE_NAME(zgees_work,ZGEES_WORK)

#define LAPACKE_sgeesx_work   LAPACKE_NAME(sgeesx_work,SGEESX_WORK)
#define LAPACKE_dgeesx_work   LAPACKE_NAME(dgeesx_work,DGEESX_WORK)
#define LAPACKE_cgeesx_work   LAPACKE_NAME(cgeesx_work,CGEESX_WORK)
#define LAPACKE_zgeesx_work   LAPACKE_NAME(zgeesx_work,ZGEESX_WORK)

#define LAPACKE_sgeev_work   LAPACKE_NAME(sgeev_work,SGEEV_WORK)
#define LAPACKE_dgeev_work   LAPACKE_NAME(dgeev_work,DGEEV_WORK)
#define LAPACKE_cgeev_work   LAPACKE_NAME(cgeev_work,CGEEV_WORK)
#define LAPACKE_zgeev_work   LAPACKE_NAME(zgeev_work,ZGEEV_WORK)

#define LAPACKE_sgeevx_work   LAPACKE_NAME(sgeevx_work,SGEEVX_WORK)
#define LAPACKE_dgeevx_work   LAPACKE_NAME(dgeevx_work,DGEEVX_WORK)
#define LAPACKE_cgeevx_work   LAPACKE_NAME(cgeevx_work,CGEEVX_WORK)
#define LAPACKE_zgeevx_work   LAPACKE_NAME(zgeevx_work,ZGEEVX_WORK)

#define LAPACKE_sgehrd_work   LAPACKE_NAME(sgehrd_work,SGEHRD_WORK)
#define LAPACKE_dgehrd_work   LAPACKE_NAME(dgehrd_work,DGEHRD_WORK)
#define LAPACKE_cgehrd_work   LAPACKE_NAME(cgehrd_work,CGEHRD_WORK)
#define LAPACKE_zgehrd_work   LAPACKE_NAME(zgehrd_work,ZGEHRD_WORK)

#define LAPACKE_sgejsv_work   LAPACKE_NAME(sgejsv_work,SGEJSV_WORK)
#define LAPACKE_dgejsv_work   LAPACKE_NAME(dgejsv_work,DGEJSV_WORK)

#define LAPACKE_sgelq2_work   LAPACKE_NAME(sgelq2_work,SGELQ2_WORK)
#define LAPACKE_dgelq2_work   LAPACKE_NAME(dgelq2_work,DGELQ2_WORK)
#define LAPACKE_cgelq2_work   LAPACKE_NAME(cgelq2_work,CGELQ2_WORK)
#define LAPACKE_zgelq2_work   LAPACKE_NAME(zgelq2_work,ZGELQ2_WORK)

#define LAPACKE_sgelqf_work   LAPACKE_NAME(sgelqf_work,SGELQF_WORK)
#define LAPACKE_dgelqf_work   LAPACKE_NAME(dgelqf_work,DGELQF_WORK)
#define LAPACKE_cgelqf_work   LAPACKE_NAME(cgelqf_work,CGELQF_WORK)
#define LAPACKE_zgelqf_work   LAPACKE_NAME(zgelqf_work,ZGELQF_WORK)

#define LAPACKE_sgels_work   LAPACKE_NAME(sgels_work,SGELS_WORK)
#define LAPACKE_dgels_work   LAPACKE_NAME(dgels_work,DGELS_WORK)
#define LAPACKE_cgels_work   LAPACKE_NAME(cgels_work,CGELS_WORK)
#define LAPACKE_zgels_work   LAPACKE_NAME(zgels_work,ZGELS_WORK)

#define LAPACKE_sgelsd_work   LAPACKE_NAME(sgelsd_work,SGELSD_WORK)
#define LAPACKE_dgelsd_work   LAPACKE_NAME(dgelsd_work,DGELSD_WORK)
#define LAPACKE_cgelsd_work   LAPACKE_NAME(cgelsd_work,CGELSD_WORK)
#define LAPACKE_zgelsd_work   LAPACKE_NAME(zgelsd_work,ZGELSD_WORK)

#define LAPACKE_sgelss_work   LAPACKE_NAME(sgelss_work,SGELSS_WORK)
#define LAPACKE_dgelss_work   LAPACKE_NAME(dgelss_work,DGELSS_WORK)
#define LAPACKE_cgelss_work   LAPACKE_NAME(cgelss_work,CGELSS_WORK)
#define LAPACKE_zgelss_work   LAPACKE_NAME(zgelss_work,ZGELSS_WORK)

#define LAPACKE_sgelsy_work   LAPACKE_NAME(sgelsy_work,SGELSY_WORK)
#define LAPACKE_dgelsy_work   LAPACKE_NAME(dgelsy_work,DGELSY_WORK)
#define LAPACKE_cgelsy_work   LAPACKE_NAME(cgelsy_work,CGELSY_WORK)
#define LAPACKE_zgelsy_work   LAPACKE_NAME(zgelsy_work,ZGELSY_WORK)

#define LAPACKE_sgeqlf_work   LAPACKE_NAME(sgeqlf_work,SGEQLF_WORK)
#define LAPACKE_dgeqlf_work   LAPACKE_NAME(dgeqlf_work,DGEQLF_WORK)
#define LAPACKE_cgeqlf_work   LAPACKE_NAME(cgeqlf_work,CGEQLF_WORK)
#define LAPACKE_zgeqlf_work   LAPACKE_NAME(zgeqlf_work,ZGEQLF_WORK)

#define LAPACKE_sgeqp3_work   LAPACKE_NAME(sgeqp3_work,SGEQP3_WORK)
#define LAPACKE_dgeqp3_work   LAPACKE_NAME(dgeqp3_work,DGEQP3_WORK)
#define LAPACKE_cgeqp3_work   LAPACKE_NAME(cgeqp3_work,CGEQP3_WORK)
#define LAPACKE_zgeqp3_work   LAPACKE_NAME(zgeqp3_work,ZGEQP3_WORK)

#define LAPACKE_sgeqpf_work   LAPACKE_NAME(sgeqpf_work,SGEQPF_WORK)
#define LAPACKE_dgeqpf_work   LAPACKE_NAME(dgeqpf_work,DGEQPF_WORK)
#define LAPACKE_cgeqpf_work   LAPACKE_NAME(cgeqpf_work,CGEQPF_WORK)
#define LAPACKE_zgeqpf_work   LAPACKE_NAME(zgeqpf_work,ZGEQPF_WORK)

#define LAPACKE_sgeqr2_work   LAPACKE_NAME(sgeqr2_work,SGEQR2_WORK)
#define LAPACKE_dgeqr2_work   LAPACKE_NAME(dgeqr2_work,DGEQR2_WORK)
#define LAPACKE_cgeqr2_work   LAPACKE_NAME(cgeqr2_work,CGEQR2_WORK)
#define LAPACKE_zgeqr2_work   LAPACKE_NAME(zgeqr2_work,ZGEQR2_WORK)

#define LAPACKE_sgeqrf_work   LAPACKE_NAME(sgeqrf_work,SGEQRF_WORK)
#define LAPACKE_dgeqrf_work   LAPACKE_NAME(dgeqrf_work,DGEQRF_WORK)
#define LAPACKE_cgeqrf_work   LAPACKE_NAME(cgeqrf_work,CGEQRF_WORK)
#define LAPACKE_zgeqrf_work   LAPACKE_NAME(zgeqrf_work,ZGEQRF_WORK)

#define LAPACKE_sgeqrfp_work   LAPACKE_NAME(sgeqrfp_work,SGEQRFP_WORK)
#define LAPACKE_dgeqrfp_work   LAPACKE_NAME(dgeqrfp_work,DGEQRFP_WORK)
#define LAPACKE_cgeqrfp_work   LAPACKE_NAME(cgeqrfp_work,CGEQRFP_WORK)
#define LAPACKE_zgeqrfp_work   LAPACKE_NAME(zgeqrfp_work,ZGEQRFP_WORK)

#define LAPACKE_sgerfs_work   LAPACKE_NAME(sgerfs_work,SGERFS_WORK)
#define LAPACKE_dgerfs_work   LAPACKE_NAME(dgerfs_work,DGERFS_WORK)
#define LAPACKE_cgerfs_work   LAPACKE_NAME(cgerfs_work,CGERFS_WORK)
#define LAPACKE_zgerfs_work   LAPACKE_NAME(zgerfs_work,ZGERFS_WORK)

#define LAPACKE_sgerfsx_work   LAPACKE_NAME(sgerfsx_work,SGERFSX_WORK)
#define LAPACKE_dgerfsx_work   LAPACKE_NAME(dgerfsx_work,DGERFSX_WORK)
#define LAPACKE_cgerfsx_work   LAPACKE_NAME(cgerfsx_work,CGERFSX_WORK)
#define LAPACKE_zgerfsx_work   LAPACKE_NAME(zgerfsx_work,ZGERFSX_WORK)

#define LAPACKE_sgerqf_work   LAPACKE_NAME(sgerqf_work,SGERQF_WORK)
#define LAPACKE_dgerqf_work   LAPACKE_NAME(dgerqf_work,DGERQF_WORK)
#define LAPACKE_cgerqf_work   LAPACKE_NAME(cgerqf_work,CGERQF_WORK)
#define LAPACKE_zgerqf_work   LAPACKE_NAME(zgerqf_work,ZGERQF_WORK)

#define LAPACKE_sgesdd_work   LAPACKE_NAME(sgesdd_work,SGESDD_WORK)
#define LAPACKE_dgesdd_work   LAPACKE_NAME(dgesdd_work,DGESDD_WORK)
#define LAPACKE_cgesdd_work   LAPACKE_NAME(cgesdd_work,CGESDD_WORK)
#define LAPACKE_zgesdd_work   LAPACKE_NAME(zgesdd_work,ZGESDD_WORK)

#define LAPACKE_sgesv_work   LAPACKE_NAME(sgesv_work,SGESV_WORK)
#define LAPACKE_dgesv_work   LAPACKE_NAME(dgesv_work,DGESV_WORK)
#define LAPACKE_cgesv_work   LAPACKE_NAME(cgesv_work,CGESV_WORK)
#define LAPACKE_zgesv_work   LAPACKE_NAME(zgesv_work,ZGESV_WORK)
#define LAPACKE_dsgesv_work   LAPACKE_NAME(dsgesv_work,DSGESV_WORK)
#define LAPACKE_zcgesv_work   LAPACKE_NAME(zcgesv_work,ZCGESV_WORK)

#define LAPACKE_sgesvd_work   LAPACKE_NAME(sgesvd_work,SGESVD_WORK)
#define LAPACKE_dgesvd_work   LAPACKE_NAME(dgesvd_work,DGESVD_WORK)
#define LAPACKE_cgesvd_work   LAPACKE_NAME(cgesvd_work,CGESVD_WORK)
#define LAPACKE_zgesvd_work   LAPACKE_NAME(zgesvd_work,ZGESVD_WORK)

#define LAPACKE_sgesvj_work   LAPACKE_NAME(sgesvj_work,SGESVJ_WORK)
#define LAPACKE_dgesvj_work   LAPACKE_NAME(dgesvj_work,DGESVJ_WORK)

#define LAPACKE_sgesvx_work   LAPACKE_NAME(sgesvx_work,SGESVX_WORK)
#define LAPACKE_dgesvx_work   LAPACKE_NAME(dgesvx_work,DGESVX_WORK)
#define LAPACKE_cgesvx_work   LAPACKE_NAME(cgesvx_work,CGESVX_WORK)
#define LAPACKE_zgesvx_work   LAPACKE_NAME(zgesvx_work,ZGESVX_WORK)

#define LAPACKE_sgesvxx_work   LAPACKE_NAME(sgesvxx_work,SGESVXX_WORK)
#define LAPACKE_dgesvxx_work   LAPACKE_NAME(dgesvxx_work,DGESVXX_WORK)
#define LAPACKE_cgesvxx_work   LAPACKE_NAME(cgesvxx_work,CGESVXX_WORK)
#define LAPACKE_zgesvxx_work   LAPACKE_NAME(zgesvxx_work,ZGESVXX_WORK)

#define LAPACKE_sgetf2_work   LAPACKE_NAME(sgetf2_work,SGETF2_WORK)
#define LAPACKE_dgetf2_work   LAPACKE_NAME(dgetf2_work,DGETF2_WORK)
#define LAPACKE_cgetf2_work   LAPACKE_NAME(cgetf2_work,CGETF2_WORK)
#define LAPACKE_zgetf2_work   LAPACKE_NAME(zgetf2_work,ZGETF2_WORK)

#define LAPACKE_sgetrf_work   LAPACKE_NAME(sgetrf_work,SGETRF_WORK)
#define LAPACKE_dgetrf_work   LAPACKE_NAME(dgetrf_work,DGETRF_WORK)
#define LAPACKE_cgetrf_work   LAPACKE_NAME(cgetrf_work,CGETRF_WORK)
#define LAPACKE_zgetrf_work   LAPACKE_NAME(zgetrf_work,ZGETRF_WORK)

#define LAPACKE_sgetri_work   LAPACKE_NAME(sgetri_work,SGETRI_WORK)
#define LAPACKE_dgetri_work   LAPACKE_NAME(dgetri_work,DGETRI_WORK)
#define LAPACKE_cgetri_work   LAPACKE_NAME(cgetri_work,CGETRI_WORK)
#define LAPACKE_zgetri_work   LAPACKE_NAME(zgetri_work,ZGETRI_WORK)

#define LAPACKE_sgetrs_work   LAPACKE_NAME(sgetrs_work,SGETRS_WORK)
#define LAPACKE_dgetrs_work   LAPACKE_NAME(dgetrs_work,DGETRS_WORK)
#define LAPACKE_cgetrs_work   LAPACKE_NAME(cgetrs_work,CGETRS_WORK)
#define LAPACKE_zgetrs_work   LAPACKE_NAME(zgetrs_work,ZGETRS_WORK)

#define LAPACKE_sggbak_work   LAPACKE_NAME(sggbak_work,SGGBAK_WORK)
#define LAPACKE_dggbak_work   LAPACKE_NAME(dggbak_work,DGGBAK_WORK)
#define LAPACKE_cggbak_work   LAPACKE_NAME(cggbak_work,CGGBAK_WORK)
#define LAPACKE_zggbak_work   LAPACKE_NAME(zggbak_work,ZGGBAK_WORK)

#define LAPACKE_sggbal_work   LAPACKE_NAME(sggbal_work,SGGBAL_WORK)
#define LAPACKE_dggbal_work   LAPACKE_NAME(dggbal_work,DGGBAL_WORK)
#define LAPACKE_cggbal_work   LAPACKE_NAME(cggbal_work,CGGBAL_WORK)
#define LAPACKE_zggbal_work   LAPACKE_NAME(zggbal_work,ZGGBAL_WORK)

#define LAPACKE_sgges_work   LAPACKE_NAME(sgges_work,SGGES_WORK)
#define LAPACKE_dgges_work   LAPACKE_NAME(dgges_work,DGGES_WORK)
#define LAPACKE_cgges_work   LAPACKE_NAME(cgges_work,CGGES_WORK)
#define LAPACKE_zgges_work   LAPACKE_NAME(zgges_work,ZGGES_WORK)

#define LAPACKE_sggesx_work   LAPACKE_NAME(sggesx_work,SGGESX_WORK)
#define LAPACKE_dggesx_work   LAPACKE_NAME(dggesx_work,DGGESX_WORK)
#define LAPACKE_cggesx_work   LAPACKE_NAME(cggesx_work,CGGESX_WORK)
#define LAPACKE_zggesx_work   LAPACKE_NAME(zggesx_work,ZGGESX_WORK)

#define LAPACKE_sggev_work   LAPACKE_NAME(sggev_work,SGGEV_WORK)
#define LAPACKE_dggev_work   LAPACKE_NAME(dggev_work,DGGEV_WORK)
#define LAPACKE_cggev_work   LAPACKE_NAME(cggev_work,CGGEV_WORK)
#define LAPACKE_zggev_work   LAPACKE_NAME(zggev_work,ZGGEV_WORK)

#define LAPACKE_sggevx_work   LAPACKE_NAME(sggevx_work,SGGEVX_WORK)
#define LAPACKE_dggevx_work   LAPACKE_NAME(dggevx_work,DGGEVX_WORK)
#define LAPACKE_cggevx_work   LAPACKE_NAME(cggevx_work,CGGEVX_WORK)
#define LAPACKE_zggevx_work   LAPACKE_NAME(zggevx_work,ZGGEVX_WORK)

#define LAPACKE_sggglm_work   LAPACKE_NAME(sggglm_work,SGGGLM_WORK)
#define LAPACKE_dggglm_work   LAPACKE_NAME(dggglm_work,DGGGLM_WORK)
#define LAPACKE_cggglm_work   LAPACKE_NAME(cggglm_work,CGGGLM_WORK)
#define LAPACKE_zggglm_work   LAPACKE_NAME(zggglm_work,ZGGGLM_WORK)

#define LAPACKE_sgghrd_work   LAPACKE_NAME(sgghrd_work,SGGHRD_WORK)
#define LAPACKE_dgghrd_work   LAPACKE_NAME(dgghrd_work,DGGHRD_WORK)
#define LAPACKE_cgghrd_work   LAPACKE_NAME(cgghrd_work,CGGHRD_WORK)
#define LAPACKE_zgghrd_work   LAPACKE_NAME(zgghrd_work,ZGGHRD_WORK)

#define LAPACKE_sgglse_work   LAPACKE_NAME(sgglse_work,SGGLSE_WORK)
#define LAPACKE_dgglse_work   LAPACKE_NAME(dgglse_work,DGGLSE_WORK)
#define LAPACKE_cgglse_work   LAPACKE_NAME(cgglse_work,CGGLSE_WORK)
#define LAPACKE_zgglse_work   LAPACKE_NAME(zgglse_work,ZGGLSE_WORK)

#define LAPACKE_sggqrf_work   LAPACKE_NAME(sggqrf_work,SGGQRF_WORK)
#define LAPACKE_dggqrf_work   LAPACKE_NAME(dggqrf_work,DGGQRF_WORK)
#define LAPACKE_cggqrf_work   LAPACKE_NAME(cggqrf_work,CGGQRF_WORK)
#define LAPACKE_zggqrf_work   LAPACKE_NAME(zggqrf_work,ZGGQRF_WORK)

#define LAPACKE_sggrqf_work   LAPACKE_NAME(sggrqf_work,SGGRQF_WORK)
#define LAPACKE_dggrqf_work   LAPACKE_NAME(dggrqf_work,DGGRQF_WORK)
#define LAPACKE_cggrqf_work   LAPACKE_NAME(cggrqf_work,CGGRQF_WORK)
#define LAPACKE_zggrqf_work   LAPACKE_NAME(zggrqf_work,ZGGRQF_WORK)

#define LAPACKE_sggsvd_work   LAPACKE_NAME(sggsvd_work,SGGSVD_WORK)
#define LAPACKE_dggsvd_work   LAPACKE_NAME(dggsvd_work,DGGSVD_WORK)
#define LAPACKE_cggsvd_work   LAPACKE_NAME(cggsvd_work,CGGSVD_WORK)
#define LAPACKE_zggsvd_work   LAPACKE_NAME(zggsvd_work,ZGGSVD_WORK)

#define LAPACKE_sggsvp_work   LAPACKE_NAME(sggsvp_work,SGGSVP_WORK)
#define LAPACKE_dggsvp_work   LAPACKE_NAME(dggsvp_work,DGGSVP_WORK)
#define LAPACKE_cggsvp_work   LAPACKE_NAME(cggsvp_work,CGGSVP_WORK)
#define LAPACKE_zggsvp_work   LAPACKE_NAME(zggsvp_work,ZGGSVP_WORK)

#define LAPACKE_sgtcon_work   LAPACKE_NAME(sgtcon_work,SGTCON_WORK)
#define LAPACKE_dgtcon_work   LAPACKE_NAME(dgtcon_work,DGTCON_WORK)
#define LAPACKE_cgtcon_work   LAPACKE_NAME(cgtcon_work,CGTCON_WORK)
#define LAPACKE_zgtcon_work   LAPACKE_NAME(zgtcon_work,ZGTCON_WORK)

#define LAPACKE_sgtrfs_work   LAPACKE_NAME(sgtrfs_work,SGTRFS_WORK)
#define LAPACKE_dgtrfs_work   LAPACKE_NAME(dgtrfs_work,DGTRFS_WORK)
#define LAPACKE_cgtrfs_work   LAPACKE_NAME(cgtrfs_work,CGTRFS_WORK)
#define LAPACKE_zgtrfs_work   LAPACKE_NAME(zgtrfs_work,ZGTRFS_WORK)

#define LAPACKE_sgtsv_work   LAPACKE_NAME(sgtsv_work,SGTSV_WORK)
#define LAPACKE_dgtsv_work   LAPACKE_NAME(dgtsv_work,DGTSV_WORK)
#define LAPACKE_cgtsv_work   LAPACKE_NAME(cgtsv_work,CGTSV_WORK)
#define LAPACKE_zgtsv_work   LAPACKE_NAME(zgtsv_work,ZGTSV_WORK)

#define LAPACKE_sgtsvx_work   LAPACKE_NAME(sgtsvx_work,SGTSVX_WORK)
#define LAPACKE_dgtsvx_work   LAPACKE_NAME(dgtsvx_work,DGTSVX_WORK)
#define LAPACKE_cgtsvx_work   LAPACKE_NAME(cgtsvx_work,CGTSVX_WORK)
#define LAPACKE_zgtsvx_work   LAPACKE_NAME(zgtsvx_work,ZGTSVX_WORK)

#define LAPACKE_sgttrf_work   LAPACKE_NAME(sgttrf_work,SGTTRF_WORK)
#define LAPACKE_dgttrf_work   LAPACKE_NAME(dgttrf_work,DGTTRF_WORK)
#define LAPACKE_cgttrf_work   LAPACKE_NAME(cgttrf_work,CGTTRF_WORK)
#define LAPACKE_zgttrf_work   LAPACKE_NAME(zgttrf_work,ZGTTRF_WORK)

#define LAPACKE_sgttrs_work   LAPACKE_NAME(sgttrs_work,SGTTRS_WORK)
#define LAPACKE_dgttrs_work   LAPACKE_NAME(dgttrs_work,DGTTRS_WORK)
#define LAPACKE_cgttrs_work   LAPACKE_NAME(cgttrs_work,CGTTRS_WORK)
#define LAPACKE_zgttrs_work   LAPACKE_NAME(zgttrs_work,ZGTTRS_WORK)

#define LAPACKE_chbev_work   LAPACKE_NAME(chbev_work,CHBEV_WORK)
#define LAPACKE_zhbev_work   LAPACKE_NAME(zhbev_work,ZHBEV_WORK)

#define LAPACKE_chbevd_work   LAPACKE_NAME(chbevd_work,CHBEVD_WORK)
#define LAPACKE_zhbevd_work   LAPACKE_NAME(zhbevd_work,ZHBEVD_WORK)

#define LAPACKE_chbevx_work   LAPACKE_NAME(chbevx_work,CHBEVX_WORK)
#define LAPACKE_zhbevx_work   LAPACKE_NAME(zhbevx_work,ZHBEVX_WORK)

#define LAPACKE_chbgst_work   LAPACKE_NAME(chbgst_work,CHBGST_WORK)
#define LAPACKE_zhbgst_work   LAPACKE_NAME(zhbgst_work,ZHBGST_WORK)

#define LAPACKE_chbgv_work   LAPACKE_NAME(chbgv_work,CHBGV_WORK)
#define LAPACKE_zhbgv_work   LAPACKE_NAME(zhbgv_work,ZHBGV_WORK)

#define LAPACKE_chbgvd_work   LAPACKE_NAME(chbgvd_work,CHBGVD_WORK)
#define LAPACKE_zhbgvd_work   LAPACKE_NAME(zhbgvd_work,ZHBGVD_WORK)

#define LAPACKE_chbgvx_work   LAPACKE_NAME(chbgvx_work,CHBGVX_WORK)
#define LAPACKE_zhbgvx_work   LAPACKE_NAME(zhbgvx_work,ZHBGVX_WORK)

#define LAPACKE_chbtrd_work   LAPACKE_NAME(chbtrd_work,CHBTRD_WORK)
#define LAPACKE_zhbtrd_work   LAPACKE_NAME(zhbtrd_work,ZHBTRD_WORK)

#define LAPACKE_checon_work   LAPACKE_NAME(checon_work,CHECON_WORK)
#define LAPACKE_zhecon_work   LAPACKE_NAME(zhecon_work,ZHECON_WORK)

#define LAPACKE_cheequb_work   LAPACKE_NAME(cheequb_work,CHEEQUB_WORK)
#define LAPACKE_zheequb_work   LAPACKE_NAME(zheequb_work,ZHEEQUB_WORK)

#define LAPACKE_cheev_work   LAPACKE_NAME(cheev_work,CHEEV_WORK)
#define LAPACKE_zheev_work   LAPACKE_NAME(zheev_work,ZHEEV_WORK)

#define LAPACKE_cheevd_work   LAPACKE_NAME(cheevd_work,CHEEVD_WORK)
#define LAPACKE_zheevd_work   LAPACKE_NAME(zheevd_work,ZHEEVD_WORK)

#define LAPACKE_cheevr_work   LAPACKE_NAME(cheevr_work,CHEEVR_WORK)
#define LAPACKE_zheevr_work   LAPACKE_NAME(zheevr_work,ZHEEVR_WORK)

#define LAPACKE_cheevx_work   LAPACKE_NAME(cheevx_work,CHEEVX_WORK)
#define LAPACKE_zheevx_work   LAPACKE_NAME(zheevx_work,ZHEEVX_WORK)

#define LAPACKE_chegst_work   LAPACKE_NAME(chegst_work,CHEGST_WORK)
#define LAPACKE_zhegst_work   LAPACKE_NAME(zhegst_work,ZHEGST_WORK)

#define LAPACKE_chegv_work   LAPACKE_NAME(chegv_work,CHEGV_WORK)
#define LAPACKE_zhegv_work   LAPACKE_NAME(zhegv_work,ZHEGV_WORK)

#define LAPACKE_chegvd_work   LAPACKE_NAME(chegvd_work,CHEGVD_WORK)
#define LAPACKE_zhegvd_work   LAPACKE_NAME(zhegvd_work,ZHEGVD_WORK)

#define LAPACKE_chegvx_work   LAPACKE_NAME(chegvx_work,CHEGVX_WORK)
#define LAPACKE_zhegvx_work   LAPACKE_NAME(zhegvx_work,ZHEGVX_WORK)

#define LAPACKE_cherfs_work   LAPACKE_NAME(cherfs_work,CHERFS_WORK)
#define LAPACKE_zherfs_work   LAPACKE_NAME(zherfs_work,ZHERFS_WORK)

#define LAPACKE_cherfsx_work   LAPACKE_NAME(cherfsx_work,CHERFSX_WORK)
#define LAPACKE_zherfsx_work   LAPACKE_NAME(zherfsx_work,ZHERFSX_WORK)

#define LAPACKE_chesv_work   LAPACKE_NAME(chesv_work,CHESV_WORK)
#define LAPACKE_zhesv_work   LAPACKE_NAME(zhesv_work,ZHESV_WORK)

#define LAPACKE_chesvx_work   LAPACKE_NAME(chesvx_work,CHESVX_WORK)
#define LAPACKE_zhesvx_work   LAPACKE_NAME(zhesvx_work,ZHESVX_WORK)

#define LAPACKE_chesvxx_work   LAPACKE_NAME(chesvxx_work,CHESVXX_WORK)
#define LAPACKE_zhesvxx_work   LAPACKE_NAME(zhesvxx_work,ZHESVXX_WORK)

#define LAPACKE_chetrd_work   LAPACKE_NAME(chetrd_work,CHETRD_WORK)
#define LAPACKE_zhetrd_work   LAPACKE_NAME(zhetrd_work,ZHETRD_WORK)

#define LAPACKE_chetrf_work   LAPACKE_NAME(chetrf_work,CHETRF_WORK)
#define LAPACKE_zhetrf_work   LAPACKE_NAME(zhetrf_work,ZHETRF_WORK)

#define LAPACKE_chetri_work   LAPACKE_NAME(chetri_work,CHETRI_WORK)
#define LAPACKE_zhetri_work   LAPACKE_NAME(zhetri_work,ZHETRI_WORK)

#define LAPACKE_chetrs_work   LAPACKE_NAME(chetrs_work,CHETRS_WORK)
#define LAPACKE_zhetrs_work   LAPACKE_NAME(zhetrs_work,ZHETRS_WORK)

#define LAPACKE_chfrk_work   LAPACKE_NAME(chfrk_work,CHFRK_WORK)
#define LAPACKE_zhfrk_work   LAPACKE_NAME(zhfrk_work,ZHFRK_WORK)

#define LAPACKE_shgeqz_work   LAPACKE_NAME(shgeqz_work,SHGEQZ_WORK)
#define LAPACKE_dhgeqz_work   LAPACKE_NAME(dhgeqz_work,DHGEQZ_WORK)
#define LAPACKE_chgeqz_work   LAPACKE_NAME(chgeqz_work,CHGEQZ_WORK)
#define LAPACKE_zhgeqz_work   LAPACKE_NAME(zhgeqz_work,ZHGEQZ_WORK)

#define LAPACKE_chpcon_work   LAPACKE_NAME(chpcon_work,CHPCON_WORK)
#define LAPACKE_zhpcon_work   LAPACKE_NAME(zhpcon_work,ZHPCON_WORK)

#define LAPACKE_chpev_work   LAPACKE_NAME(chpev_work,CHPEV_WORK)
#define LAPACKE_zhpev_work   LAPACKE_NAME(zhpev_work,ZHPEV_WORK)

#define LAPACKE_chpevd_work   LAPACKE_NAME(chpevd_work,CHPEVD_WORK)
#define LAPACKE_zhpevd_work   LAPACKE_NAME(zhpevd_work,ZHPEVD_WORK)

#define LAPACKE_chpevx_work   LAPACKE_NAME(chpevx_work,CHPEVX_WORK)
#define LAPACKE_zhpevx_work   LAPACKE_NAME(zhpevx_work,ZHPEVX_WORK)

#define LAPACKE_chpgst_work   LAPACKE_NAME(chpgst_work,CHPGST_WORK)
#define LAPACKE_zhpgst_work   LAPACKE_NAME(zhpgst_work,ZHPGST_WORK)

#define LAPACKE_chpgv_work   LAPACKE_NAME(chpgv_work,CHPGV_WORK)
#define LAPACKE_zhpgv_work   LAPACKE_NAME(zhpgv_work,ZHPGV_WORK)

#define LAPACKE_chpgvd_work   LAPACKE_NAME(chpgvd_work,CHPGVD_WORK)
#define LAPACKE_zhpgvd_work   LAPACKE_NAME(zhpgvd_work,ZHPGVD_WORK)

#define LAPACKE_chpgvx_work   LAPACKE_NAME(chpgvx_work,CHPGVX_WORK)
#define LAPACKE_zhpgvx_work   LAPACKE_NAME(zhpgvx_work,ZHPGVX_WORK)

#define LAPACKE_chprfs_work   LAPACKE_NAME(chprfs_work,CHPRFS_WORK)
#define LAPACKE_zhprfs_work   LAPACKE_NAME(zhprfs_work,ZHPRFS_WORK)

#define LAPACKE_chpsv_work   LAPACKE_NAME(chpsv_work,CHPSV_WORK)
#define LAPACKE_zhpsv_work   LAPACKE_NAME(zhpsv_work,ZHPSV_WORK)

#define LAPACKE_chpsvx_work   LAPACKE_NAME(chpsvx_work,CHPSVX_WORK)
#define LAPACKE_zhpsvx_work   LAPACKE_NAME(zhpsvx_work,ZHPSVX_WORK)

#define LAPACKE_chptrd_work   LAPACKE_NAME(chptrd_work,CHPTRD_WORK)
#define LAPACKE_zhptrd_work   LAPACKE_NAME(zhptrd_work,ZHPTRD_WORK)

#define LAPACKE_chptrf_work   LAPACKE_NAME(chptrf_work,CHPTRF_WORK)
#define LAPACKE_zhptrf_work   LAPACKE_NAME(zhptrf_work,ZHPTRF_WORK)

#define LAPACKE_chptri_work   LAPACKE_NAME(chptri_work,CHPTRI_WORK)
#define LAPACKE_zhptri_work   LAPACKE_NAME(zhptri_work,ZHPTRI_WORK)

#define LAPACKE_chptrs_work   LAPACKE_NAME(chptrs_work,CHPTRS_WORK)
#define LAPACKE_zhptrs_work   LAPACKE_NAME(zhptrs_work,ZHPTRS_WORK)

#define LAPACKE_shsein_work   LAPACKE_NAME(shsein_work,SHSEIN_WORK)
#define LAPACKE_dhsein_work   LAPACKE_NAME(dhsein_work,DHSEIN_WORK)
#define LAPACKE_chsein_work   LAPACKE_NAME(chsein_work,CHSEIN_WORK)
#define LAPACKE_zhsein_work   LAPACKE_NAME(zhsein_work,ZHSEIN_WORK)

#define LAPACKE_shseqr_work   LAPACKE_NAME(shseqr_work,SHSEQR_WORK)
#define LAPACKE_dhseqr_work   LAPACKE_NAME(dhseqr_work,DHSEQR_WORK)
#define LAPACKE_chseqr_work   LAPACKE_NAME(chseqr_work,CHSEQR_WORK)
#define LAPACKE_zhseqr_work   LAPACKE_NAME(zhseqr_work,ZHSEQR_WORK)

#define LAPACKE_clacgv_work   LAPACKE_NAME(clacgv_work,CLACGV_WORK)
#define LAPACKE_zlacgv_work   LAPACKE_NAME(zlacgv_work,ZLACGV_WORK)

#define LAPACKE_slacpy_work   LAPACKE_NAME(slacpy_work,SLACPY_WORK)
#define LAPACKE_dlacpy_work   LAPACKE_NAME(dlacpy_work,DLACPY_WORK)
#define LAPACKE_clacpy_work   LAPACKE_NAME(clacpy_work,CLACPY_WORK)
#define LAPACKE_zlacpy_work   LAPACKE_NAME(zlacpy_work,ZLACPY_WORK)

#define LAPACKE_zlag2c_work   LAPACKE_NAME(zlag2c_work,ZLAG2C_WORK)

#define LAPACKE_slag2d_work   LAPACKE_NAME(slag2d_work,SLAG2D_WORK)

#define LAPACKE_dlag2s_work   LAPACKE_NAME(dlag2s_work,DLAG2S_WORK)

#define LAPACKE_clag2z_work   LAPACKE_NAME(clag2z_work,CLAG2Z_WORK)

#define LAPACKE_slagge_work   LAPACKE_NAME(slagge_work,SLAGGE_WORK)
#define LAPACKE_dlagge_work   LAPACKE_NAME(dlagge_work,DLAGGE_WORK)
#define LAPACKE_clagge_work   LAPACKE_NAME(clagge_work,CLAGGE_WORK)
#define LAPACKE_zlagge_work   LAPACKE_NAME(zlagge_work,ZLAGGE_WORK)

#define LAPACKE_slamch_work   LAPACKE_NAME(slamch_work,SLAMCH_WORK)
#define LAPACKE_dlamch_work   LAPACKE_NAME(dlamch_work,DLAMCH_WORK)

#define LAPACKE_slange_work   LAPACKE_NAME(slange_work,SLANGE_WORK)
#define LAPACKE_dlange_work   LAPACKE_NAME(dlange_work,DLANGE_WORK)
#define LAPACKE_clange_work   LAPACKE_NAME(clange_work,CLANGE_WORK)
#define LAPACKE_zlange_work   LAPACKE_NAME(zlange_work,ZLANGE_WORK)

#define LAPACKE_clanhe_work   LAPACKE_NAME(clanhe_work,CLANHE_WORK)
#define LAPACKE_zlanhe_work   LAPACKE_NAME(zlanhe_work,ZLANHE_WORK)

#define LAPACKE_slansy_work   LAPACKE_NAME(slansy_work,SLANSY_WORK)
#define LAPACKE_dlansy_work   LAPACKE_NAME(dlansy_work,DLANSY_WORK)
#define LAPACKE_clansy_work   LAPACKE_NAME(clansy_work,CLANSY_WORK)
#define LAPACKE_zlansy_work   LAPACKE_NAME(zlansy_work,ZLANSY_WORK)

#define LAPACKE_slantr_work   LAPACKE_NAME(slantr_work,SLANTR_WORK)
#define LAPACKE_dlantr_work   LAPACKE_NAME(dlantr_work,DLANTR_WORK)
#define LAPACKE_clantr_work   LAPACKE_NAME(clantr_work,CLANTR_WORK)
#define LAPACKE_zlantr_work   LAPACKE_NAME(zlantr_work,ZLANTR_WORK)

#define LAPACKE_slarfb_work   LAPACKE_NAME(slarfb_work,SLARFB_WORK)
#define LAPACKE_dlarfb_work   LAPACKE_NAME(dlarfb_work,DLARFB_WORK)
#define LAPACKE_clarfb_work   LAPACKE_NAME(clarfb_work,CLARFB_WORK)
#define LAPACKE_zlarfb_work   LAPACKE_NAME(zlarfb_work,ZLARFB_WORK)

#define LAPACKE_slarfg_work   LAPACKE_NAME(slarfg_work,SLARFG_WORK)
#define LAPACKE_dlarfg_work   LAPACKE_NAME(dlarfg_work,DLARFG_WORK)
#define LAPACKE_clarfg_work   LAPACKE_NAME(clarfg_work,CLARFG_WORK)
#define LAPACKE_zlarfg_work   LAPACKE_NAME(zlarfg_work,ZLARFG_WORK)

#define LAPACKE_slarft_work   LAPACKE_NAME(slarft_work,SLARFT_WORK)
#define LAPACKE_dlarft_work   LAPACKE_NAME(dlarft_work,DLARFT_WORK)
#define LAPACKE_clarft_work   LAPACKE_NAME(clarft_work,CLARFT_WORK)
#define LAPACKE_zlarft_work   LAPACKE_NAME(zlarft_work,ZLARFT_WORK)

#define LAPACKE_slarfx_work   LAPACKE_NAME(slarfx_work,SLARFX_WORK)
#define LAPACKE_dlarfx_work   LAPACKE_NAME(dlarfx_work,DLARFX_WORK)
#define LAPACKE_clarfx_work   LAPACKE_NAME(clarfx_work,CLARFX_WORK)
#define LAPACKE_zlarfx_work   LAPACKE_NAME(zlarfx_work,ZLARFX_WORK)

#define LAPACKE_slarnv_work   LAPACKE_NAME(slarnv_work,SLARNV_WORK)
#define LAPACKE_dlarnv_work   LAPACKE_NAME(dlarnv_work,DLARNV_WORK)
#define LAPACKE_clarnv_work   LAPACKE_NAME(clarnv_work,CLARNV_WORK)
#define LAPACKE_zlarnv_work   LAPACKE_NAME(zlarnv_work,ZLARNV_WORK)

#define LAPACKE_slaset_work   LAPACKE_NAME(slaset_work,SLASET_WORK)
#define LAPACKE_dlaset_work   LAPACKE_NAME(dlaset_work,DLASET_WORK)
#define LAPACKE_claset_work   LAPACKE_NAME(claset_work,CLASET_WORK)
#define LAPACKE_zlaset_work   LAPACKE_NAME(zlaset_work,ZLASET_WORK)

#define LAPACKE_slasrt_work   LAPACKE_NAME(slasrt_work,SLASRT_WORK)
#define LAPACKE_dlasrt_work   LAPACKE_NAME(dlasrt_work,DLASRT_WORK)

#define LAPACKE_slaswp_work   LAPACKE_NAME(slaswp_work,SLASWP_WORK)
#define LAPACKE_dlaswp_work   LAPACKE_NAME(dlaswp_work,DLASWP_WORK)
#define LAPACKE_claswp_work   LAPACKE_NAME(claswp_work,CLASWP_WORK)
#define LAPACKE_zlaswp_work   LAPACKE_NAME(zlaswp_work,ZLASWP_WORK)

#define LAPACKE_slatms_work   LAPACKE_NAME(slatms_work,SLATMS_WORK)
#define LAPACKE_dlatms_work   LAPACKE_NAME(dlatms_work,DLATMS_WORK)
#define LAPACKE_clatms_work   LAPACKE_NAME(clatms_work,CLATMS_WORK)
#define LAPACKE_zlatms_work   LAPACKE_NAME(zlatms_work,ZLATMS_WORK)

#define LAPACKE_slauum_work   LAPACKE_NAME(slauum_work,SLAUUM_WORK)
#define LAPACKE_dlauum_work   LAPACKE_NAME(dlauum_work,DLAUUM_WORK)
#define LAPACKE_clauum_work   LAPACKE_NAME(clauum_work,CLAUUM_WORK)
#define LAPACKE_zlauum_work   LAPACKE_NAME(zlauum_work,ZLAUUM_WORK)

#define LAPACKE_sopgtr_work   LAPACKE_NAME(sopgtr_work,SOPGTR_WORK)
#define LAPACKE_dopgtr_work   LAPACKE_NAME(dopgtr_work,DOPGTR_WORK)

#define LAPACKE_sopmtr_work   LAPACKE_NAME(sopmtr_work,SOPMTR_WORK)
#define LAPACKE_dopmtr_work   LAPACKE_NAME(dopmtr_work,DOPMTR_WORK)

#define LAPACKE_sorgbr_work   LAPACKE_NAME(sorgbr_work,SORGBR_WORK)
#define LAPACKE_dorgbr_work   LAPACKE_NAME(dorgbr_work,DORGBR_WORK)

#define LAPACKE_sorghr_work   LAPACKE_NAME(sorghr_work,SORGHR_WORK)
#define LAPACKE_dorghr_work   LAPACKE_NAME(dorghr_work,DORGHR_WORK)

#define LAPACKE_sorglq_work   LAPACKE_NAME(sorglq_work,SORGLQ_WORK)
#define LAPACKE_dorglq_work   LAPACKE_NAME(dorglq_work,DORGLQ_WORK)

#define LAPACKE_sorgql_work   LAPACKE_NAME(sorgql_work,SORGQL_WORK)
#define LAPACKE_dorgql_work   LAPACKE_NAME(dorgql_work,DORGQL_WORK)

#define LAPACKE_sorgqr_work   LAPACKE_NAME(sorgqr_work,SORGQR_WORK)
#define LAPACKE_dorgqr_work   LAPACKE_NAME(dorgqr_work,DORGQR_WORK)

#define LAPACKE_sorgrq_work   LAPACKE_NAME(sorgrq_work,SORGRQ_WORK)
#define LAPACKE_dorgrq_work   LAPACKE_NAME(dorgrq_work,DORGRQ_WORK)

#define LAPACKE_sorgtr_work   LAPACKE_NAME(sorgtr_work,SORGTR_WORK)
#define LAPACKE_dorgtr_work   LAPACKE_NAME(dorgtr_work,DORGTR_WORK)

#define LAPACKE_sormbr_work   LAPACKE_NAME(sormbr_work,SORMBR_WORK)
#define LAPACKE_dormbr_work   LAPACKE_NAME(dormbr_work,DORMBR_WORK)

#define LAPACKE_sormhr_work   LAPACKE_NAME(sormhr_work,SORMHR_WORK)
#define LAPACKE_dormhr_work   LAPACKE_NAME(dormhr_work,DORMHR_WORK)

#define LAPACKE_sormlq_work   LAPACKE_NAME(sormlq_work,SORMLQ_WORK)
#define LAPACKE_dormlq_work   LAPACKE_NAME(dormlq_work,DORMLQ_WORK)

#define LAPACKE_sormql_work   LAPACKE_NAME(sormql_work,SORMQL_WORK)
#define LAPACKE_dormql_work   LAPACKE_NAME(dormql_work,DORMQL_WORK)

#define LAPACKE_sormqr_work   LAPACKE_NAME(sormqr_work,SORMQR_WORK)
#define LAPACKE_dormqr_work   LAPACKE_NAME(dormqr_work,DORMQR_WORK)

#define LAPACKE_sormrq_work   LAPACKE_NAME(sormrq_work,SORMRQ_WORK)
#define LAPACKE_dormrq_work   LAPACKE_NAME(dormrq_work,DORMRQ_WORK)

#define LAPACKE_sormrz_work   LAPACKE_NAME(sormrz_work,SORMRZ_WORK)
#define LAPACKE_dormrz_work   LAPACKE_NAME(dormrz_work,DORMRZ_WORK)

#define LAPACKE_sormtr_work   LAPACKE_NAME(sormtr_work,SORMTR_WORK)
#define LAPACKE_dormtr_work   LAPACKE_NAME(dormtr_work,DORMTR_WORK)

#define LAPACKE_spbcon_work   LAPACKE_NAME(spbcon_work,SPBCON_WORK)
#define LAPACKE_dpbcon_work   LAPACKE_NAME(dpbcon_work,DPBCON_WORK)
#define LAPACKE_cpbcon_work   LAPACKE_NAME(cpbcon_work,CPBCON_WORK)
#define LAPACKE_zpbcon_work   LAPACKE_NAME(zpbcon_work,ZPBCON_WORK)

#define LAPACKE_spbequ_work   LAPACKE_NAME(spbequ_work,SPBEQU_WORK)
#define LAPACKE_dpbequ_work   LAPACKE_NAME(dpbequ_work,DPBEQU_WORK)
#define LAPACKE_cpbequ_work   LAPACKE_NAME(cpbequ_work,CPBEQU_WORK)
#define LAPACKE_zpbequ_work   LAPACKE_NAME(zpbequ_work,ZPBEQU_WORK)

#define LAPACKE_spbrfs_work   LAPACKE_NAME(spbrfs_work,SPBRFS_WORK)
#define LAPACKE_dpbrfs_work   LAPACKE_NAME(dpbrfs_work,DPBRFS_WORK)
#define LAPACKE_cpbrfs_work   LAPACKE_NAME(cpbrfs_work,CPBRFS_WORK)
#define LAPACKE_zpbrfs_work   LAPACKE_NAME(zpbrfs_work,ZPBRFS_WORK)

#define LAPACKE_spbstf_work   LAPACKE_NAME(spbstf_work,SPBSTF_WORK)
#define LAPACKE_dpbstf_work   LAPACKE_NAME(dpbstf_work,DPBSTF_WORK)
#define LAPACKE_cpbstf_work   LAPACKE_NAME(cpbstf_work,CPBSTF_WORK)
#define LAPACKE_zpbstf_work   LAPACKE_NAME(zpbstf_work,ZPBSTF_WORK)

#define LAPACKE_spbsv_work   LAPACKE_NAME(spbsv_work,SPBSV_WORK)
#define LAPACKE_dpbsv_work   LAPACKE_NAME(dpbsv_work,DPBSV_WORK)
#define LAPACKE_cpbsv_work   LAPACKE_NAME(cpbsv_work,CPBSV_WORK)
#define LAPACKE_zpbsv_work   LAPACKE_NAME(zpbsv_work,ZPBSV_WORK)

#define LAPACKE_spbsvx_work   LAPACKE_NAME(spbsvx_work,SPBSVX_WORK)
#define LAPACKE_dpbsvx_work   LAPACKE_NAME(dpbsvx_work,DPBSVX_WORK)
#define LAPACKE_cpbsvx_work   LAPACKE_NAME(cpbsvx_work,CPBSVX_WORK)
#define LAPACKE_zpbsvx_work   LAPACKE_NAME(zpbsvx_work,ZPBSVX_WORK)

#define LAPACKE_spbtrf_work   LAPACKE_NAME(spbtrf_work,SPBTRF_WORK)
#define LAPACKE_dpbtrf_work   LAPACKE_NAME(dpbtrf_work,DPBTRF_WORK)
#define LAPACKE_cpbtrf_work   LAPACKE_NAME(cpbtrf_work,CPBTRF_WORK)
#define LAPACKE_zpbtrf_work   LAPACKE_NAME(zpbtrf_work,ZPBTRF_WORK)

#define LAPACKE_spbtrs_work   LAPACKE_NAME(spbtrs_work,SPBTRS_WORK)
#define LAPACKE_dpbtrs_work   LAPACKE_NAME(dpbtrs_work,DPBTRS_WORK)
#define LAPACKE_cpbtrs_work   LAPACKE_NAME(cpbtrs_work,CPBTRS_WORK)
#define LAPACKE_zpbtrs_work   LAPACKE_NAME(zpbtrs_work,ZPBTRS_WORK)

#define LAPACKE_spftrf_work   LAPACKE_NAME(spftrf_work,SPFTRF_WORK)
#define LAPACKE_dpftrf_work   LAPACKE_NAME(dpftrf_work,DPFTRF_WORK)
#define LAPACKE_cpftrf_work   LAPACKE_NAME(cpftrf_work,CPFTRF_WORK)
#define LAPACKE_zpftrf_work   LAPACKE_NAME(zpftrf_work,ZPFTRF_WORK)

#define LAPACKE_spftri_work   LAPACKE_NAME(spftri_work,SPFTRI_WORK)
#define LAPACKE_dpftri_work   LAPACKE_NAME(dpftri_work,DPFTRI_WORK)
#define LAPACKE_cpftri_work   LAPACKE_NAME(cpftri_work,CPFTRI_WORK)
#define LAPACKE_zpftri_work   LAPACKE_NAME(zpftri_work,ZPFTRI_WORK)

#define LAPACKE_spftrs_work   LAPACKE_NAME(spftrs_work,SPFTRS_WORK)
#define LAPACKE_dpftrs_work   LAPACKE_NAME(dpftrs_work,DPFTRS_WORK)
#define LAPACKE_cpftrs_work   LAPACKE_NAME(cpftrs_work,CPFTRS_WORK)
#define LAPACKE_zpftrs_work   LAPACKE_NAME(zpftrs_work,ZPFTRS_WORK)

#define LAPACKE_spocon_work   LAPACKE_NAME(spocon_work,SPOCON_WORK)
#define LAPACKE_dpocon_work   LAPACKE_NAME(dpocon_work,DPOCON_WORK)
#define LAPACKE_cpocon_work   LAPACKE_NAME(cpocon_work,CPOCON_WORK)
#define LAPACKE_zpocon_work   LAPACKE_NAME(zpocon_work,ZPOCON_WORK)

#define LAPACKE_spoequ_work   LAPACKE_NAME(spoequ_work,SPOEQU_WORK)
#define LAPACKE_dpoequ_work   LAPACKE_NAME(dpoequ_work,DPOEQU_WORK)
#define LAPACKE_cpoequ_work   LAPACKE_NAME(cpoequ_work,CPOEQU_WORK)
#define LAPACKE_zpoequ_work   LAPACKE_NAME(zpoequ_work,ZPOEQU_WORK)

#define LAPACKE_spoequb_work   LAPACKE_NAME(spoequb_work,SPOEQUB_WORK)
#define LAPACKE_dpoequb_work   LAPACKE_NAME(dpoequb_work,DPOEQUB_WORK)
#define LAPACKE_cpoequb_work   LAPACKE_NAME(cpoequb_work,CPOEQUB_WORK)
#define LAPACKE_zpoequb_work   LAPACKE_NAME(zpoequb_work,ZPOEQUB_WORK)

#define LAPACKE_sporfs_work   LAPACKE_NAME(sporfs_work,SPORFS_WORK)
#define LAPACKE_dporfs_work   LAPACKE_NAME(dporfs_work,DPORFS_WORK)
#define LAPACKE_cporfs_work   LAPACKE_NAME(cporfs_work,CPORFS_WORK)
#define LAPACKE_zporfs_work   LAPACKE_NAME(zporfs_work,ZPORFS_WORK)

#define LAPACKE_sporfsx_work   LAPACKE_NAME(sporfsx_work,SPORFSX_WORK)
#define LAPACKE_dporfsx_work   LAPACKE_NAME(dporfsx_work,DPORFSX_WORK)
#define LAPACKE_cporfsx_work   LAPACKE_NAME(cporfsx_work,CPORFSX_WORK)
#define LAPACKE_zporfsx_work   LAPACKE_NAME(zporfsx_work,ZPORFSX_WORK)

#define LAPACKE_sposv_work   LAPACKE_NAME(sposv_work,SPOSV_WORK)
#define LAPACKE_dposv_work   LAPACKE_NAME(dposv_work,DPOSV_WORK)
#define LAPACKE_cposv_work   LAPACKE_NAME(cposv_work,CPOSV_WORK)
#define LAPACKE_zposv_work   LAPACKE_NAME(zposv_work,ZPOSV_WORK)
#define LAPACKE_dsposv_work   LAPACKE_NAME(dsposv_work,DSPOSV_WORK)
#define LAPACKE_zcposv_work   LAPACKE_NAME(zcposv_work,ZCPOSV_WORK)

#define LAPACKE_sposvx_work   LAPACKE_NAME(sposvx_work,SPOSVX_WORK)
#define LAPACKE_dposvx_work   LAPACKE_NAME(dposvx_work,DPOSVX_WORK)
#define LAPACKE_cposvx_work   LAPACKE_NAME(cposvx_work,CPOSVX_WORK)
#define LAPACKE_zposvx_work   LAPACKE_NAME(zposvx_work,ZPOSVX_WORK)

#define LAPACKE_sposvxx_work   LAPACKE_NAME(sposvxx_work,SPOSVXX_WORK)
#define LAPACKE_dposvxx_work   LAPACKE_NAME(dposvxx_work,DPOSVXX_WORK)
#define LAPACKE_cposvxx_work   LAPACKE_NAME(cposvxx_work,CPOSVXX_WORK)
#define LAPACKE_zposvxx_work   LAPACKE_NAME(zposvxx_work,ZPOSVXX_WORK)

#define LAPACKE_spotrf_work   LAPACKE_NAME(spotrf_work,SPOTRF_WORK)
#define LAPACKE_dpotrf_work   LAPACKE_NAME(dpotrf_work,DPOTRF_WORK)
#define LAPACKE_cpotrf_work   LAPACKE_NAME(cpotrf_work,CPOTRF_WORK)
#define LAPACKE_zpotrf_work   LAPACKE_NAME(zpotrf_work,ZPOTRF_WORK)

#define LAPACKE_spotri_work   LAPACKE_NAME(spotri_work,SPOTRI_WORK)
#define LAPACKE_dpotri_work   LAPACKE_NAME(dpotri_work,DPOTRI_WORK)
#define LAPACKE_cpotri_work   LAPACKE_NAME(cpotri_work,CPOTRI_WORK)
#define LAPACKE_zpotri_work   LAPACKE_NAME(zpotri_work,ZPOTRI_WORK)

#define LAPACKE_spotrs_work   LAPACKE_NAME(spotrs_work,SPOTRS_WORK)
#define LAPACKE_dpotrs_work   LAPACKE_NAME(dpotrs_work,DPOTRS_WORK)
#define LAPACKE_cpotrs_work   LAPACKE_NAME(cpotrs_work,CPOTRS_WORK)
#define LAPACKE_zpotrs_work   LAPACKE_NAME(zpotrs_work,ZPOTRS_WORK)

#define LAPACKE_sppcon_work   LAPACKE_NAME(sppcon_work,SPPCON_WORK)
#define LAPACKE_dppcon_work   LAPACKE_NAME(dppcon_work,DPPCON_WORK)
#define LAPACKE_cppcon_work   LAPACKE_NAME(cppcon_work,CPPCON_WORK)
#define LAPACKE_zppcon_work   LAPACKE_NAME(zppcon_work,ZPPCON_WORK)

#define LAPACKE_sppequ_work   LAPACKE_NAME(sppequ_work,SPPEQU_WORK)
#define LAPACKE_dppequ_work   LAPACKE_NAME(dppequ_work,DPPEQU_WORK)
#define LAPACKE_cppequ_work   LAPACKE_NAME(cppequ_work,CPPEQU_WORK)
#define LAPACKE_zppequ_work   LAPACKE_NAME(zppequ_work,ZPPEQU_WORK)

#define LAPACKE_spprfs_work   LAPACKE_NAME(spprfs_work,SPPRFS_WORK)
#define LAPACKE_dpprfs_work   LAPACKE_NAME(dpprfs_work,DPPRFS_WORK)
#define LAPACKE_cpprfs_work   LAPACKE_NAME(cpprfs_work,CPPRFS_WORK)
#define LAPACKE_zpprfs_work   LAPACKE_NAME(zpprfs_work,ZPPRFS_WORK)

#define LAPACKE_sppsv_work   LAPACKE_NAME(sppsv_work,SPPSV_WORK)
#define LAPACKE_dppsv_work   LAPACKE_NAME(dppsv_work,DPPSV_WORK)
#define LAPACKE_cppsv_work   LAPACKE_NAME(cppsv_work,CPPSV_WORK)
#define LAPACKE_zppsv_work   LAPACKE_NAME(zppsv_work,ZPPSV_WORK)

#define LAPACKE_sppsvx_work   LAPACKE_NAME(sppsvx_work,SPPSVX_WORK)
#define LAPACKE_dppsvx_work   LAPACKE_NAME(dppsvx_work,DPPSVX_WORK)
#define LAPACKE_cppsvx_work   LAPACKE_NAME(cppsvx_work,CPPSVX_WORK)
#define LAPACKE_zppsvx_work   LAPACKE_NAME(zppsvx_work,ZPPSVX_WORK)

#define LAPACKE_spptrf_work   LAPACKE_NAME(spptrf_work,SPPTRF_WORK)
#define LAPACKE_dpptrf_work   LAPACKE_NAME(dpptrf_work,DPPTRF_WORK)
#define LAPACKE_cpptrf_work   LAPACKE_NAME(cpptrf_work,CPPTRF_WORK)
#define LAPACKE_zpptrf_work   LAPACKE_NAME(zpptrf_work,ZPPTRF_WORK)

#define LAPACKE_spptri_work   LAPACKE_NAME(spptri_work,SPPTRI_WORK)
#define LAPACKE_dpptri_work   LAPACKE_NAME(dpptri_work,DPPTRI_WORK)
#define LAPACKE_cpptri_work   LAPACKE_NAME(cpptri_work,CPPTRI_WORK)
#define LAPACKE_zpptri_work   LAPACKE_NAME(zpptri_work,ZPPTRI_WORK)

#define LAPACKE_spptrs_work   LAPACKE_NAME(spptrs_work,SPPTRS_WORK)
#define LAPACKE_dpptrs_work   LAPACKE_NAME(dpptrs_work,DPPTRS_WORK)
#define LAPACKE_cpptrs_work   LAPACKE_NAME(cpptrs_work,CPPTRS_WORK)
#define LAPACKE_zpptrs_work   LAPACKE_NAME(zpptrs_work,ZPPTRS_WORK)

#define LAPACKE_spstrf_work   LAPACKE_NAME(spstrf_work,SPSTRF_WORK)
#define LAPACKE_dpstrf_work   LAPACKE_NAME(dpstrf_work,DPSTRF_WORK)
#define LAPACKE_cpstrf_work   LAPACKE_NAME(cpstrf_work,CPSTRF_WORK)
#define LAPACKE_zpstrf_work   LAPACKE_NAME(zpstrf_work,ZPSTRF_WORK)

#define LAPACKE_sptcon_work   LAPACKE_NAME(sptcon_work,SPTCON_WORK)
#define LAPACKE_dptcon_work   LAPACKE_NAME(dptcon_work,DPTCON_WORK)
#define LAPACKE_cptcon_work   LAPACKE_NAME(cptcon_work,CPTCON_WORK)
#define LAPACKE_zptcon_work   LAPACKE_NAME(zptcon_work,ZPTCON_WORK)

#define LAPACKE_spteqr_work   LAPACKE_NAME(spteqr_work,SPTEQR_WORK)
#define LAPACKE_dpteqr_work   LAPACKE_NAME(dpteqr_work,DPTEQR_WORK)
#define LAPACKE_cpteqr_work   LAPACKE_NAME(cpteqr_work,CPTEQR_WORK)
#define LAPACKE_zpteqr_work   LAPACKE_NAME(zpteqr_work,ZPTEQR_WORK)

#define LAPACKE_sptrfs_work   LAPACKE_NAME(sptrfs_work,SPTRFS_WORK)
#define LAPACKE_dptrfs_work   LAPACKE_NAME(dptrfs_work,DPTRFS_WORK)
#define LAPACKE_cptrfs_work   LAPACKE_NAME(cptrfs_work,CPTRFS_WORK)
#define LAPACKE_zptrfs_work   LAPACKE_NAME(zptrfs_work,ZPTRFS_WORK)

#define LAPACKE_sptsv_work   LAPACKE_NAME(sptsv_work,SPTSV_WORK)
#define LAPACKE_dptsv_work   LAPACKE_NAME(dptsv_work,DPTSV_WORK)
#define LAPACKE_cptsv_work   LAPACKE_NAME(cptsv_work,CPTSV_WORK)
#define LAPACKE_zptsv_work   LAPACKE_NAME(zptsv_work,ZPTSV_WORK)

#define LAPACKE_sptsvx_work   LAPACKE_NAME(sptsvx_work,SPTSVX_WORK)
#define LAPACKE_dptsvx_work   LAPACKE_NAME(dptsvx_work,DPTSVX_WORK)
#define LAPACKE_cptsvx_work   LAPACKE_NAME(cptsvx_work,CPTSVX_WORK)
#define LAPACKE_zptsvx_work   LAPACKE_NAME(zptsvx_work,ZPTSVX_WORK)

#define LAPACKE_spttrf_work   LAPACKE_NAME(spttrf_work,SPTTRF_WORK)
#define LAPACKE_dpttrf_work   LAPACKE_NAME(dpttrf_work,DPTTRF_WORK)
#define LAPACKE_cpttrf_work   LAPACKE_NAME(cpttrf_work,CPTTRF_WORK)
#define LAPACKE_zpttrf_work   LAPACKE_NAME(zpttrf_work,ZPTTRF_WORK)

#define LAPACKE_spttrs_work   LAPACKE_NAME(spttrs_work,SPTTRS_WORK)
#define LAPACKE_dpttrs_work   LAPACKE_NAME(dpttrs_work,DPTTRS_WORK)
#define LAPACKE_cpttrs_work   LAPACKE_NAME(cpttrs_work,CPTTRS_WORK)
#define LAPACKE_zpttrs_work   LAPACKE_NAME(zpttrs_work,ZPTTRS_WORK)

#define LAPACKE_ssbev_work   LAPACKE_NAME(ssbev_work,SSBEV_WORK)
#define LAPACKE_dsbev_work   LAPACKE_NAME(dsbev_work,DSBEV_WORK)

#define LAPACKE_ssbevd_work   LAPACKE_NAME(ssbevd_work,SSBEVD_WORK)
#define LAPACKE_dsbevd_work   LAPACKE_NAME(dsbevd_work,DSBEVD_WORK)

#define LAPACKE_ssbevx_work   LAPACKE_NAME(ssbevx_work,SSBEVX_WORK)
#define LAPACKE_dsbevx_work   LAPACKE_NAME(dsbevx_work,DSBEVX_WORK)

#define LAPACKE_ssbgst_work   LAPACKE_NAME(ssbgst_work,SSBGST_WORK)
#define LAPACKE_dsbgst_work   LAPACKE_NAME(dsbgst_work,DSBGST_WORK)

#define LAPACKE_ssbgv_work   LAPACKE_NAME(ssbgv_work,SSBGV_WORK)
#define LAPACKE_dsbgv_work   LAPACKE_NAME(dsbgv_work,DSBGV_WORK)

#define LAPACKE_ssbgvd_work   LAPACKE_NAME(ssbgvd_work,SSBGVD_WORK)
#define LAPACKE_dsbgvd_work   LAPACKE_NAME(dsbgvd_work,DSBGVD_WORK)

#define LAPACKE_ssbgvx_work   LAPACKE_NAME(ssbgvx_work,SSBGVX_WORK)
#define LAPACKE_dsbgvx_work   LAPACKE_NAME(dsbgvx_work,DSBGVX_WORK)

#define LAPACKE_ssbtrd_work   LAPACKE_NAME(ssbtrd_work,SSBTRD_WORK)
#define LAPACKE_dsbtrd_work   LAPACKE_NAME(dsbtrd_work,DSBTRD_WORK)

#define LAPACKE_ssfrk_work   LAPACKE_NAME(ssfrk_work,SSFRK_WORK)
#define LAPACKE_dsfrk_work   LAPACKE_NAME(dsfrk_work,DSFRK_WORK)

#define LAPACKE_sspcon_work   LAPACKE_NAME(sspcon_work,SSPCON_WORK)
#define LAPACKE_dspcon_work   LAPACKE_NAME(dspcon_work,DSPCON_WORK)
#define LAPACKE_cspcon_work   LAPACKE_NAME(cspcon_work,CSPCON_WORK)
#define LAPACKE_zspcon_work   LAPACKE_NAME(zspcon_work,ZSPCON_WORK)

#define LAPACKE_sspev_work   LAPACKE_NAME(sspev_work,SSPEV_WORK)
#define LAPACKE_dspev_work   LAPACKE_NAME(dspev_work,DSPEV_WORK)

#define LAPACKE_sspevd_work   LAPACKE_NAME(sspevd_work,SSPEVD_WORK)
#define LAPACKE_dspevd_work   LAPACKE_NAME(dspevd_work,DSPEVD_WORK)

#define LAPACKE_sspevx_work   LAPACKE_NAME(sspevx_work,SSPEVX_WORK)
#define LAPACKE_dspevx_work   LAPACKE_NAME(dspevx_work,DSPEVX_WORK)

#define LAPACKE_sspgst_work   LAPACKE_NAME(sspgst_work,SSPGST_WORK)
#define LAPACKE_dspgst_work   LAPACKE_NAME(dspgst_work,DSPGST_WORK)

#define LAPACKE_sspgv_work   LAPACKE_NAME(sspgv_work,SSPGV_WORK)
#define LAPACKE_dspgv_work   LAPACKE_NAME(dspgv_work,DSPGV_WORK)

#define LAPACKE_sspgvd_work   LAPACKE_NAME(sspgvd_work,SSPGVD_WORK)
#define LAPACKE_dspgvd_work   LAPACKE_NAME(dspgvd_work,DSPGVD_WORK)

#define LAPACKE_sspgvx_work   LAPACKE_NAME(sspgvx_work,SSPGVX_WORK)
#define LAPACKE_dspgvx_work   LAPACKE_NAME(dspgvx_work,DSPGVX_WORK)

#define LAPACKE_ssprfs_work   LAPACKE_NAME(ssprfs_work,SSPRFS_WORK)
#define LAPACKE_dsprfs_work   LAPACKE_NAME(dsprfs_work,DSPRFS_WORK)
#define LAPACKE_csprfs_work   LAPACKE_NAME(csprfs_work,CSPRFS_WORK)
#define LAPACKE_zsprfs_work   LAPACKE_NAME(zsprfs_work,ZSPRFS_WORK)

#define LAPACKE_sspsv_work   LAPACKE_NAME(sspsv_work,SSPSV_WORK)
#define LAPACKE_dspsv_work   LAPACKE_NAME(dspsv_work,DSPSV_WORK)
#define LAPACKE_cspsv_work   LAPACKE_NAME(cspsv_work,CSPSV_WORK)
#define LAPACKE_zspsv_work   LAPACKE_NAME(zspsv_work,ZSPSV_WORK)

#define LAPACKE_sspsvx_work   LAPACKE_NAME(sspsvx_work,SSPSVX_WORK)
#define LAPACKE_dspsvx_work   LAPACKE_NAME(dspsvx_work,DSPSVX_WORK)
#define LAPACKE_cspsvx_work   LAPACKE_NAME(cspsvx_work,CSPSVX_WORK)
#define LAPACKE_zspsvx_work   LAPACKE_NAME(zspsvx_work,ZSPSVX_WORK)

#define LAPACKE_ssptrd_work   LAPACKE_NAME(ssptrd_work,SSPTRD_WORK)
#define LAPACKE_dsptrd_work   LAPACKE_NAME(dsptrd_work,DSPTRD_WORK)

#define LAPACKE_ssptrf_work   LAPACKE_NAME(ssptrf_work,SSPTRF_WORK)
#define LAPACKE_dsptrf_work   LAPACKE_NAME(dsptrf_work,DSPTRF_WORK)
#define LAPACKE_csptrf_work   LAPACKE_NAME(csptrf_work,CSPTRF_WORK)
#define LAPACKE_zsptrf_work   LAPACKE_NAME(zsptrf_work,ZSPTRF_WORK)

#define LAPACKE_ssptri_work   LAPACKE_NAME(ssptri_work,SSPTRI_WORK)
#define LAPACKE_dsptri_work   LAPACKE_NAME(dsptri_work,DSPTRI_WORK)
#define LAPACKE_csptri_work   LAPACKE_NAME(csptri_work,CSPTRI_WORK)
#define LAPACKE_zsptri_work   LAPACKE_NAME(zsptri_work,ZSPTRI_WORK)

#define LAPACKE_ssptrs_work   LAPACKE_NAME(ssptrs_work,SSPTRS_WORK)
#define LAPACKE_dsptrs_work   LAPACKE_NAME(dsptrs_work,DSPTRS_WORK)
#define LAPACKE_csptrs_work   LAPACKE_NAME(csptrs_work,CSPTRS_WORK)
#define LAPACKE_zsptrs_work   LAPACKE_NAME(zsptrs_work,ZSPTRS_WORK)

#define LAPACKE_sstebz_work   LAPACKE_NAME(sstebz_work,SSTEBZ_WORK)
#define LAPACKE_dstebz_work   LAPACKE_NAME(dstebz_work,DSTEBZ_WORK)

#define LAPACKE_sstedc_work   LAPACKE_NAME(sstedc_work,SSTEDC_WORK)
#define LAPACKE_dstedc_work   LAPACKE_NAME(dstedc_work,DSTEDC_WORK)
#define LAPACKE_cstedc_work   LAPACKE_NAME(cstedc_work,CSTEDC_WORK)
#define LAPACKE_zstedc_work   LAPACKE_NAME(zstedc_work,ZSTEDC_WORK)

#define LAPACKE_sstegr_work   LAPACKE_NAME(sstegr_work,SSTEGR_WORK)
#define LAPACKE_dstegr_work   LAPACKE_NAME(dstegr_work,DSTEGR_WORK)
#define LAPACKE_cstegr_work   LAPACKE_NAME(cstegr_work,CSTEGR_WORK)
#define LAPACKE_zstegr_work   LAPACKE_NAME(zstegr_work,ZSTEGR_WORK)

#define LAPACKE_sstein_work   LAPACKE_NAME(sstein_work,SSTEIN_WORK)
#define LAPACKE_dstein_work   LAPACKE_NAME(dstein_work,DSTEIN_WORK)
#define LAPACKE_cstein_work   LAPACKE_NAME(cstein_work,CSTEIN_WORK)
#define LAPACKE_zstein_work   LAPACKE_NAME(zstein_work,ZSTEIN_WORK)

#define LAPACKE_sstemr_work   LAPACKE_NAME(sstemr_work,SSTEMR_WORK)
#define LAPACKE_dstemr_work   LAPACKE_NAME(dstemr_work,DSTEMR_WORK)
#define LAPACKE_cstemr_work   LAPACKE_NAME(cstemr_work,CSTEMR_WORK)
#define LAPACKE_zstemr_work   LAPACKE_NAME(zstemr_work,ZSTEMR_WORK)

#define LAPACKE_ssteqr_work   LAPACKE_NAME(ssteqr_work,SSTEQR_WORK)
#define LAPACKE_dsteqr_work   LAPACKE_NAME(dsteqr_work,DSTEQR_WORK)
#define LAPACKE_csteqr_work   LAPACKE_NAME(csteqr_work,CSTEQR_WORK)
#define LAPACKE_zsteqr_work   LAPACKE_NAME(zsteqr_work,ZSTEQR_WORK)

#define LAPACKE_ssterf_work   LAPACKE_NAME(ssterf_work,SSTERF_WORK)
#define LAPACKE_dsterf_work   LAPACKE_NAME(dsterf_work,DSTERF_WORK)

#define LAPACKE_sstev_work   LAPACKE_NAME(sstev_work,SSTEV_WORK)
#define LAPACKE_dstev_work   LAPACKE_NAME(dstev_work,DSTEV_WORK)

#define LAPACKE_sstevd_work   LAPACKE_NAME(sstevd_work,SSTEVD_WORK)
#define LAPACKE_dstevd_work   LAPACKE_NAME(dstevd_work,DSTEVD_WORK)

#define LAPACKE_sstevr_work   LAPACKE_NAME(sstevr_work,SSTEVR_WORK)
#define LAPACKE_dstevr_work   LAPACKE_NAME(dstevr_work,DSTEVR_WORK)

#define LAPACKE_sstevx_work   LAPACKE_NAME(sstevx_work,SSTEVX_WORK)
#define LAPACKE_dstevx_work   LAPACKE_NAME(dstevx_work,DSTEVX_WORK)

#define LAPACKE_ssycon_work   LAPACKE_NAME(ssycon_work,SSYCON_WORK)
#define LAPACKE_dsycon_work   LAPACKE_NAME(dsycon_work,DSYCON_WORK)
#define LAPACKE_csycon_work   LAPACKE_NAME(csycon_work,CSYCON_WORK)
#define LAPACKE_zsycon_work   LAPACKE_NAME(zsycon_work,ZSYCON_WORK)

#define LAPACKE_ssyequb_work   LAPACKE_NAME(ssyequb_work,SSYEQUB_WORK)
#define LAPACKE_dsyequb_work   LAPACKE_NAME(dsyequb_work,DSYEQUB_WORK)
#define LAPACKE_csyequb_work   LAPACKE_NAME(csyequb_work,CSYEQUB_WORK)
#define LAPACKE_zsyequb_work   LAPACKE_NAME(zsyequb_work,ZSYEQUB_WORK)

#define LAPACKE_ssyev_work   LAPACKE_NAME(ssyev_work,SSYEV_WORK)
#define LAPACKE_dsyev_work   LAPACKE_NAME(dsyev_work,DSYEV_WORK)

#define LAPACKE_ssyevd_work   LAPACKE_NAME(ssyevd_work,SSYEVD_WORK)
#define LAPACKE_dsyevd_work   LAPACKE_NAME(dsyevd_work,DSYEVD_WORK)

#define LAPACKE_ssyevr_work   LAPACKE_NAME(ssyevr_work,SSYEVR_WORK)
#define LAPACKE_dsyevr_work   LAPACKE_NAME(dsyevr_work,DSYEVR_WORK)

#define LAPACKE_ssyevx_work   LAPACKE_NAME(ssyevx_work,SSYEVX_WORK)
#define LAPACKE_dsyevx_work   LAPACKE_NAME(dsyevx_work,DSYEVX_WORK)

#define LAPACKE_ssygst_work   LAPACKE_NAME(ssygst_work,SSYGST_WORK)
#define LAPACKE_dsygst_work   LAPACKE_NAME(dsygst_work,DSYGST_WORK)

#define LAPACKE_ssygv_work   LAPACKE_NAME(ssygv_work,SSYGV_WORK)
#define LAPACKE_dsygv_work   LAPACKE_NAME(dsygv_work,DSYGV_WORK)

#define LAPACKE_ssygvd_work   LAPACKE_NAME(ssygvd_work,SSYGVD_WORK)
#define LAPACKE_dsygvd_work   LAPACKE_NAME(dsygvd_work,DSYGVD_WORK)

#define LAPACKE_ssygvx_work   LAPACKE_NAME(ssygvx_work,SSYGVX_WORK)
#define LAPACKE_dsygvx_work   LAPACKE_NAME(dsygvx_work,DSYGVX_WORK)

#define LAPACKE_ssyrfs_work   LAPACKE_NAME(ssyrfs_work,SSYRFS_WORK)
#define LAPACKE_dsyrfs_work   LAPACKE_NAME(dsyrfs_work,DSYRFS_WORK)
#define LAPACKE_csyrfs_work   LAPACKE_NAME(csyrfs_work,CSYRFS_WORK)
#define LAPACKE_zsyrfs_work   LAPACKE_NAME(zsyrfs_work,ZSYRFS_WORK)

#define LAPACKE_ssyrfsx_work   LAPACKE_NAME(ssyrfsx_work,SSYRFSX_WORK)
#define LAPACKE_dsyrfsx_work   LAPACKE_NAME(dsyrfsx_work,DSYRFSX_WORK)
#define LAPACKE_csyrfsx_work   LAPACKE_NAME(csyrfsx_work,CSYRFSX_WORK)
#define LAPACKE_zsyrfsx_work   LAPACKE_NAME(zsyrfsx_work,ZSYRFSX_WORK)

#define LAPACKE_ssysv_work   LAPACKE_NAME(ssysv_work,SSYSV_WORK)
#define LAPACKE_dsysv_work   LAPACKE_NAME(dsysv_work,DSYSV_WORK)
#define LAPACKE_csysv_work   LAPACKE_NAME(csysv_work,CSYSV_WORK)
#define LAPACKE_zsysv_work   LAPACKE_NAME(zsysv_work,ZSYSV_WORK)

#define LAPACKE_ssysvx_work   LAPACKE_NAME(ssysvx_work,SSYSVX_WORK)
#define LAPACKE_dsysvx_work   LAPACKE_NAME(dsysvx_work,DSYSVX_WORK)
#define LAPACKE_csysvx_work   LAPACKE_NAME(csysvx_work,CSYSVX_WORK)
#define LAPACKE_zsysvx_work   LAPACKE_NAME(zsysvx_work,ZSYSVX_WORK)

#define LAPACKE_ssysvxx_work   LAPACKE_NAME(ssysvxx_work,SSYSVXX_WORK)
#define LAPACKE_dsysvxx_work   LAPACKE_NAME(dsysvxx_work,DSYSVXX_WORK)
#define LAPACKE_csysvxx_work   LAPACKE_NAME(csysvxx_work,CSYSVXX_WORK)
#define LAPACKE_zsysvxx_work   LAPACKE_NAME(zsysvxx_work,ZSYSVXX_WORK)

#define LAPACKE_ssytrd_work   LAPACKE_NAME(ssytrd_work,SSYTRD_WORK)
#define LAPACKE_dsytrd_work   LAPACKE_NAME(dsytrd_work,DSYTRD_WORK)

#define LAPACKE_ssytrf_work   LAPACKE_NAME(ssytrf_work,SSYTRF_WORK)
#define LAPACKE_dsytrf_work   LAPACKE_NAME(dsytrf_work,DSYTRF_WORK)
#define LAPACKE_csytrf_work   LAPACKE_NAME(csytrf_work,CSYTRF_WORK)
#define LAPACKE_zsytrf_work   LAPACKE_NAME(zsytrf_work,ZSYTRF_WORK)

#define LAPACKE_ssytri_work   LAPACKE_NAME(ssytri_work,SSYTRI_WORK)
#define LAPACKE_dsytri_work   LAPACKE_NAME(dsytri_work,DSYTRI_WORK)
#define LAPACKE_csytri_work   LAPACKE_NAME(csytri_work,CSYTRI_WORK)
#define LAPACKE_zsytri_work   LAPACKE_NAME(zsytri_work,ZSYTRI_WORK)

#define LAPACKE_ssytrs_work   LAPACKE_NAME(ssytrs_work,SSYTRS_WORK)
#define LAPACKE_dsytrs_work   LAPACKE_NAME(dsytrs_work,DSYTRS_WORK)
#define LAPACKE_csytrs_work   LAPACKE_NAME(csytrs_work,CSYTRS_WORK)
#define LAPACKE_zsytrs_work   LAPACKE_NAME(zsytrs_work,ZSYTRS_WORK)

#define LAPACKE_stbcon_work   LAPACKE_NAME(stbcon_work,STBCON_WORK)
#define LAPACKE_dtbcon_work   LAPACKE_NAME(dtbcon_work,DTBCON_WORK)
#define LAPACKE_ctbcon_work   LAPACKE_NAME(ctbcon_work,CTBCON_WORK)
#define LAPACKE_ztbcon_work   LAPACKE_NAME(ztbcon_work,ZTBCON_WORK)

#define LAPACKE_stbrfs_work   LAPACKE_NAME(stbrfs_work,STBRFS_WORK)
#define LAPACKE_dtbrfs_work   LAPACKE_NAME(dtbrfs_work,DTBRFS_WORK)
#define LAPACKE_ctbrfs_work   LAPACKE_NAME(ctbrfs_work,CTBRFS_WORK)
#define LAPACKE_ztbrfs_work   LAPACKE_NAME(ztbrfs_work,ZTBRFS_WORK)

#define LAPACKE_stbtrs_work   LAPACKE_NAME(stbtrs_work,STBTRS_WORK)
#define LAPACKE_dtbtrs_work   LAPACKE_NAME(dtbtrs_work,DTBTRS_WORK)
#define LAPACKE_ctbtrs_work   LAPACKE_NAME(ctbtrs_work,CTBTRS_WORK)
#define LAPACKE_ztbtrs_work   LAPACKE_NAME(ztbtrs_work,ZTBTRS_WORK)

#define LAPACKE_stfsm_work   LAPACKE_NAME(stfsm_work,STFSM_WORK)
#define LAPACKE_dtfsm_work   LAPACKE_NAME(dtfsm_work,DTFSM_WORK)
#define LAPACKE_ctfsm_work   LAPACKE_NAME(ctfsm_work,CTFSM_WORK)
#define LAPACKE_ztfsm_work   LAPACKE_NAME(ztfsm_work,ZTFSM_WORK)

#define LAPACKE_stftri_work   LAPACKE_NAME(stftri_work,STFTRI_WORK)
#define LAPACKE_dtftri_work   LAPACKE_NAME(dtftri_work,DTFTRI_WORK)
#define LAPACKE_ctftri_work   LAPACKE_NAME(ctftri_work,CTFTRI_WORK)
#define LAPACKE_ztftri_work   LAPACKE_NAME(ztftri_work,ZTFTRI_WORK)

#define LAPACKE_stfttp_work   LAPACKE_NAME(stfttp_work,STFTTP_WORK)
#define LAPACKE_dtfttp_work   LAPACKE_NAME(dtfttp_work,DTFTTP_WORK)
#define LAPACKE_ctfttp_work   LAPACKE_NAME(ctfttp_work,CTFTTP_WORK)
#define LAPACKE_ztfttp_work   LAPACKE_NAME(ztfttp_work,ZTFTTP_WORK)

#define LAPACKE_stfttr_work   LAPACKE_NAME(stfttr_work,STFTTR_WORK)
#define LAPACKE_dtfttr_work   LAPACKE_NAME(dtfttr_work,DTFTTR_WORK)
#define LAPACKE_ctfttr_work   LAPACKE_NAME(ctfttr_work,CTFTTR_WORK)
#define LAPACKE_ztfttr_work   LAPACKE_NAME(ztfttr_work,ZTFTTR_WORK)

#define LAPACKE_stgevc_work   LAPACKE_NAME(stgevc_work,STGEVC_WORK)
#define LAPACKE_dtgevc_work   LAPACKE_NAME(dtgevc_work,DTGEVC_WORK)
#define LAPACKE_ctgevc_work   LAPACKE_NAME(ctgevc_work,CTGEVC_WORK)
#define LAPACKE_ztgevc_work   LAPACKE_NAME(ztgevc_work,ZTGEVC_WORK)

#define LAPACKE_stgexc_work   LAPACKE_NAME(stgexc_work,STGEXC_WORK)
#define LAPACKE_dtgexc_work   LAPACKE_NAME(dtgexc_work,DTGEXC_WORK)
#define LAPACKE_ctgexc_work   LAPACKE_NAME(ctgexc_work,CTGEXC_WORK)
#define LAPACKE_ztgexc_work   LAPACKE_NAME(ztgexc_work,ZTGEXC_WORK)

#define LAPACKE_stgsen_work   LAPACKE_NAME(stgsen_work,STGSEN_WORK)
#define LAPACKE_dtgsen_work   LAPACKE_NAME(dtgsen_work,DTGSEN_WORK)
#define LAPACKE_ctgsen_work   LAPACKE_NAME(ctgsen_work,CTGSEN_WORK)
#define LAPACKE_ztgsen_work   LAPACKE_NAME(ztgsen_work,ZTGSEN_WORK)

#define LAPACKE_stgsja_work   LAPACKE_NAME(stgsja_work,STGSJA_WORK)
#define LAPACKE_dtgsja_work   LAPACKE_NAME(dtgsja_work,DTGSJA_WORK)
#define LAPACKE_ctgsja_work   LAPACKE_NAME(ctgsja_work,CTGSJA_WORK)
#define LAPACKE_ztgsja_work   LAPACKE_NAME(ztgsja_work,ZTGSJA_WORK)

#define LAPACKE_stgsna_work   LAPACKE_NAME(stgsna_work,STGSNA_WORK)
#define LAPACKE_dtgsna_work   LAPACKE_NAME(dtgsna_work,DTGSNA_WORK)
#define LAPACKE_ctgsna_work   LAPACKE_NAME(ctgsna_work,CTGSNA_WORK)
#define LAPACKE_ztgsna_work   LAPACKE_NAME(ztgsna_work,ZTGSNA_WORK)

#define LAPACKE_stgsyl_work   LAPACKE_NAME(stgsyl_work,STGSYL_WORK)
#define LAPACKE_dtgsyl_work   LAPACKE_NAME(dtgsyl_work,DTGSYL_WORK)
#define LAPACKE_ctgsyl_work   LAPACKE_NAME(ctgsyl_work,CTGSYL_WORK)
#define LAPACKE_ztgsyl_work   LAPACKE_NAME(ztgsyl_work,ZTGSYL_WORK)

#define LAPACKE_stpcon_work   LAPACKE_NAME(stpcon_work,STPCON_WORK)
#define LAPACKE_dtpcon_work   LAPACKE_NAME(dtpcon_work,DTPCON_WORK)
#define LAPACKE_ctpcon_work   LAPACKE_NAME(ctpcon_work,CTPCON_WORK)
#define LAPACKE_ztpcon_work   LAPACKE_NAME(ztpcon_work,ZTPCON_WORK)

#define LAPACKE_stprfs_work   LAPACKE_NAME(stprfs_work,STPRFS_WORK)
#define LAPACKE_dtprfs_work   LAPACKE_NAME(dtprfs_work,DTPRFS_WORK)
#define LAPACKE_ctprfs_work   LAPACKE_NAME(ctprfs_work,CTPRFS_WORK)
#define LAPACKE_ztprfs_work   LAPACKE_NAME(ztprfs_work,ZTPRFS_WORK)

#define LAPACKE_stptri_work   LAPACKE_NAME(stptri_work,STPTRI_WORK)
#define LAPACKE_dtptri_work   LAPACKE_NAME(dtptri_work,DTPTRI_WORK)
#define LAPACKE_ctptri_work   LAPACKE_NAME(ctptri_work,CTPTRI_WORK)
#define LAPACKE_ztptri_work   LAPACKE_NAME(ztptri_work,ZTPTRI_WORK)

#define LAPACKE_stptrs_work   LAPACKE_NAME(stptrs_work,STPTRS_WORK)
#define LAPACKE_dtptrs_work   LAPACKE_NAME(dtptrs_work,DTPTRS_WORK)
#define LAPACKE_ctptrs_work   LAPACKE_NAME(ctptrs_work,CTPTRS_WORK)
#define LAPACKE_ztptrs_work   LAPACKE_NAME(ztptrs_work,ZTPTRS_WORK)

#define LAPACKE_stpttf_work   LAPACKE_NAME(stpttf_work,STPTTF_WORK)
#define LAPACKE_dtpttf_work   LAPACKE_NAME(dtpttf_work,DTPTTF_WORK)
#define LAPACKE_ctpttf_work   LAPACKE_NAME(ctpttf_work,CTPTTF_WORK)
#define LAPACKE_ztpttf_work   LAPACKE_NAME(ztpttf_work,ZTPTTF_WORK)

#define LAPACKE_stpttr_work   LAPACKE_NAME(stpttr_work,STPTTR_WORK)
#define LAPACKE_dtpttr_work   LAPACKE_NAME(dtpttr_work,DTPTTR_WORK)
#define LAPACKE_ctpttr_work   LAPACKE_NAME(ctpttr_work,CTPTTR_WORK)
#define LAPACKE_ztpttr_work   LAPACKE_NAME(ztpttr_work,ZTPTTR_WORK)

#define LAPACKE_strcon_work   LAPACKE_NAME(strcon_work,STRCON_WORK)
#define LAPACKE_dtrcon_work   LAPACKE_NAME(dtrcon_work,DTRCON_WORK)
#define LAPACKE_ctrcon_work   LAPACKE_NAME(ctrcon_work,CTRCON_WORK)
#define LAPACKE_ztrcon_work   LAPACKE_NAME(ztrcon_work,ZTRCON_WORK)

#define LAPACKE_strevc_work   LAPACKE_NAME(strevc_work,STREVC_WORK)
#define LAPACKE_dtrevc_work   LAPACKE_NAME(dtrevc_work,DTREVC_WORK)
#define LAPACKE_ctrevc_work   LAPACKE_NAME(ctrevc_work,CTREVC_WORK)
#define LAPACKE_ztrevc_work   LAPACKE_NAME(ztrevc_work,ZTREVC_WORK)

#define LAPACKE_strexc_work   LAPACKE_NAME(strexc_work,STREXC_WORK)
#define LAPACKE_dtrexc_work   LAPACKE_NAME(dtrexc_work,DTREXC_WORK)
#define LAPACKE_ctrexc_work   LAPACKE_NAME(ctrexc_work,CTREXC_WORK)
#define LAPACKE_ztrexc_work   LAPACKE_NAME(ztrexc_work,ZTREXC_WORK)

#define LAPACKE_strrfs_work   LAPACKE_NAME(strrfs_work,STRRFS_WORK)
#define LAPACKE_dtrrfs_work   LAPACKE_NAME(dtrrfs_work,DTRRFS_WORK)
#define LAPACKE_ctrrfs_work   LAPACKE_NAME(ctrrfs_work,CTRRFS_WORK)
#define LAPACKE_ztrrfs_work   LAPACKE_NAME(ztrrfs_work,ZTRRFS_WORK)

#define LAPACKE_strsen_work   LAPACKE_NAME(strsen_work,STRSEN_WORK)
#define LAPACKE_dtrsen_work   LAPACKE_NAME(dtrsen_work,DTRSEN_WORK)
#define LAPACKE_ctrsen_work   LAPACKE_NAME(ctrsen_work,CTRSEN_WORK)
#define LAPACKE_ztrsen_work   LAPACKE_NAME(ztrsen_work,ZTRSEN_WORK)

#define LAPACKE_strsna_work   LAPACKE_NAME(strsna_work,STRSNA_WORK)
#define LAPACKE_dtrsna_work   LAPACKE_NAME(dtrsna_work,DTRSNA_WORK)
#define LAPACKE_ctrsna_work   LAPACKE_NAME(ctrsna_work,CTRSNA_WORK)
#define LAPACKE_ztrsna_work   LAPACKE_NAME(ztrsna_work,ZTRSNA_WORK)

#define LAPACKE_strsyl_work   LAPACKE_NAME(strsyl_work,STRSYL_WORK)
#define LAPACKE_dtrsyl_work   LAPACKE_NAME(dtrsyl_work,DTRSYL_WORK)
#define LAPACKE_ctrsyl_work   LAPACKE_NAME(ctrsyl_work,CTRSYL_WORK)
#define LAPACKE_ztrsyl_work   LAPACKE_NAME(ztrsyl_work,ZTRSYL_WORK)

#define LAPACKE_strtri_work   LAPACKE_NAME(strtri_work,STRTRI_WORK)
#define LAPACKE_dtrtri_work   LAPACKE_NAME(dtrtri_work,DTRTRI_WORK)
#define LAPACKE_ctrtri_work   LAPACKE_NAME(ctrtri_work,CTRTRI_WORK)
#define LAPACKE_ztrtri_work   LAPACKE_NAME(ztrtri_work,ZTRTRI_WORK)

#define LAPACKE_strtrs_work   LAPACKE_NAME(strtrs_work,STRTRS_WORK)
#define LAPACKE_dtrtrs_work   LAPACKE_NAME(dtrtrs_work,DTRTRS_WORK)
#define LAPACKE_ctrtrs_work   LAPACKE_NAME(ctrtrs_work,CTRTRS_WORK)
#define LAPACKE_ztrtrs_work   LAPACKE_NAME(ztrtrs_work,ZTRTRS_WORK)

#define LAPACKE_strttf_work   LAPACKE_NAME(strttf_work,STRTTF_WORK)
#define LAPACKE_dtrttf_work   LAPACKE_NAME(dtrttf_work,DTRTTF_WORK)
#define LAPACKE_ctrttf_work   LAPACKE_NAME(ctrttf_work,CTRTTF_WORK)
#define LAPACKE_ztrttf_work   LAPACKE_NAME(ztrttf_work,ZTRTTF_WORK)

#define LAPACKE_strttp_work   LAPACKE_NAME(strttp_work,STRTTP_WORK)
#define LAPACKE_dtrttp_work   LAPACKE_NAME(dtrttp_work,DTRTTP_WORK)
#define LAPACKE_ctrttp_work   LAPACKE_NAME(ctrttp_work,CTRTTP_WORK)
#define LAPACKE_ztrttp_work   LAPACKE_NAME(ztrttp_work,ZTRTTP_WORK)

#define LAPACKE_stzrzf_work   LAPACKE_NAME(stzrzf_work,STZRZF_WORK)
#define LAPACKE_dtzrzf_work   LAPACKE_NAME(dtzrzf_work,DTZRZF_WORK)
#define LAPACKE_ctzrzf_work   LAPACKE_NAME(ctzrzf_work,CTZRZF_WORK)
#define LAPACKE_ztzrzf_work   LAPACKE_NAME(ztzrzf_work,ZTZRZF_WORK)

#define LAPACKE_cungbr_work   LAPACKE_NAME(cungbr_work,CUNGBR_WORK)
#define LAPACKE_zungbr_work   LAPACKE_NAME(zungbr_work,ZUNGBR_WORK)

#define LAPACKE_cunghr_work   LAPACKE_NAME(cunghr_work,CUNGHR_WORK)
#define LAPACKE_zunghr_work   LAPACKE_NAME(zunghr_work,ZUNGHR_WORK)

#define LAPACKE_cunglq_work   LAPACKE_NAME(cunglq_work,CUNGLQ_WORK)
#define LAPACKE_zunglq_work   LAPACKE_NAME(zunglq_work,ZUNGLQ_WORK)

#define LAPACKE_cungql_work   LAPACKE_NAME(cungql_work,CUNGQL_WORK)
#define LAPACKE_zungql_work   LAPACKE_NAME(zungql_work,ZUNGQL_WORK)

#define LAPACKE_cungqr_work   LAPACKE_NAME(cungqr_work,CUNGQR_WORK)
#define LAPACKE_zungqr_work   LAPACKE_NAME(zungqr_work,ZUNGQR_WORK)

#define LAPACKE_cungrq_work   LAPACKE_NAME(cungrq_work,CUNGRQ_WORK)
#define LAPACKE_zungrq_work   LAPACKE_NAME(zungrq_work,ZUNGRQ_WORK)

#define LAPACKE_cungtr_work   LAPACKE_NAME(cungtr_work,CUNGTR_WORK)
#define LAPACKE_zungtr_work   LAPACKE_NAME(zungtr_work,ZUNGTR_WORK)

#define LAPACKE_cunmbr_work   LAPACKE_NAME(cunmbr_work,CUNMBR_WORK)
#define LAPACKE_zunmbr_work   LAPACKE_NAME(zunmbr_work,ZUNMBR_WORK)

#define LAPACKE_cunmhr_work   LAPACKE_NAME(cunmhr_work,CUNMHR_WORK)
#define LAPACKE_zunmhr_work   LAPACKE_NAME(zunmhr_work,ZUNMHR_WORK)

#define LAPACKE_cunmlq_work   LAPACKE_NAME(cunmlq_work,CUNMLQ_WORK)
#define LAPACKE_zunmlq_work   LAPACKE_NAME(zunmlq_work,ZUNMLQ_WORK)

#define LAPACKE_cunmql_work   LAPACKE_NAME(cunmql_work,CUNMQL_WORK)
#define LAPACKE_zunmql_work   LAPACKE_NAME(zunmql_work,ZUNMQL_WORK)

#define LAPACKE_cunmqr_work   LAPACKE_NAME(cunmqr_work,CUNMQR_WORK)
#define LAPACKE_zunmqr_work   LAPACKE_NAME(zunmqr_work,ZUNMQR_WORK)

#define LAPACKE_cunmrq_work   LAPACKE_NAME(cunmrq_work,CUNMRQ_WORK)
#define LAPACKE_zunmrq_work   LAPACKE_NAME(zunmrq_work,ZUNMRQ_WORK)

#define LAPACKE_cunmrz_work   LAPACKE_NAME(cunmrz_work,CUNMRZ_WORK)
#define LAPACKE_zunmrz_work   LAPACKE_NAME(zunmrz_work,ZUNMRZ_WORK)

#define LAPACKE_cunmtr_work   LAPACKE_NAME(cunmtr_work,CUNMTR_WORK)
#define LAPACKE_zunmtr_work   LAPACKE_NAME(zunmtr_work,ZUNMTR_WORK)

#define LAPACKE_cupgtr_work   LAPACKE_NAME(cupgtr_work,CUPGTR_WORK)
#define LAPACKE_zupgtr_work   LAPACKE_NAME(zupgtr_work,ZUPGTR_WORK)

#define LAPACKE_cupmtr_work   LAPACKE_NAME(cupmtr_work,CUPMTR_WORK)
#define LAPACKE_zupmtr_work   LAPACKE_NAME(zupmtr_work,ZUPMTR_WORK)

#define LAPACKE_claghe   LAPACKE_NAME(claghe,CLAGHE)
#define LAPACKE_zlaghe   LAPACKE_NAME(zlaghe,ZLAGHE)

#define LAPACKE_slagsy   LAPACKE_NAME(slagsy,SLAGSY)
#define LAPACKE_dlagsy   LAPACKE_NAME(dlagsy,DLAGSY)
#define LAPACKE_clagsy   LAPACKE_NAME(clagsy,CLAGSY)
#define LAPACKE_zlagsy   LAPACKE_NAME(zlagsy,ZLAGSY)

#define LAPACKE_slapmr   LAPACKE_NAME(slapmr,SLAPMR)
#define LAPACKE_dlapmr   LAPACKE_NAME(dlapmr,DLAPMR)
#define LAPACKE_clapmr   LAPACKE_NAME(clapmr,CLAPMR)
#define LAPACKE_zlapmr   LAPACKE_NAME(zlapmr,ZLAPMR)

#define LAPACKE_slapy2   LAPACKE_NAME(slapy2,SLAPY2)
#define LAPACKE_dlapy2   LAPACKE_NAME(dlapy2,DLAPY2)

#define LAPACKE_slapy3   LAPACKE_NAME(slapy3,SLAPY3)
#define LAPACKE_dlapy3   LAPACKE_NAME(dlapy3,DLAPY3)

#define LAPACKE_slartgp   LAPACKE_NAME(slartgp,SLARTGP)
#define LAPACKE_dlartgp   LAPACKE_NAME(dlartgp,DLARTGP)

#define LAPACKE_slartgs   LAPACKE_NAME(slartgs,SLARTGS)
#define LAPACKE_dlartgs   LAPACKE_NAME(dlartgs,DLARTGS)

//LAPACK 3.3.0
#define LAPACKE_cbbcsd_work LAPACKE_NAME(cbbcsd_work,CBBCSD_WORK)
#define LAPACKE_cheswapr_work LAPACKE_NAME(cheswapr_work,CHESWAPR_WORK)
#define LAPACKE_chetri2_work LAPACKE_NAME(chetri2_work,CHETRI2_WORK)
#define LAPACKE_chetri2x_work LAPACKE_NAME(chetri2x_work,CHETRI2X_WORK)
#define LAPACKE_chetrs2_work LAPACKE_NAME(chetrs2_work,CHETRS2_WORK)
#define LAPACKE_csyconv_work LAPACKE_NAME(csyconv_work,CSYCONV_WORK)
#define LAPACKE_csyswapr_work LAPACKE_NAME(csyswapr_work,CSYSWAPR_WORK)
#define LAPACKE_csytri2_work LAPACKE_NAME(csytri2_work,CSYTRI2_WORK)
#define LAPACKE_csytri2x_work LAPACKE_NAME(csytri2x_work,CSYTRI2X_WORK)
#define LAPACKE_csytrs2_work LAPACKE_NAME(csytrs2_work,CSYTRS2_WORK)
#define LAPACKE_cunbdb_work LAPACKE_NAME(cunbdb_work,CUNBDB_WORK)
#define LAPACKE_cuncsd_work LAPACKE_NAME(cuncsd_work,CUNCSD_WORK)
#define LAPACKE_dbbcsd_work LAPACKE_NAME(dbbcsd_work,DBBCSD_WORK)
#define LAPACKE_dorbdb_work LAPACKE_NAME(dorbdb_work,DORBDB_WORK)
#define LAPACKE_dorcsd_work LAPACKE_NAME(dorcsd_work,DORCSD_WORK)
#define LAPACKE_dsyconv_work LAPACKE_NAME(dsyconv_work,DSYCONV_WORK)
#define LAPACKE_dsyswapr_work LAPACKE_NAME(dsyswapr_work,DSYSWAPR_WORK)
#define LAPACKE_dsytri2_work LAPACKE_NAME(dsytri2_work,DSYTRI2_WORK)
#define LAPACKE_dsytri2x_work LAPACKE_NAME(dsytri2x_work,DSYTRI2X_WORK)
#define LAPACKE_dsytrs2_work LAPACKE_NAME(dsytrs2_work,DSYTRS2_WORK)
#define LAPACKE_sbbcsd_work LAPACKE_NAME(sbbcsd_work,SBBCSD_WORK)
#define LAPACKE_sorbdb_work LAPACKE_NAME(sorbdb_work,SORBDB_WORK)
#define LAPACKE_sorcsd_work LAPACKE_NAME(sorcsd_work,SORCSD_WORK)
#define LAPACKE_ssyconv_work LAPACKE_NAME(ssyconv_work,SSYCONV_WORK)
#define LAPACKE_ssyswapr_work LAPACKE_NAME(ssyswapr_work,SSYSWAPR_WORK)
#define LAPACKE_ssytri2_work LAPACKE_NAME(ssytri2_work,SSYTRI2_WORK)
#define LAPACKE_ssytri2x_work LAPACKE_NAME(ssytri2x_work,SSYTRI2X_WORK)
#define LAPACKE_ssytrs2_work LAPACKE_NAME(ssytrs2_work,SSYTRS2_WORK)
#define LAPACKE_zbbcsd_work LAPACKE_NAME(zbbcsd_work,ZBBCSD_WORK)
#define LAPACKE_zheswapr_work LAPACKE_NAME(zheswapr_work,ZHESWAPR_WORK)
#define LAPACKE_zhetri2_work LAPACKE_NAME(zhetri2_work,ZHETRI2_WORK)
#define LAPACKE_zhetri2x_work LAPACKE_NAME(zhetri2x_work,ZHETRI2X_WORK)
#define LAPACKE_zhetrs2_work LAPACKE_NAME(zhetrs2_work,ZHETRS2_WORK)
#define LAPACKE_zsyconv_work LAPACKE_NAME(zsyconv_work,ZSYCONV_WORK)
#define LAPACKE_zsyswapr_work LAPACKE_NAME(zsyswapr_work,ZSYSWAPR_WORK)
#define LAPACKE_zsytri2_work LAPACKE_NAME(zsytri2_work,ZSYTRI2_WORK)
#define LAPACKE_zsytri2x_work LAPACKE_NAME(zsytri2x_work,ZSYTRI2X_WORK)
#define LAPACKE_zsytrs2_work LAPACKE_NAME(zsytrs2_work,ZSYTRS2_WORK)
#define LAPACKE_zunbdb_work LAPACKE_NAME(zunbdb_work,ZUNBDB_WORK)
#define LAPACKE_zuncsd_work LAPACKE_NAME(zuncsd_work,ZUNCSD_WORK)

//LAPACK 3.4.0
#define LAPACKE_sgemqrt   LAPACKE_NAME(sgemqrt,SGEMQRT)
#define LAPACKE_dgemqrt   LAPACKE_NAME(dgemqrt,DGEMQRT)
#define LAPACKE_cgemqrt   LAPACKE_NAME(cgemqrt,CGEMQRT)
#define LAPACKE_zgemqrt   LAPACKE_NAME(zgemqrt,ZGEMQRT)

#define LAPACKE_sgeqrt   LAPACKE_NAME(sgeqrt,SGEQRT)
#define LAPACKE_dgeqrt   LAPACKE_NAME(dgeqrt,DGEQRT)
#define LAPACKE_cgeqrt   LAPACKE_NAME(cgeqrt,CGEQRT)
#define LAPACKE_zgeqrt   LAPACKE_NAME(zgeqrt,ZGEQRT)

#define LAPACKE_sgeqrt2   LAPACKE_NAME(sgeqrt2,SGEQRT2)
#define LAPACKE_dgeqrt2   LAPACKE_NAME(dgeqrt2,DGEQRT2)
#define LAPACKE_cgeqrt2   LAPACKE_NAME(cgeqrt2,CGEQRT2)
#define LAPACKE_zgeqrt2   LAPACKE_NAME(zgeqrt2,ZGEQRT2)

#define LAPACKE_sgeqrt3   LAPACKE_NAME(sgeqrt3,SGEQRT3)
#define LAPACKE_dgeqrt3   LAPACKE_NAME(dgeqrt3,DGEQRT3)
#define LAPACKE_cgeqrt3   LAPACKE_NAME(cgeqrt3,CGEQRT3)
#define LAPACKE_zgeqrt3   LAPACKE_NAME(zgeqrt3,ZGEQRT3)

#define LAPACKE_stpmqrt   LAPACKE_NAME(stpmqrt,STPMQRT)
#define LAPACKE_dtpmqrt   LAPACKE_NAME(dtpmqrt,DTPMQRT)
#define LAPACKE_ctpmqrt   LAPACKE_NAME(ctpmqrt,CTPMQRT)
#define LAPACKE_ztpmqrt   LAPACKE_NAME(ztpmqrt,ZTPMQRT)

#define LAPACKE_dtpqrt   LAPACKE_NAME(dtpqrt,DTPQRT)
#define LAPACKE_ctpqrt   LAPACKE_NAME(ctpqrt,CTPQRT)
#define LAPACKE_ztpqrt   LAPACKE_NAME(ztpqrt,ZTPQRT)

#define LAPACKE_stpqrt2   LAPACKE_NAME(stpqrt2,STPQRT2)
#define LAPACKE_dtpqrt2   LAPACKE_NAME(dtpqrt2,DTPQRT2)
#define LAPACKE_ctpqrt2   LAPACKE_NAME(ctpqrt2,CTPQRT2)
#define LAPACKE_ztpqrt2   LAPACKE_NAME(ztpqrt2,ZTPQRT2)

#define LAPACKE_stprfb   LAPACKE_NAME(stprfb,STPRFB)
#define LAPACKE_dtprfb   LAPACKE_NAME(dtprfb,DTPRFB)
#define LAPACKE_ctprfb   LAPACKE_NAME(ctprfb,CTPRFB)
#define LAPACKE_ztprfb   LAPACKE_NAME(ztprfb,ZTPRFB)

#define LAPACKE_sgemqrt_work   LAPACKE_NAME(sgemqrt_work,SGEMQRT_WORK)
#define LAPACKE_dgemqrt_work   LAPACKE_NAME(dgemqrt_work,DGEMQRT_WORK)
#define LAPACKE_cgemqrt_work   LAPACKE_NAME(cgemqrt_work,CGEMQRT_WORK)
#define LAPACKE_zgemqrt_work   LAPACKE_NAME(zgemqrt_work,ZGEMQRT_WORK)

#define LAPACKE_sgeqrt_work   LAPACKE_NAME(sgeqrt_work,SGEQRT_WORK)
#define LAPACKE_dgeqrt_work   LAPACKE_NAME(dgeqrt_work,DGEQRT_WORK)
#define LAPACKE_cgeqrt_work   LAPACKE_NAME(cgeqrt_work,CGEQRT_WORK)
#define LAPACKE_zgeqrt_work   LAPACKE_NAME(zgeqrt_work,ZGEQRT_WORK)

#define LAPACKE_sgeqrt2_work   LAPACKE_NAME(sgeqrt2_work,SGEQRT2_WORK)
#define LAPACKE_dgeqrt2_work   LAPACKE_NAME(dgeqrt2_work,DGEQRT2_WORK)
#define LAPACKE_cgeqrt2_work   LAPACKE_NAME(cgeqrt2_work,CGEQRT2_WORK)
#define LAPACKE_zgeqrt2_work   LAPACKE_NAME(zgeqrt2_work,ZGEQRT2_WORK)

#define LAPACKE_sgeqrt3_work   LAPACKE_NAME(sgeqrt3_work,SGEQRT3_WORK)
#define LAPACKE_dgeqrt3_work   LAPACKE_NAME(dgeqrt3_work,DGEQRT3_WORK)
#define LAPACKE_cgeqrt3_work   LAPACKE_NAME(cgeqrt3_work,CGEQRT3_WORK)
#define LAPACKE_zgeqrt3_work   LAPACKE_NAME(zgeqrt3_work,ZGEQRT3_WORK)

#define LAPACKE_stpmqrt_work   LAPACKE_NAME(stpmqrt_work,STPMQRT_WORK)
#define LAPACKE_dtpmqrt_work   LAPACKE_NAME(dtpmqrt_work,DTPMQRT_WORK)
#define LAPACKE_ctpmqrt_work   LAPACKE_NAME(ctpmqrt_work,CTPMQRT_WORK)
#define LAPACKE_ztpmqrt_work   LAPACKE_NAME(ztpmqrt_work,ZTPMQRT_WORK)

#define LAPACKE_dtpqrt_work   LAPACKE_NAME(dtpqrt_work,DTPQRT_WORK)
#define LAPACKE_ctpqrt_work   LAPACKE_NAME(ctpqrt_work,CTPQRT_WORK)
#define LAPACKE_ztpqrt_work   LAPACKE_NAME(ztpqrt_work,ZTPQRT_WORK)

#define LAPACKE_stpqrt2_work   LAPACKE_NAME(stpqrt2_work,STPQRT2_WORK)
#define LAPACKE_dtpqrt2_work   LAPACKE_NAME(dtpqrt2_work,DTPQRT2_WORK)
#define LAPACKE_ctpqrt2_work   LAPACKE_NAME(ctpqrt2_work,CTPQRT2_WORK)
#define LAPACKE_ztpqrt2_work   LAPACKE_NAME(ztpqrt2_work,ZTPQRT2_WORK)

#define LAPACKE_stprfb_work   LAPACKE_NAME(stprfb_work,STPRFB_WORK)
#define LAPACKE_dtprfb_work   LAPACKE_NAME(dtprfb_work,DTPRFB_WORK)
#define LAPACKE_ctprfb_work   LAPACKE_NAME(ctprfb_work,CTPRFB_WORK)
#define LAPACKE_ztprfb_work   LAPACKE_NAME(ztprfb_work,ZTPRFB_WORK)

#define LAPACKE_cgb_trans LAPACKE_NAME(cgb_trans,CGB_TRANS)
#define LAPACKE_cge_trans LAPACKE_NAME(cge_trans,CGE_TRANS)
#define LAPACKE_cgg_trans LAPACKE_NAME(cgg_trans,CGG_TRANS)
#define LAPACKE_chb_trans LAPACKE_NAME(chb_trans,CHB_TRANS)
#define LAPACKE_che_trans LAPACKE_NAME(che_trans,CHE_TRANS)
#define LAPACKE_chp_trans LAPACKE_NAME(chp_trans,CHP_TRANS)
#define LAPACKE_chs_trans LAPACKE_NAME(chs_trans,CHS_TRANS)
#define LAPACKE_cpb_trans LAPACKE_NAME(cpb_trans,CPB_TRANS)
#define LAPACKE_cpf_trans LAPACKE_NAME(cpf_trans,CPF_TRANS)
#define LAPACKE_cpo_trans LAPACKE_NAME(cpo_trans,CPO_TRANS)
#define LAPACKE_cpp_trans LAPACKE_NAME(cpp_trans,CPP_TRANS)
#define LAPACKE_csp_trans LAPACKE_NAME(csp_trans,CSP_TRANS)
#define LAPACKE_csy_trans LAPACKE_NAME(csy_trans,CSY_TRANS)
#define LAPACKE_ctb_trans LAPACKE_NAME(ctb_trans,CTB_TRANS)
#define LAPACKE_ctf_trans LAPACKE_NAME(ctf_trans,CTF_TRANS)
#define LAPACKE_ctp_trans LAPACKE_NAME(ctp_trans,CTP_TRANS)
#define LAPACKE_ctr_trans LAPACKE_NAME(ctr_trans,CTR_TRANS)
#define LAPACKE_dgb_trans LAPACKE_NAME(dgb_trans,DGB_TRANS)
#define LAPACKE_dge_trans LAPACKE_NAME(dge_trans,DGE_TRANS)
#define LAPACKE_dgg_trans LAPACKE_NAME(dgg_trans,DGG_TRANS)
#define LAPACKE_dhs_trans LAPACKE_NAME(dhs_trans,DHS_TRANS)
#define LAPACKE_dpb_trans LAPACKE_NAME(dpb_trans,DPB_TRANS)
#define LAPACKE_dpf_trans LAPACKE_NAME(dpf_trans,DPF_TRANS)
#define LAPACKE_dpo_trans LAPACKE_NAME(dpo_trans,DPO_TRANS)
#define LAPACKE_dpp_trans LAPACKE_NAME(dpp_trans,DPP_TRANS)
#define LAPACKE_dsb_trans LAPACKE_NAME(dsb_trans,DSB_TRANS)
#define LAPACKE_dsp_trans LAPACKE_NAME(dsp_trans,DSP_TRANS)
#define LAPACKE_dsy_trans LAPACKE_NAME(dsy_trans,DSY_TRANS)
#define LAPACKE_dtb_trans LAPACKE_NAME(dtb_trans,DTB_TRANS)
#define LAPACKE_dtf_trans LAPACKE_NAME(dtf_trans,DTF_TRANS)
#define LAPACKE_dtp_trans LAPACKE_NAME(dtp_trans,DTP_TRANS)
#define LAPACKE_dtr_trans LAPACKE_NAME(dtr_trans,DTR_TRANS)
#define LAPACKE_sgb_trans LAPACKE_NAME(sgb_trans,SGB_TRANS)
#define LAPACKE_sge_trans LAPACKE_NAME(sge_trans,SGE_TRANS)
#define LAPACKE_sgg_trans LAPACKE_NAME(sgg_trans,SGG_TRANS)
#define LAPACKE_shs_trans LAPACKE_NAME(shs_trans,SHS_TRANS)
#define LAPACKE_spb_trans LAPACKE_NAME(spb_trans,SPB_TRANS)
#define LAPACKE_spf_trans LAPACKE_NAME(spf_trans,SPF_TRANS)
#define LAPACKE_spo_trans LAPACKE_NAME(spo_trans,SPO_TRANS)
#define LAPACKE_spp_trans LAPACKE_NAME(spp_trans,SPP_TRANS)
#define LAPACKE_ssb_trans LAPACKE_NAME(ssb_trans,SSB_TRANS)
#define LAPACKE_ssp_trans LAPACKE_NAME(ssp_trans,SSP_TRANS)
#define LAPACKE_ssy_trans LAPACKE_NAME(ssy_trans,SSY_TRANS)
#define LAPACKE_stb_trans LAPACKE_NAME(stb_trans,STB_TRANS)
#define LAPACKE_stf_trans LAPACKE_NAME(stf_trans,STF_TRANS)
#define LAPACKE_stp_trans LAPACKE_NAME(stp_trans,STP_TRANS)
#define LAPACKE_str_trans LAPACKE_NAME(str_trans,STR_TRANS)
#define LAPACKE_zgb_trans LAPACKE_NAME(zgb_trans,ZGB_TRANS)
#define LAPACKE_zge_trans LAPACKE_NAME(zge_trans,ZGE_TRANS)
#define LAPACKE_zgg_trans LAPACKE_NAME(zgg_trans,ZGG_TRANS)
#define LAPACKE_zhb_trans LAPACKE_NAME(zhb_trans,ZHB_TRANS)
#define LAPACKE_zhe_trans LAPACKE_NAME(zhe_trans,ZHE_TRANS)
#define LAPACKE_zhp_trans LAPACKE_NAME(zhp_trans,ZHP_TRANS)
#define LAPACKE_zhs_trans LAPACKE_NAME(zhs_trans,ZHS_TRANS)
#define LAPACKE_zpb_trans LAPACKE_NAME(zpb_trans,ZPB_TRANS)
#define LAPACKE_zpf_trans LAPACKE_NAME(zpf_trans,ZPF_TRANS)
#define LAPACKE_zpo_trans LAPACKE_NAME(zpo_trans,ZPO_TRANS)
#define LAPACKE_zpp_trans LAPACKE_NAME(zpp_trans,ZPP_TRANS)
#define LAPACKE_zsp_trans LAPACKE_NAME(zsp_trans,ZSP_TRANS)
#define LAPACKE_zsy_trans LAPACKE_NAME(zsy_trans,ZSY_TRANS)
#define LAPACKE_ztb_trans LAPACKE_NAME(ztb_trans,ZTB_TRANS)
#define LAPACKE_ztf_trans LAPACKE_NAME(ztf_trans,ZTF_TRANS)
#define LAPACKE_ztp_trans LAPACKE_NAME(ztp_trans,ZTP_TRANS)
#define LAPACKE_ztr_trans LAPACKE_NAME(ztr_trans,ZTR_TRANS)

#define LAPACKE_c_nancheck LAPACKE_NAME(c_nancheck,C_NANCHECK)
#define LAPACKE_d_nancheck LAPACKE_NAME(d_nancheck,D_NANCHECK)
#define LAPACKE_s_nancheck LAPACKE_NAME(s_nancheck,S_NANCHECK)
#define LAPACKE_z_nancheck LAPACKE_NAME(z_nancheck,Z_NANCHECK)
#define LAPACKE_cgb_nancheck LAPACKE_NAME(cgb_nancheck,CGB_NANCHECK)
#define LAPACKE_cge_nancheck LAPACKE_NAME(cge_nancheck,CGE_NANCHECK)
#define LAPACKE_cgg_nancheck LAPACKE_NAME(cgg_nancheck,CGG_NANCHECK)
#define LAPACKE_cgt_nancheck LAPACKE_NAME(cgt_nancheck,CGT_NANCHECK)
#define LAPACKE_chb_nancheck LAPACKE_NAME(chb_nancheck,CHB_NANCHECK)
#define LAPACKE_che_nancheck LAPACKE_NAME(che_nancheck,CHE_NANCHECK)
#define LAPACKE_chp_nancheck LAPACKE_NAME(chp_nancheck,CHP_NANCHECK)
#define LAPACKE_chs_nancheck LAPACKE_NAME(chs_nancheck,CHS_NANCHECK)
#define LAPACKE_cpb_nancheck LAPACKE_NAME(cpb_nancheck,CPB_NANCHECK)
#define LAPACKE_cpf_nancheck LAPACKE_NAME(cpf_nancheck,CPF_NANCHECK)
#define LAPACKE_cpo_nancheck LAPACKE_NAME(cpo_nancheck,CPO_NANCHECK)
#define LAPACKE_cpp_nancheck LAPACKE_NAME(cpp_nancheck,CPP_NANCHECK)
#define LAPACKE_cpt_nancheck LAPACKE_NAME(cpt_nancheck,CPT_NANCHECK)
#define LAPACKE_csp_nancheck LAPACKE_NAME(csp_nancheck,CSP_NANCHECK)
#define LAPACKE_cst_nancheck LAPACKE_NAME(cst_nancheck,CST_NANCHECK)
#define LAPACKE_csy_nancheck LAPACKE_NAME(csy_nancheck,CSY_NANCHECK)
#define LAPACKE_ctb_nancheck LAPACKE_NAME(ctb_nancheck,CTB_NANCHECK)
#define LAPACKE_ctf_nancheck LAPACKE_NAME(ctf_nancheck,CTF_NANCHECK)
#define LAPACKE_ctp_nancheck LAPACKE_NAME(ctp_nancheck,CTP_NANCHECK)
#define LAPACKE_ctr_nancheck LAPACKE_NAME(ctr_nancheck,CTR_NANCHECK)
#define LAPACKE_dgb_nancheck LAPACKE_NAME(dgb_nancheck,DGB_NANCHECK)
#define LAPACKE_dge_nancheck LAPACKE_NAME(dge_nancheck,DGE_NANCHECK)
#define LAPACKE_dgg_nancheck LAPACKE_NAME(dgg_nancheck,DGG_NANCHECK)
#define LAPACKE_dgt_nancheck LAPACKE_NAME(dgt_nancheck,DGT_NANCHECK)
#define LAPACKE_dhs_nancheck LAPACKE_NAME(dhs_nancheck,DHS_NANCHECK)
#define LAPACKE_dpb_nancheck LAPACKE_NAME(dpb_nancheck,DPB_NANCHECK)
#define LAPACKE_dpf_nancheck LAPACKE_NAME(dpf_nancheck,DPF_NANCHECK)
#define LAPACKE_dpo_nancheck LAPACKE_NAME(dpo_nancheck,DPO_NANCHECK)
#define LAPACKE_dpp_nancheck LAPACKE_NAME(dpp_nancheck,DPP_NANCHECK)
#define LAPACKE_dpt_nancheck LAPACKE_NAME(dpt_nancheck,DPT_NANCHECK)
#define LAPACKE_dsb_nancheck LAPACKE_NAME(dsb_nancheck,DSB_NANCHECK)
#define LAPACKE_dsp_nancheck LAPACKE_NAME(dsp_nancheck,DSP_NANCHECK)
#define LAPACKE_dst_nancheck LAPACKE_NAME(dst_nancheck,DST_NANCHECK)
#define LAPACKE_dsy_nancheck LAPACKE_NAME(dsy_nancheck,DSY_NANCHECK)
#define LAPACKE_dtb_nancheck LAPACKE_NAME(dtb_nancheck,DTB_NANCHECK)
#define LAPACKE_dtf_nancheck LAPACKE_NAME(dtf_nancheck,DTF_NANCHECK)
#define LAPACKE_dtp_nancheck LAPACKE_NAME(dtp_nancheck,DTP_NANCHECK)
#define LAPACKE_dtr_nancheck LAPACKE_NAME(dtr_nancheck,DTR_NANCHECK)
#define LAPACKE_sgb_nancheck LAPACKE_NAME(sgb_nancheck,SGB_NANCHECK)
#define LAPACKE_sge_nancheck LAPACKE_NAME(sge_nancheck,SGE_NANCHECK)
#define LAPACKE_sgg_nancheck LAPACKE_NAME(sgg_nancheck,SGG_NANCHECK)
#define LAPACKE_sgt_nancheck LAPACKE_NAME(sgt_nancheck,SGT_NANCHECK)
#define LAPACKE_shs_nancheck LAPACKE_NAME(shs_nancheck,SHS_NANCHECK)
#define LAPACKE_spb_nancheck LAPACKE_NAME(spb_nancheck,SPB_NANCHECK)
#define LAPACKE_spf_nancheck LAPACKE_NAME(spf_nancheck,SPF_NANCHECK)
#define LAPACKE_spo_nancheck LAPACKE_NAME(spo_nancheck,SPO_NANCHECK)
#define LAPACKE_spp_nancheck LAPACKE_NAME(spp_nancheck,SPP_NANCHECK)
#define LAPACKE_spt_nancheck LAPACKE_NAME(spt_nancheck,SPT_NANCHECK)
#define LAPACKE_ssb_nancheck LAPACKE_NAME(ssb_nancheck,SSB_NANCHECK)
#define LAPACKE_ssp_nancheck LAPACKE_NAME(ssp_nancheck,SSP_NANCHECK)
#define LAPACKE_sst_nancheck LAPACKE_NAME(sst_nancheck,SST_NANCHECK)
#define LAPACKE_ssy_nancheck LAPACKE_NAME(ssy_nancheck,SSY_NANCHECK)
#define LAPACKE_stb_nancheck LAPACKE_NAME(stb_nancheck,STB_NANCHECK)
#define LAPACKE_stf_nancheck LAPACKE_NAME(stf_nancheck,STF_NANCHECK)
#define LAPACKE_stp_nancheck LAPACKE_NAME(stp_nancheck,STP_NANCHECK)
#define LAPACKE_str_nancheck LAPACKE_NAME(str_nancheck,STR_NANCHECK)
#define LAPACKE_zgb_nancheck LAPACKE_NAME(zgb_nancheck,ZGB_NANCHECK)
#define LAPACKE_zge_nancheck LAPACKE_NAME(zge_nancheck,ZGE_NANCHECK)
#define LAPACKE_zgg_nancheck LAPACKE_NAME(zgg_nancheck,ZGG_NANCHECK)
#define LAPACKE_zgt_nancheck LAPACKE_NAME(zgt_nancheck,ZGT_NANCHECK)
#define LAPACKE_zhb_nancheck LAPACKE_NAME(zhb_nancheck,ZHB_NANCHECK)
#define LAPACKE_zhe_nancheck LAPACKE_NAME(zhe_nancheck,ZHE_NANCHECK)
#define LAPACKE_zhp_nancheck LAPACKE_NAME(zhp_nancheck,ZHP_NANCHECK)
#define LAPACKE_zhs_nancheck LAPACKE_NAME(zhs_nancheck,ZHS_NANCHECK)
#define LAPACKE_zpb_nancheck LAPACKE_NAME(zpb_nancheck,ZPB_NANCHECK)
#define LAPACKE_zpf_nancheck LAPACKE_NAME(zpf_nancheck,ZPF_NANCHECK)
#define LAPACKE_zpo_nancheck LAPACKE_NAME(zpo_nancheck,ZPO_NANCHECK)
#define LAPACKE_zpp_nancheck LAPACKE_NAME(zpp_nancheck,ZPP_NANCHECK)
#define LAPACKE_zpt_nancheck LAPACKE_NAME(zpt_nancheck,ZPT_NANCHECK)
#define LAPACKE_zsp_nancheck LAPACKE_NAME(zsp_nancheck,ZSP_NANCHECK)
#define LAPACKE_zst_nancheck LAPACKE_NAME(zst_nancheck,ZST_NANCHECK)
#define LAPACKE_zsy_nancheck LAPACKE_NAME(zsy_nancheck,ZSY_NANCHECK)
#define LAPACKE_ztb_nancheck LAPACKE_NAME(ztb_nancheck,ZTB_NANCHECK)
#define LAPACKE_ztf_nancheck LAPACKE_NAME(ztf_nancheck,ZTF_NANCHECK)
#define LAPACKE_ztp_nancheck LAPACKE_NAME(ztp_nancheck,ZTP_NANCHECK)
#define LAPACKE_ztr_nancheck LAPACKE_NAME(ztr_nancheck,ZTR_NANCHECK)

#endif /* LAPACK_NAME_PATTERN_MC */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LAPACKE_CONFIG_H_ */
