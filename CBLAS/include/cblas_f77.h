/*
 * cblas_f77.h
 * Written by Keita Teranishi
 *
 * Updated by Jeff Horner
 * Merged cblas_f77.h and cblas_fortran_header.h
 */

#ifndef CBLAS_F77_H
#define CBLAS_F77_H

#include <stdarg.h>

/* It seems all current Fortran compilers put strlen at end.
*  Some historical compilers put strlen after the str argument
*  or make the str argument into a struct. */
#define BLAS_FORTRAN_STRLEN_END

#ifdef CRAY
   #include <fortran.h>
   #define F77_CHAR _fcd
   #define C2F_CHAR(a) ( _cptofcd( (a), 1 ) )
   #define C2F_STR(a, i) ( _cptofcd( (a), (i) ) )
   #define F77_STRLEN(a) (_fcdlen)
#endif

#ifdef WeirdNEC
   #define F77_INT long
#else
   #define F77_INT int
#endif

#ifdef  F77_CHAR
   #define FCHAR F77_CHAR
#else
   #define FCHAR char *
#endif

#define FINT const F77_INT *
#define FINT2 F77_INT *

/*
 * Level 1 BLAS
 */

#define F77_xerbla_base 		F77_GLOBAL(xerbla,XERBLA)
#define F77_srotg_base 		F77_GLOBAL(srotg,SROTG)
#define F77_srotmg_base 		F77_GLOBAL(srotmg,SROTMG)
#define F77_srot_base 		F77_GLOBAL(srot,SROT)
#define F77_srotm_base 		F77_GLOBAL(srotm,SROTM)
#define F77_drotg_base 		F77_GLOBAL(drotg,DROTG)
#define F77_drotmg_base 		F77_GLOBAL(drotmg,DROTMG)
#define F77_drot_base 		F77_GLOBAL(drot,DROT)
#define F77_drotm_base 		F77_GLOBAL(drotm,DROTM)
#define F77_sswap_base 		F77_GLOBAL(sswap,SSWAP)
#define F77_scopy_base 		F77_GLOBAL(scopy,SCOPY)
#define F77_saxpy_base 		F77_GLOBAL(saxpy,SAXPY)
#define F77_isamax_sub_base 	F77_GLOBAL(isamaxsub,ISAMAXSUB)
#define F77_dswap_base 		F77_GLOBAL(dswap,DSWAP)
#define F77_dcopy_base 		F77_GLOBAL(dcopy,DCOPY)
#define F77_daxpy_base 		F77_GLOBAL(daxpy,DAXPY)
#define F77_idamax_sub_base 	F77_GLOBAL(idamaxsub,IDAMAXSUB)
#define F77_cswap_base 		F77_GLOBAL(cswap,CSWAP)
#define F77_ccopy_base 		F77_GLOBAL(ccopy,CCOPY)
#define F77_caxpy_base 		F77_GLOBAL(caxpy,CAXPY)
#define F77_icamax_sub_base 	F77_GLOBAL(icamaxsub,ICAMAXSUB)
#define F77_zswap_base 		F77_GLOBAL(zswap,ZSWAP)
#define F77_zcopy_base 		F77_GLOBAL(zcopy,ZCOPY)
#define F77_zaxpy_base 		F77_GLOBAL(zaxpy,ZAXPY)
#define F77_izamax_sub_base 	F77_GLOBAL(izamaxsub,IZAMAXSUB)
#define F77_sdot_sub_base 	F77_GLOBAL(sdotsub,SDOTSUB)
#define F77_ddot_sub_base 	F77_GLOBAL(ddotsub,DDOTSUB)
#define F77_dsdot_sub_base 	F77_GLOBAL(dsdotsub,DSDOTSUB)
#define F77_sscal_base 		F77_GLOBAL(sscal,SSCAL)
#define F77_dscal_base 		F77_GLOBAL(dscal,DSCAL)
#define F77_cscal_base 		F77_GLOBAL(cscal,CSCAL)
#define F77_zscal_base 		F77_GLOBAL(zscal,ZSCAL)
#define F77_csscal_base 		F77_GLOBAL(csscal,CSSCAL)
#define F77_zdscal_base 		F77_GLOBAL(zdscal,ZDSCAL)
#define F77_cdotu_sub_base 	F77_GLOBAL(cdotusub,CDOTUSUB)
#define F77_cdotc_sub_base 	F77_GLOBAL(cdotcsub,CDOTCSUB)
#define F77_zdotu_sub_base 	F77_GLOBAL(zdotusub,ZDOTUSUB)
#define F77_zdotc_sub_base 	F77_GLOBAL(zdotcsub,ZDOTCSUB)
#define F77_snrm2_sub_base 	F77_GLOBAL(snrm2sub,SNRM2SUB)
#define F77_sasum_sub_base 	F77_GLOBAL(sasumsub,SASUMSUB)
#define F77_dnrm2_sub_base 	F77_GLOBAL(dnrm2sub,DNRM2SUB)
#define F77_dasum_sub_base 	F77_GLOBAL(dasumsub,DASUMSUB)
#define F77_scnrm2_sub_base 	F77_GLOBAL(scnrm2sub,SCNRM2SUB)
#define F77_scasum_sub_base 	F77_GLOBAL(scasumsub,SCASUMSUB)
#define F77_dznrm2_sub_base 	F77_GLOBAL(dznrm2sub,DZNRM2SUB)
#define F77_dzasum_sub_base 	F77_GLOBAL(dzasumsub,DZASUMSUB)
#define F77_sdsdot_sub_base 	F77_GLOBAL(sdsdotsub,SDSDOTSUB)
/*
 * Level 2 BLAS
 */
#define F77_ssymv_base 		F77_GLOBAL(ssymv,SSYMV)
#define F77_ssbmv_base 		F77_GLOBAL(ssbmv,SSBMV)
#define F77_sspmv_base 		F77_GLOBAL(sspmv,SSPMV)
#define F77_sger_base 		F77_GLOBAL(sger,SGER)
#define F77_ssyr_base 		F77_GLOBAL(ssyr,SSYR)
#define F77_sspr_base 		F77_GLOBAL(sspr,SSPR)
#define F77_ssyr2_base 		F77_GLOBAL(ssyr2,SSYR2)
#define F77_sspr2_base 		F77_GLOBAL(sspr2,SSPR2)
#define F77_dsymv_base 		F77_GLOBAL(dsymv,DSYMV)
#define F77_dsbmv_base 		F77_GLOBAL(dsbmv,DSBMV)
#define F77_dspmv_base 		F77_GLOBAL(dspmv,DSPMV)
#define F77_dger_base 		F77_GLOBAL(dger,DGER)
#define F77_dsyr_base 		F77_GLOBAL(dsyr,DSYR)
#define F77_dspr_base 		F77_GLOBAL(dspr,DSPR)
#define F77_dsyr2_base 		F77_GLOBAL(dsyr2,DSYR2)
#define F77_dspr2_base 		F77_GLOBAL(dspr2,DSPR2)
#define F77_chemv_base 		F77_GLOBAL(chemv,CHEMV)
#define F77_chbmv_base 		F77_GLOBAL(chbmv,CHBMV)
#define F77_chpmv_base 		F77_GLOBAL(chpmv,CHPMV)
#define F77_cgeru_base 		F77_GLOBAL(cgeru,CGERU)
#define F77_cgerc_base 		F77_GLOBAL(cgerc,CGERC)
#define F77_cher_base 		F77_GLOBAL(cher,CHER)
#define F77_chpr_base 		F77_GLOBAL(chpr,CHPR)
#define F77_cher2_base 		F77_GLOBAL(cher2,CHER2)
#define F77_chpr2_base 		F77_GLOBAL(chpr2,CHPR2)
#define F77_zhemv_base 		F77_GLOBAL(zhemv,ZHEMV)
#define F77_zhbmv_base 		F77_GLOBAL(zhbmv,ZHBMV)
#define F77_zhpmv_base 		F77_GLOBAL(zhpmv,ZHPMV)
#define F77_zgeru_base 		F77_GLOBAL(zgeru,ZGERU)
#define F77_zgerc_base 		F77_GLOBAL(zgerc,ZGERC)
#define F77_zher_base 		F77_GLOBAL(zher,ZHER)
#define F77_zhpr_base 		F77_GLOBAL(zhpr,ZHPR)
#define F77_zher2_base 		F77_GLOBAL(zher2,ZHER2)
#define F77_zhpr2_base 		F77_GLOBAL(zhpr2,ZHPR2)
#define F77_sgemv_base 		F77_GLOBAL(sgemv,SGEMV)
#define F77_sgbmv_base 		F77_GLOBAL(sgbmv,SGBMV)
#define F77_strmv_base 		F77_GLOBAL(strmv,STRMV)
#define F77_stbmv_base 		F77_GLOBAL(stbmv,STBMV)
#define F77_stpmv_base 		F77_GLOBAL(stpmv,STPMV)
#define F77_strsv_base 		F77_GLOBAL(strsv,STRSV)
#define F77_stbsv_base 		F77_GLOBAL(stbsv,STBSV)
#define F77_stpsv_base 		F77_GLOBAL(stpsv,STPSV)
#define F77_dgemv_base 		F77_GLOBAL(dgemv,DGEMV)
#define F77_dgbmv_base 		F77_GLOBAL(dgbmv,DGBMV)
#define F77_dtrmv_base 		F77_GLOBAL(dtrmv,DTRMV)
#define F77_dtbmv_base 		F77_GLOBAL(dtbmv,DTBMV)
#define F77_dtpmv_base 		F77_GLOBAL(dtpmv,DTPMV)
#define F77_dtrsv_base 		F77_GLOBAL(dtrsv,DTRSV)
#define F77_dtbsv_base 		F77_GLOBAL(dtbsv,DTBSV)
#define F77_dtpsv_base 		F77_GLOBAL(dtpsv,DTPSV)
#define F77_cgemv_base 		F77_GLOBAL(cgemv,CGEMV)
#define F77_cgbmv_base 		F77_GLOBAL(cgbmv,CGBMV)
#define F77_ctrmv_base 		F77_GLOBAL(ctrmv,CTRMV)
#define F77_ctbmv_base 		F77_GLOBAL(ctbmv,CTBMV)
#define F77_ctpmv_base 		F77_GLOBAL(ctpmv,CTPMV)
#define F77_ctrsv_base 		F77_GLOBAL(ctrsv,CTRSV)
#define F77_ctbsv_base 		F77_GLOBAL(ctbsv,CTBSV)
#define F77_ctpsv_base 		F77_GLOBAL(ctpsv,CTPSV)
#define F77_zgemv_base 		F77_GLOBAL(zgemv,ZGEMV)
#define F77_zgbmv_base 		F77_GLOBAL(zgbmv,ZGBMV)
#define F77_ztrmv_base 		F77_GLOBAL(ztrmv,ZTRMV)
#define F77_ztbmv_base 		F77_GLOBAL(ztbmv,ZTBMV)
#define F77_ztpmv_base 		F77_GLOBAL(ztpmv,ZTPMV)
#define F77_ztrsv_base 		F77_GLOBAL(ztrsv,ZTRSV)
#define F77_ztbsv_base 		F77_GLOBAL(ztbsv,ZTBSV)
#define F77_ztpsv_base 		F77_GLOBAL(ztpsv,ZTPSV)
/*
 * Level 3 BLAS
 */
#define F77_chemm_base 		F77_GLOBAL(chemm,CHEMM)
#define F77_cherk_base 		F77_GLOBAL(cherk,CHERK)
#define F77_cher2k_base 		F77_GLOBAL(cher2k,CHER2K)
#define F77_zhemm_base 		F77_GLOBAL(zhemm,ZHEMM)
#define F77_zherk_base 		F77_GLOBAL(zherk,ZHERK)
#define F77_zher2k_base 		F77_GLOBAL(zher2k,ZHER2K)
#define F77_sgemm_base 		F77_GLOBAL(sgemm,SGEMM)
#define F77_ssymm_base 		F77_GLOBAL(ssymm,SSYMM)
#define F77_ssyrk_base 		F77_GLOBAL(ssyrk,SSYRK)
#define F77_ssyr2k_base 		F77_GLOBAL(ssyr2k,SSYR2K)
#define F77_strmm_base 		F77_GLOBAL(strmm,STRMM)
#define F77_strsm_base 		F77_GLOBAL(strsm,STRSM)
#define F77_dgemm_base 		F77_GLOBAL(dgemm,DGEMM)
#define F77_dsymm_base 		F77_GLOBAL(dsymm,DSYMM)
#define F77_dsyrk_base 		F77_GLOBAL(dsyrk,DSYRK)
#define F77_dsyr2k_base 		F77_GLOBAL(dsyr2k,DSYR2K)
#define F77_dtrmm_base 		F77_GLOBAL(dtrmm,DTRMM)
#define F77_dtrsm_base 		F77_GLOBAL(dtrsm,DTRSM)
#define F77_cgemm_base 		F77_GLOBAL(cgemm,CGEMM)
#define F77_csymm_base 		F77_GLOBAL(csymm,CSYMM)
#define F77_csyrk_base 		F77_GLOBAL(csyrk,CSYRK)
#define F77_csyr2k_base 		F77_GLOBAL(csyr2k,CSYR2K)
#define F77_ctrmm_base 		F77_GLOBAL(ctrmm,CTRMM)
#define F77_ctrsm_base 		F77_GLOBAL(ctrsm,CTRSM)
#define F77_zgemm_base 		F77_GLOBAL(zgemm,ZGEMM)
#define F77_zsymm_base 		F77_GLOBAL(zsymm,ZSYMM)
#define F77_zsyrk_base 		F77_GLOBAL(zsyrk,ZSYRK)
#define F77_zsyr2k_base 		F77_GLOBAL(zsyr2k,ZSYR2K)
#define F77_ztrmm_base 		F77_GLOBAL(ztrmm,ZTRMM)
#define F77_ztrsm_base 		F77_GLOBAL(ztrsm,ZTRSM)

/*
 * Level 1 Fortran variadic definitions
 */

/* Single Precision */

   #define F77_srot(...) F77_srot_base(__VA_ARGS__)
   #define F77_srotg(...) F77_srotg_base(__VA_ARGS__)
   #define F77_srotm(...) F77_srotm_base(__VA_ARGS__)
   #define F77_srotmg(...) F77_srotmg_base(__VA_ARGS__)
   #define F77_sswap(...) F77_sswap_base(__VA_ARGS__)
   #define F77_scopy(...) F77_scopy_base(__VA_ARGS__)
   #define F77_saxpy(...) F77_saxpy_base(__VA_ARGS__)
   #define F77_sdot_sub(...) F77_sdot_sub_base(__VA_ARGS__)
   #define F77_sdsdot_sub(...) F77_sdsdot_sub_base(__VA_ARGS__)
   #define F77_sscal(...) F77_sscal_base(__VA_ARGS__)
   #define F77_snrm2_sub(...) F77_snrm2_sub_base(__VA_ARGS__)
   #define F77_sasum_sub(...) F77_sasum_sub_base(__VA_ARGS__)
   #define F77_isamax_sub(...) F77_isamax_sub_base(__VA_ARGS__)

/* Double Precision */

   #define F77_drot(...) F77_drot_base(__VA_ARGS__)
   #define F77_drotg(...) F77_drotg_base(__VA_ARGS__)
   #define F77_drotm(...) F77_drotm_base(__VA_ARGS__)
   #define F77_drotmg(...) F77_drotmg_base(__VA_ARGS__)
   #define F77_dswap(...) F77_dswap_base(__VA_ARGS__)
   #define F77_dcopy(...) F77_dcopy_base(__VA_ARGS__)
   #define F77_daxpy(...) F77_daxpy_base(__VA_ARGS__)
   #define F77_dswap(...) F77_dswap_base(__VA_ARGS__)
   #define F77_dsdot_sub(...) F77_dsdot_sub_base(__VA_ARGS__)
   #define F77_ddot_sub(...) F77_ddot_sub_base(__VA_ARGS__)
   #define F77_dscal(...) F77_dscal_base(__VA_ARGS__)
   #define F77_dnrm2_sub(...) F77_dnrm2_sub_base(__VA_ARGS__)
   #define F77_dasum_sub(...) F77_dasum_sub_base(__VA_ARGS__)
   #define F77_idamax_sub(...) F77_idamax_sub_base(__VA_ARGS__)

/* Single Complex Precision */

   #define F77_cswap(...) F77_cswap_base(__VA_ARGS__)
   #define F77_ccopy(...) F77_ccopy_base(__VA_ARGS__)
   #define F77_caxpy(...) F77_caxpy_base(__VA_ARGS__)
   #define F77_cswap(...) F77_cswap_base(__VA_ARGS__)
   #define F77_cdotc_sub(...) F77_cdotc_sub_base(__VA_ARGS__)
   #define F77_cdotu_sub(...) F77_cdotu_sub_base(__VA_ARGS__)
   #define F77_cscal(...) F77_cscal_base(__VA_ARGS__)
   #define F77_icamax_sub(...) F77_icamax_sub_base(__VA_ARGS__)
   #define F77_csscal(...) F77_csscal_base(__VA_ARGS__)
   #define F77_scnrm2_sub(...) F77_scnrm2_sub_base(__VA_ARGS__)
   #define F77_scasum_sub(...) F77_scasum_sub_base(__VA_ARGS__)

/* Double Complex Precision */

   #define F77_zswap(...) F77_zswap_base(__VA_ARGS__)
   #define F77_zcopy(...) F77_zcopy_base(__VA_ARGS__)
   #define F77_zaxpy(...) F77_zaxpy_base(__VA_ARGS__)
   #define F77_zswap(...) F77_zswap_base(__VA_ARGS__)
   #define F77_zdotc_sub(...) F77_zdotc_sub_base(__VA_ARGS__)
   #define F77_zdotu_sub(...) F77_zdotu_sub_base(__VA_ARGS__)
   #define F77_zdscal(...) F77_zdscal_base(__VA_ARGS__)
   #define F77_zscal(...) F77_zscal_base(__VA_ARGS__)
   #define F77_dznrm2_sub(...) F77_dznrm2_sub_base(__VA_ARGS__)
   #define F77_dzasum_sub(...) F77_dzasum_sub_base(__VA_ARGS__)
   #define F77_izamax_sub(...) F77_izamax_sub_base(__VA_ARGS__)

/*
 * Level 2 Fortran variadic definitions without FCHAR
 */

   #define F77_sger(...) F77_sger_base(__VA_ARGS__)
   #define F77_dger(...) F77_dger_base(__VA_ARGS__)
   #define F77_cgerc(...) F77_cgerc_base(__VA_ARGS__)
   #define F77_cgeru(...) F77_cgeru_base(__VA_ARGS__)
   #define F77_zgerc(...) F77_zgerc_base(__VA_ARGS__)
   #define F77_zgeru(...) F77_zgeru_base(__VA_ARGS__)

#ifdef BLAS_FORTRAN_STRLEN_END

   /*
   * Level 2 Fortran variadic definitions with BLAS_FORTRAN_STRLEN_END
   */

   /* Single Precision */

   #define F77_sgemv(...) F77_sgemv_base(__VA_ARGS__, 1)
   #define F77_sgbmv(...) F77_sgbmv_base(__VA_ARGS__, 1)
   #define F77_ssymv(...) F77_ssymv_base(__VA_ARGS__, 1)
   #define F77_ssbmv(...) F77_ssbmv_base(__VA_ARGS__, 1)
   #define F77_sspmv(...) F77_sspmv_base(__VA_ARGS__, 1)
   #define F77_strmv(...) F77_strmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_stbmv(...) F77_stbmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_strsv(...) F77_strsv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_stbsv(...) F77_stbsv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_stpmv(...) F77_stpmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_stpsv(...) F77_stpsv_base(__VA_ARGS__, 1, 1, 1)
      #define F77_ssyr(...) F77_ssyr_base(__VA_ARGS__, 1)
   #define F77_sspr(...) F77_sspr_base(__VA_ARGS__, 1)
   #define F77_sspr2(...) F77_sspr2_base(__VA_ARGS__, 1)
   #define F77_ssyr2(...) F77_ssyr2_base(__VA_ARGS__, 1)

   /* Double Precision */

   #define F77_dgemv(...) F77_dgemv_base(__VA_ARGS__, 1)
   #define F77_dgbmv(...) F77_dgbmv_base(__VA_ARGS__, 1)
   #define F77_dsymv(...) F77_dsymv_base(__VA_ARGS__, 1)
   #define F77_dsbmv(...) F77_dsbmv_base(__VA_ARGS__, 1)
   #define F77_dspmv(...) F77_dspmv_base(__VA_ARGS__, 1)
   #define F77_dtrmv(...) F77_dtrmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_dtbmv(...) F77_dtbmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_dtrsv(...) F77_dtrsv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_dtbsv(...) F77_dtbsv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_dtpmv(...) F77_dtpmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_dtpsv(...) F77_dtpsv_base(__VA_ARGS__, 1, 1, 1)
      #define F77_dsyr(...) F77_dsyr_base(__VA_ARGS__, 1)
   #define F77_dspr(...) F77_dspr_base(__VA_ARGS__, 1)
   #define F77_dspr2(...) F77_dspr2_base(__VA_ARGS__, 1)
   #define F77_dsyr2(...) F77_dsyr2_base(__VA_ARGS__, 1)

      /* Single Complex Precision */

   #define F77_cgemv(...) F77_cgemv_base(__VA_ARGS__, 1)
   #define F77_cgbmv(...) F77_cgbmv_base(__VA_ARGS__, 1)
   #define F77_chemv(...) F77_chemv_base(__VA_ARGS__, 1)
   #define F77_chbmv(...) F77_chbmv_base(__VA_ARGS__, 1)
   #define F77_chpmv(...) F77_chpmv_base(__VA_ARGS__, 1)
   #define F77_ctrmv(...) F77_ctrmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ctbmv(...) F77_ctbmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ctpmv(...) F77_ctpmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ctrsv(...) F77_ctrsv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ctbsv(...) F77_ctbsv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ctpsv(...) F77_ctpsv_base(__VA_ARGS__, 1, 1, 1)
         #define F77_cher(...) F77_cher_base(__VA_ARGS__, 1)
   #define F77_cher2(...) F77_cher2_base(__VA_ARGS__, 1)
   #define F77_chpr(...) F77_chpr_base(__VA_ARGS__, 1)
   #define F77_chpr2(...) F77_chpr2_base(__VA_ARGS__, 1)

   /* Double Complex Precision */

   #define F77_zgemv(...) F77_zgemv_base(__VA_ARGS__, 1)
   #define F77_zgbmv(...) F77_zgbmv_base(__VA_ARGS__, 1)
   #define F77_zhemv(...) F77_zhemv_base(__VA_ARGS__, 1)
   #define F77_zhbmv(...) F77_zhbmv_base(__VA_ARGS__, 1)
   #define F77_zhpmv(...) F77_zhpmv_base(__VA_ARGS__, 1)
   #define F77_ztrmv(...) F77_ztrmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ztbmv(...) F77_ztbmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ztpmv(...) F77_ztpmv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ztrsv(...) F77_ztrsv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ztbsv(...) F77_ztbsv_base(__VA_ARGS__, 1, 1, 1)
   #define F77_ztpsv(...) F77_ztpsv_base(__VA_ARGS__, 1, 1, 1)
         #define F77_zher(...) F77_zher_base(__VA_ARGS__, 1)
   #define F77_zher2(...) F77_zher2_base(__VA_ARGS__, 1)
   #define F77_zhpr(...) F77_zhpr_base(__VA_ARGS__, 1)
   #define F77_zhpr2(...) F77_zhpr2_base(__VA_ARGS__, 1)

   /*
   * Level 3 Fortran variadic definitions with BLAS_FORTRAN_STRLEN_END
   */

   /* Single Precision */

   #define F77_sgemm(...) F77_sgemm_base(__VA_ARGS__, 1, 1)
   #define F77_ssymm(...) F77_ssymm_base(__VA_ARGS__, 1, 1)
   #define F77_ssyrk(...) F77_ssyrk_base(__VA_ARGS__, 1, 1)
   #define F77_ssyr2k(...) F77_ssyr2k_base(__VA_ARGS__, 1, 1)
   #define F77_strmm(...) F77_strmm_base(__VA_ARGS__, 1, 1, 1, 1)
   #define F77_strsm(...) F77_strsm_base(__VA_ARGS__, 1, 1, 1, 1)

   /* Double Precision */

   #define F77_dgemm(...) F77_dgemm_base(__VA_ARGS__, 1, 1)
   #define F77_dsymm(...) F77_dsymm_base(__VA_ARGS__, 1, 1)
   #define F77_dsyrk(...) F77_dsyrk_base(__VA_ARGS__, 1, 1)
   #define F77_dsyr2k(...) F77_dsyr2k_base(__VA_ARGS__, 1, 1)
   #define F77_dtrmm(...) F77_dtrmm_base(__VA_ARGS__, 1, 1, 1, 1)
   #define F77_dtrsm(...) F77_dtrsm_base(__VA_ARGS__, 1, 1, 1, 1)

   /* Single Complex Precision */

   #define F77_cgemm(...) F77_cgemm_base(__VA_ARGS__, 1, 1)
   #define F77_csymm(...) F77_csymm_base(__VA_ARGS__, 1, 1)
   #define F77_chemm(...) F77_chemm_base(__VA_ARGS__, 1, 1)
   #define F77_csyrk(...) F77_csyrk_base(__VA_ARGS__, 1, 1)
   #define F77_cherk(...) F77_cherk_base(__VA_ARGS__, 1, 1)
   #define F77_csyr2k(...) F77_csyr2k_base(__VA_ARGS__, 1, 1)
   #define F77_cher2k(...) F77_cher2k_base(__VA_ARGS__, 1, 1)
   #define F77_ctrmm(...) F77_ctrmm_base(__VA_ARGS__, 1, 1, 1, 1)
   #define F77_ctrsm(...) F77_ctrsm_base(__VA_ARGS__, 1, 1, 1, 1)

   /* Double Complex Precision */

   #define F77_zgemm(...) F77_zgemm_base(__VA_ARGS__, 1, 1)
   #define F77_zsymm(...) F77_zsymm_base(__VA_ARGS__, 1, 1)
   #define F77_zhemm(...) F77_zhemm_base(__VA_ARGS__, 1, 1)
   #define F77_zsyrk(...) F77_zsyrk_base(__VA_ARGS__, 1, 1)
   #define F77_zherk(...) F77_zherk_base(__VA_ARGS__, 1, 1)
   #define F77_zsyr2k(...) F77_zsyr2k_base(__VA_ARGS__, 1, 1)
   #define F77_zher2k(...) F77_zher2k_base(__VA_ARGS__, 1, 1)
   #define F77_ztrmm(...) F77_ztrmm_base(__VA_ARGS__, 1, 1, 1, 1)
   #define F77_ztrsm(...) F77_ztrsm_base(__VA_ARGS__, 1, 1, 1, 1)

#else
                              
   /*
   * Level 2 Fortran variadic definitions without BLAS_FORTRAN_STRLEN_END
   */

   /* Single Precision */

   #define F77_sgemv(...) F77_sgemv_base(__VA_ARGS__)
   #define F77_sgbmv(...) F77_sgbmv_base(__VA_ARGS__)
   #define F77_ssymv(...) F77_ssymv_base(__VA_ARGS__)
   #define F77_ssbmv(...) F77_ssbmv_base(__VA_ARGS__)
   #define F77_sspmv(...) F77_sspmv_base(__VA_ARGS__)
   #define F77_strmv(...) F77_strmv_base(__VA_ARGS__)
   #define F77_stbmv(...) F77_stbmv_base(__VA_ARGS__)
   #define F77_strsv(...) F77_strsv_base(__VA_ARGS__)
   #define F77_stbsv(...) F77_stbsv_base(__VA_ARGS__)
   #define F77_stpmv(...) F77_stpmv_base(__VA_ARGS__)
   #define F77_stpsv(...) F77_stpsv_base(__VA_ARGS__)
      #define F77_ssyr(...) F77_ssyr_base(__VA_ARGS__)
   #define F77_sspr(...) F77_sspr_base(__VA_ARGS__)
   #define F77_sspr2(...) F77_sspr2_base(__VA_ARGS__)
   #define F77_ssyr2(...) F77_ssyr2_base(__VA_ARGS__)

   /* Double Precision */

   #define F77_dgemv(...) F77_dgemv_base(__VA_ARGS__)
   #define F77_dgbmv(...) F77_dgbmv_base(__VA_ARGS__)
   #define F77_dsymv(...) F77_dsymv_base(__VA_ARGS__)
   #define F77_dsbmv(...) F77_dsbmv_base(__VA_ARGS__)
   #define F77_dspmv(...) F77_dspmv_base(__VA_ARGS__)
   #define F77_dtrmv(...) F77_dtrmv_base(__VA_ARGS__)
   #define F77_dtbmv(...) F77_dtbmv_base(__VA_ARGS__)
   #define F77_dtrsv(...) F77_dtrsv_base(__VA_ARGS__)
   #define F77_dtbsv(...) F77_dtbsv_base(__VA_ARGS__)
   #define F77_dtpmv(...) F77_dtpmv_base(__VA_ARGS__)
   #define F77_dtpsv(...) F77_dtpsv_base(__VA_ARGS__)
      #define F77_dsyr(...) F77_dsyr_base(__VA_ARGS__)
   #define F77_dspr(...) F77_dspr_base(__VA_ARGS__)
   #define F77_dspr2(...) F77_dspr2_base(__VA_ARGS__)
   #define F77_dsyr2(...) F77_dsyr2_base(__VA_ARGS__)

   /* Single Complex Precision */

   #define F77_cgemv(...) F77_cgemv_base(__VA_ARGS__)
   #define F77_cgbmv(...) F77_cgbmv_base(__VA_ARGS__)
   #define F77_chemv(...) F77_chemv_base(__VA_ARGS__)
   #define F77_chbmv(...) F77_chbmv_base(__VA_ARGS__)
   #define F77_chpmv(...) F77_chpmv_base(__VA_ARGS__)
   #define F77_ctrmv(...) F77_ctrmv_base(__VA_ARGS__)
   #define F77_ctbmv(...) F77_ctbmv_base(__VA_ARGS__)
   #define F77_ctpmv(...) F77_ctpmv_base(__VA_ARGS__)
   #define F77_ctrsv(...) F77_ctrsv_base(__VA_ARGS__)
   #define F77_ctbsv(...) F77_ctbsv_base(__VA_ARGS__)
   #define F77_ctpsv(...) F77_ctpsv_base(__VA_ARGS__)
         #define F77_cher(...) F77_cher_base(__VA_ARGS__)
   #define F77_cher2(...) F77_cher2_base(__VA_ARGS__)
   #define F77_chpr(...) F77_chpr_base(__VA_ARGS__)
   #define F77_chpr2(...) F77_chpr2_base(__VA_ARGS__)

   /* Double Complex Precision */

   #define F77_zgemv(...) F77_zgemv_base(__VA_ARGS__)
   #define F77_zgbmv(...) F77_zgbmv_base(__VA_ARGS__)
   #define F77_zhemv(...) F77_zhemv_base(__VA_ARGS__)
   #define F77_zhbmv(...) F77_zhbmv_base(__VA_ARGS__)
   #define F77_zhpmv(...) F77_zhpmv_base(__VA_ARGS__)
   #define F77_ztrmv(...) F77_ztrmv_base(__VA_ARGS__)
   #define F77_ztbmv(...) F77_ztbmv_base(__VA_ARGS__)
   #define F77_ztpmv(...) F77_ztpmv_base(__VA_ARGS__)
   #define F77_ztrsv(...) F77_ztrsv_base(__VA_ARGS__)
   #define F77_ztbsv(...) F77_ztbsv_base(__VA_ARGS__)
   #define F77_ztpsv(...) F77_ztpsv_base(__VA_ARGS__)
         #define F77_zher(...) F77_zher_base(__VA_ARGS__)
   #define F77_zher2(...) F77_zher2_base(__VA_ARGS__)
   #define F77_zhpr(...) F77_zhpr_base(__VA_ARGS__)
   #define F77_zhpr2(...) F77_zhpr2_base(__VA_ARGS__)

   /*
   * Level 3 Fortran variadic definitions without BLAS_FORTRAN_STRLEN_END
   */

   /* Single Precision */

   #define F77_sgemm(...) F77_sgemm_base(__VA_ARGS__)
   #define F77_ssymm(...) F77_ssymm_base(__VA_ARGS__)
   #define F77_ssyrk(...) F77_ssyrk_base(__VA_ARGS__)
   #define F77_ssyr2k(...) F77_ssyr2k_base(__VA_ARGS__)
   #define F77_strmm(...) F77_strmm_base(__VA_ARGS__)
   #define F77_strsm(...) F77_strsm_base(__VA_ARGS__)

   /* Double Precision */

   #define F77_dgemm(...) F77_dgemm_base(__VA_ARGS__)
   #define F77_dsymm(...) F77_dsymm_base(__VA_ARGS__)
   #define F77_dsyrk(...) F77_dsyrk_base(__VA_ARGS__)
   #define F77_dsyr2k(...) F77_dsyr2k_base(__VA_ARGS__)
   #define F77_dtrmm(...) F77_dtrmm_base(__VA_ARGS__)
   #define F77_dtrsm(...) F77_dtrsm_base(__VA_ARGS__)

   /* Single Complex Precision */

   #define F77_cgemm(...) F77_cgemm_base(__VA_ARGS__)
   #define F77_csymm(...) F77_csymm_base(__VA_ARGS__)
   #define F77_chemm(...) F77_chemm_base(__VA_ARGS__)
   #define F77_csyrk(...) F77_csyrk_base(__VA_ARGS__)
   #define F77_cherk(...) F77_cherk_base(__VA_ARGS__)
   #define F77_csyr2k(...) F77_csyr2k_base(__VA_ARGS__)
   #define F77_cher2k(...) F77_cher2k_base(__VA_ARGS__)
   #define F77_ctrmm(...) F77_ctrmm_base(__VA_ARGS__)
   #define F77_ctrsm(...) F77_ctrsm_base(__VA_ARGS__)

   /* Double Complex Precision */

   #define F77_zgemm(...) F77_zgemm_base(__VA_ARGS__)
   #define F77_zsymm(...) F77_zsymm_base(__VA_ARGS__)
   #define F77_zhemm(...) F77_zhemm_base(__VA_ARGS__)
   #define F77_zsyrk(...) F77_zsyrk_base(__VA_ARGS__)
   #define F77_zherk(...) F77_zherk_base(__VA_ARGS__)
   #define F77_zsyr2k(...) F77_zsyr2k_base(__VA_ARGS__)
   #define F77_zher2k(...) F77_zher2k_base(__VA_ARGS__)
   #define F77_ztrmm(...) F77_ztrmm_base(__VA_ARGS__)
   #define F77_ztrsm(...) F77_ztrsm_base(__VA_ARGS__)

#endif

/*
 * Base function prototypes
 */

#ifdef __cplusplus
extern "C" {
#endif

void F77_xerbla(FCHAR, void *);
void F77_xerbla_base(FCHAR, void *
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
/*
 * Level 1 Fortran Prototypes
 */

/* Single Precision */

   void F77_srot_base(FINT, float *, FINT, float *, FINT, const float *, const float *);
   void F77_srotg_base(float *,float *,float *,float *);
   void F77_srotm_base( FINT, float *, FINT, float *, FINT, const float *);
   void F77_srotmg_base(float *,float *,float *,const float *, float *);
   void F77_sswap_base( FINT, float *, FINT, float *, FINT);
   void F77_scopy_base( FINT, const float *, FINT, float *, FINT);
   void F77_saxpy_base( FINT, const float *, const float *, FINT, float *, FINT);
   void F77_sdot_sub_base(FINT, const float *, FINT, const float *, FINT, float *);
   void F77_sdsdot_sub_base( FINT, const float *, const float *, FINT, const float *, FINT, float *);
   void F77_sscal_base( FINT, const float *, float *, FINT);
   void F77_snrm2_sub_base( FINT, const float *, FINT, float *);
   void F77_sasum_sub_base( FINT, const float *, FINT, float *);
   void F77_isamax_sub_base( FINT, const float * , FINT, FINT2);

/* Double Precision */

   void F77_drot_base(FINT, double *, FINT, double *, FINT, const double *, const double *);
   void F77_drotg_base(double *,double *,double *,double *);
   void F77_drotm_base( FINT, double *, FINT, double *, FINT, const double *);
   void F77_drotmg_base(double *,double *,double *,const double *, double *);
   void F77_dswap_base( FINT, double *, FINT, double *, FINT);
   void F77_dcopy_base( FINT, const double *, FINT, double *, FINT);
   void F77_daxpy_base( FINT, const double *, const double *, FINT, double *, FINT);
   void F77_dswap_base( FINT, double *, FINT, double *, FINT);
   void F77_dsdot_sub_base(FINT, const float *, FINT, const float *, FINT, double *);
   void F77_ddot_sub_base( FINT, const double *, FINT, const double *, FINT, double *);
   void F77_dscal_base( FINT, const double *, double *, FINT);
   void F77_dnrm2_sub_base( FINT, const double *, FINT, double *);
   void F77_dasum_sub_base( FINT, const double *, FINT, double *);
   void F77_idamax_sub_base( FINT, const double * , FINT, FINT2);

/* Single Complex Precision */

   void F77_cswap_base( FINT, void *, FINT, void *, FINT);
   void F77_ccopy_base( FINT, const void *, FINT, void *, FINT);
   void F77_caxpy_base( FINT, const void *, const void *, FINT, void *, FINT);
   void F77_cswap_base( FINT, void *, FINT, void *, FINT);
   void F77_cdotc_sub_base( FINT, const void *, FINT, const void *, FINT, void *);
   void F77_cdotu_sub_base( FINT, const void *, FINT, const void *, FINT, void *);
   void F77_cscal_base( FINT, const void *, void *, FINT);
   void F77_icamax_sub_base( FINT, const void *, FINT, FINT2);
   void F77_csscal_base( FINT, const float *, void *, FINT);
   void F77_scnrm2_sub_base( FINT, const void *, FINT, float *);
   void F77_scasum_sub_base( FINT, const void *, FINT, float *);

/* Double Complex Precision */

   void F77_zswap_base( FINT, void *, FINT, void *, FINT);
   void F77_zcopy_base( FINT, const void *, FINT, void *, FINT);
   void F77_zaxpy_base( FINT, const void *, const void *, FINT, void *, FINT);
   void F77_zswap_base( FINT, void *, FINT, void *, FINT);
   void F77_zdotc_sub_base( FINT, const void *, FINT, const void *, FINT, void *);
   void F77_zdotu_sub_base( FINT, const void *, FINT, const void *, FINT, void *);
   void F77_zdscal_base( FINT, const double *, void *, FINT);
   void F77_zscal_base( FINT, const void *, void *, FINT);
   void F77_dznrm2_sub_base( FINT, const void *, FINT, double *);
   void F77_dzasum_sub_base( FINT, const void *, FINT, double *);
   void F77_izamax_sub_base( FINT, const void *, FINT, FINT2);

/*
 * Level 2 Fortran Prototypes
 */

/* Single Precision */

   void F77_sgemv_base(FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_sgbmv_base(FCHAR, FINT, FINT, FINT, FINT, const float *,  const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_ssymv_base(FCHAR, FINT, const float *, const float *, FINT, const float *,  FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_ssbmv_base(FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_sspmv_base(FCHAR, FINT, const float *, const float *, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_strmv_base(FCHAR, FCHAR, FCHAR, FINT, const float *, FINT, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_stbmv_base(FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, FINT, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_strsv_base(FCHAR, FCHAR, FCHAR, FINT, const float *, FINT, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_stbsv_base(FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, FINT, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_stpmv_base(FCHAR, FCHAR, FCHAR, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_stpsv_base(FCHAR, FCHAR, FCHAR, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_sger_base( FINT, FINT, const float *, const float *, FINT, const float *, FINT, float *, FINT);
   void F77_ssyr_base(FCHAR, FINT, const float *, const float *, FINT, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_sspr_base(FCHAR, FINT, const float *, const float *, FINT, float *
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_sspr2_base(FCHAR, FINT, const float *, const float *, FINT, const float *, FINT,  float *
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_ssyr2_base(FCHAR, FINT, const float *, const float *, FINT, const float *, FINT,  float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );

/* Double Precision */

   void F77_dgemv_base(FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_dgbmv_base(FCHAR, FINT, FINT, FINT, FINT, const double *,  const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_dsymv_base(FCHAR, FINT, const double *, const double *, FINT, const double *,  FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_dsbmv_base(FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_dspmv_base(FCHAR, FINT, const double *, const double *, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_dtrmv_base(FCHAR, FCHAR, FCHAR, FINT, const double *, FINT, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_dtbmv_base(FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, FINT, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_dtrsv_base(FCHAR, FCHAR, FCHAR, FINT, const double *, FINT, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_dtbsv_base(FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, FINT, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_dtpmv_base(FCHAR, FCHAR, FCHAR, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_dtpsv_base(FCHAR, FCHAR, FCHAR, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_dger_base( FINT, FINT, const double *, const double *, FINT, const double *, FINT, double *, FINT);
   void F77_dsyr_base(FCHAR, FINT, const double *, const double *, FINT, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_dspr_base(FCHAR, FINT, const double *, const double *, FINT, double *
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_dspr2_base(FCHAR, FINT, const double *, const double *, FINT, const double *, FINT,  double *
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_dsyr2_base(FCHAR, FINT, const double *, const double *, FINT, const double *, FINT,  double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );

/* Single Complex Precision */

   void F77_cgemv_base(FCHAR, FINT, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_cgbmv_base(FCHAR, FINT, FINT, FINT, FINT, const void *,  const void *, FINT, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_chemv_base(FCHAR, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_chbmv_base(FCHAR, FINT, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_chpmv_base(FCHAR, FINT, const void *, const void *, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_ctrmv_base(FCHAR, FCHAR, FCHAR, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ctbmv_base(FCHAR, FCHAR, FCHAR, FINT, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ctpmv_base(FCHAR, FCHAR, FCHAR, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ctrsv_base(FCHAR, FCHAR, FCHAR, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ctbsv_base(FCHAR, FCHAR, FCHAR, FINT, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ctpsv_base(FCHAR, FCHAR, FCHAR, FINT, const void *, void *,FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_cgerc_base( FINT, FINT, const void *, const void *, FINT, const void *, FINT, void *, FINT);
   void F77_cgeru_base( FINT, FINT, const void *, const void *, FINT, const void *, FINT, void *,  FINT);
   void F77_cher_base(FCHAR, FINT, const float *, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_cher2_base(FCHAR, FINT, const void *, const void *, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_chpr_base(FCHAR, FINT, const float *, const void *, FINT, void *
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_chpr2_base(FCHAR, FINT, const float *, const void *, FINT, const void *, FINT, void *
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );

/* Double Complex Precision */

   void F77_zgemv_base(FCHAR, FINT, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_zgbmv_base(FCHAR, FINT, FINT, FINT, FINT, const void *,  const void *, FINT, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_zhemv_base(FCHAR, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_zhbmv_base(FCHAR, FINT, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_zhpmv_base(FCHAR, FINT, const void *, const void *, const void *, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_ztrmv_base(FCHAR, FCHAR, FCHAR, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ztbmv_base(FCHAR, FCHAR, FCHAR, FINT, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ztpmv_base(FCHAR, FCHAR, FCHAR, FINT, const void *, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ztrsv_base(FCHAR, FCHAR, FCHAR, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ztbsv_base(FCHAR, FCHAR, FCHAR, FINT, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_ztpsv_base(FCHAR, FCHAR, FCHAR, FINT, const void *, void *,FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t
   #endif
   );
   void F77_zgerc_base( FINT, FINT, const void *, const void *, FINT, const void *, FINT, void *, FINT);
   void F77_zgeru_base( FINT, FINT, const void *, const void *, FINT, const void *, FINT, void *,  FINT);
   void F77_zher_base(FCHAR, FINT, const double *, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_zher2_base(FCHAR, FINT, const void *, const void *, FINT, const void *, FINT, void *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_zhpr_base(FCHAR, FINT, const double *, const void *, FINT, void *
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );
   void F77_zhpr2_base(FCHAR, FINT, const double *, const void *, FINT, const void *, FINT, void *
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t
   #endif
   );

/*
 * Level 3 Fortran Prototypes
 */

/* Single Precision */

   void F77_sgemm_base(FCHAR, FCHAR, FINT, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_ssymm_base(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_ssyrk_base(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_ssyr2k_base(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_strmm_base(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t, size_t
   #endif
   );
   void F77_strsm_base(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t, size_t
   #endif
   );

/* Double Precision */

   void F77_dgemm_base(FCHAR, FCHAR, FINT, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_dsymm_base(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_dsyrk_base(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_dsyr2k_base(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_dtrmm_base(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t, size_t
   #endif
   );
   void F77_dtrsm_base(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t, size_t
   #endif
   );

/* Single Complex Precision */

   void F77_cgemm_base(FCHAR, FCHAR, FINT, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_csymm_base(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_chemm_base(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_csyrk_base(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_cherk_base(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_csyr2k_base(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_cher2k_base(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_ctrmm_base(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t, size_t
   #endif
   );
   void F77_ctrsm_base(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, float *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t, size_t
   #endif
   );

/* Double Complex Precision */

   void F77_zgemm_base(FCHAR, FCHAR, FINT, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_zsymm_base(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_zhemm_base(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_zsyrk_base(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_zherk_base(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_zsyr2k_base(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_zher2k_base(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t
   #endif
   );
   void F77_ztrmm_base(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t, size_t
   #endif
   );
   void F77_ztrsm_base(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT
   #ifdef BLAS_FORTRAN_STRLEN_END
      , size_t, size_t, size_t, size_t
   #endif
   );

#ifdef __cplusplus
}
#endif

#endif /*  CBLAS_F77_H */
