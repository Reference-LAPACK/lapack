/*
 * cblas_test.h
 * Written by Keita Teranishi
 */
#ifndef CBLAS_TEST_H
#define CBLAS_TEST_H
#include "cblas.h"
#include "cblas_globals.h"
#include "cblas_mangling.h"

#ifndef F77_GLOBAL_SUFFIX
#define F77_GLOBAL_SUFFIX(a,b) F77_GLOBAL_SUFFIX_(API_SUFFIX(a),API_SUFFIX(b))
#define F77_GLOBAL_SUFFIX_(a,b) F77_GLOBAL(a,b)
#endif

/* It seems all current Fortran compilers put strlen at end.
*  Some historical compilers put strlen after the str argument
*  or make the str argument into a struct. */
#define BLAS_FORTRAN_STRLEN_END

#ifndef FORTRAN_STRLEN
  #define FORTRAN_STRLEN size_t
#endif

#define  TRUE           1
#define  PASSED         1
#define  TEST_ROW_MJR	1

#define  FALSE          0
#define  FAILED         0
#define  TEST_COL_MJR	0

#define  INVALID       -1
#define  UNDEFINED     -1

typedef struct { float real; float imag; } CBLAS_TEST_COMPLEX;
typedef struct { double real; double imag; } CBLAS_TEST_ZOMPLEX;

#define F77_xerbla 		F77_GLOBAL_SUFFIX(xerbla,XERBLA)
/*
 * Level 1 BLAS
 */
#define F77_srotg 		F77_GLOBAL_SUFFIX(srotgtest,SROTGTEST)
#define F77_srotmg 		F77_GLOBAL_SUFFIX(srotmgtest,SROTMGTEST)
#define F77_srot 		F77_GLOBAL_SUFFIX(srottest,SROTTEST)
#define F77_srotm 		F77_GLOBAL_SUFFIX(srotmtest,SROTMTEST)
#define F77_drotg 		F77_GLOBAL_SUFFIX(drotgtest,DROTGTEST)
#define F77_drotmg 		F77_GLOBAL_SUFFIX(drotmgtest,DROTMGTEST)
#define F77_drot 		F77_GLOBAL_SUFFIX(drottest,DROTTEST)
#define F77_drotm 		F77_GLOBAL_SUFFIX(drotmtest,DROTMTEST)
#define F77_sswap 		F77_GLOBAL_SUFFIX(sswaptest,SSWAPTEST)
#define F77_scopy 		F77_GLOBAL_SUFFIX(scopytest,SCOPYTEST)
#define F77_saxpy 		F77_GLOBAL_SUFFIX(saxpytest,SAXPYTEST)
#define F77_saxpby 		F77_GLOBAL_SUFFIX(saxpbytest,SAXPBYTEST)
#define F77_isamax 		F77_GLOBAL_SUFFIX(isamaxtest,ISAMAXTEST)
#define F77_dswap 		F77_GLOBAL_SUFFIX(dswaptest,DSWAPTEST)
#define F77_dcopy 		F77_GLOBAL_SUFFIX(dcopytest,DCOPYTEST)
#define F77_daxpy 		F77_GLOBAL_SUFFIX(daxpytest,DAXPYTEST)
#define F77_daxpby 		F77_GLOBAL_SUFFIX(daxpbytest,DAXPBYTEST)
#define F77_idamax 		F77_GLOBAL_SUFFIX(idamaxtest,IDAMAXTEST)
#define F77_cswap 		F77_GLOBAL_SUFFIX(cswaptest,CSWAPTEST)
#define F77_ccopy 		F77_GLOBAL_SUFFIX(ccopytest,CCOPYTEST)
#define F77_caxpy 		F77_GLOBAL_SUFFIX(caxpytest,CAXPYTEST)
#define F77_caxpby 		F77_GLOBAL_SUFFIX(caxpbytest,CAXPBYTEST)
#define F77_icamax 		F77_GLOBAL_SUFFIX(icamaxtest,ICAMAXTEST)
#define F77_zswap 		F77_GLOBAL_SUFFIX(zswaptest,ZSWAPTEST)
#define F77_zcopy 		F77_GLOBAL_SUFFIX(zcopytest,ZCOPYTEST)
#define F77_zaxpy 		F77_GLOBAL_SUFFIX(zaxpytest,ZAXPYTEST)
#define F77_zaxpby 		F77_GLOBAL_SUFFIX(zaxpbytest,ZAXPBYTEST)
#define F77_izamax 		F77_GLOBAL_SUFFIX(izamaxtest,IZAMAXTEST)
#define F77_sdot 		F77_GLOBAL_SUFFIX(sdottest,SDOTTEST)
#define F77_ddot 		F77_GLOBAL_SUFFIX(ddottest,DDOTTEST)
#define F77_dsdot 		F77_GLOBAL_SUFFIX(dsdottest,DSDOTTEST)
#define F77_sscal 		F77_GLOBAL_SUFFIX(sscaltest,SSCALTEST)
#define F77_dscal 		F77_GLOBAL_SUFFIX(dscaltest,DSCALTEST)
#define F77_cscal 		F77_GLOBAL_SUFFIX(cscaltest,CSCALTEST)
#define F77_zscal 		F77_GLOBAL_SUFFIX(zscaltest,ZSCALTEST)
#define F77_csscal 		F77_GLOBAL_SUFFIX(csscaltest,CSSCALTEST)
#define F77_zdscal 		F77_GLOBAL_SUFFIX(zdscaltest,ZDSCALTEST)
#define F77_cdotu 		F77_GLOBAL_SUFFIX(cdotutest,CDOTUTEST)
#define F77_cdotc 		F77_GLOBAL_SUFFIX(cdotctest,CDOTCTEST)
#define F77_zdotu 		F77_GLOBAL_SUFFIX(zdotutest,ZDOTUTEST)
#define F77_zdotc 		F77_GLOBAL_SUFFIX(zdotctest,ZDOTCTEST)
#define F77_snrm2 		F77_GLOBAL_SUFFIX(snrm2test,SNRM2TEST)
#define F77_sasum 		F77_GLOBAL_SUFFIX(sasumtest,SASUMTEST)
#define F77_dnrm2 		F77_GLOBAL_SUFFIX(dnrm2test,DNRM2TEST)
#define F77_dasum 		F77_GLOBAL_SUFFIX(dasumtest,DASUMTEST)
#define F77_scnrm2 		F77_GLOBAL_SUFFIX(scnrm2test,SCNRM2TEST)
#define F77_scasum 		F77_GLOBAL_SUFFIX(scasumtest,SCASUMTEST)
#define F77_dznrm2 		F77_GLOBAL_SUFFIX(dznrm2test,DZNRM2TEST)
#define F77_dzasum 		F77_GLOBAL_SUFFIX(dzasumtest,DZASUMTEST)
#define F77_sdsdot 		F77_GLOBAL_SUFFIX(sdsdottest, SDSDOTTEST)
/*
 * Level 2 BLAS
 */
#define F77_s2chke 		F77_GLOBAL_SUFFIX(cs2chke,CS2CHKE)
#define F77_d2chke 		F77_GLOBAL_SUFFIX(cd2chke,CD2CHKE)
#define F77_c2chke 		F77_GLOBAL_SUFFIX(cc2chke,CC2CHKE)
#define F77_z2chke 		F77_GLOBAL_SUFFIX(cz2chke,CZ2CHKE)
#define F77_ssymv 		F77_GLOBAL_SUFFIX(cssymv,CSSYMV)
#define F77_ssbmv 		F77_GLOBAL_SUFFIX(cssbmv,CSSBMV)
#define F77_sspmv 		F77_GLOBAL_SUFFIX(csspmv,CSSPMV)
#define F77_sskewsymv F77_GLOBAL_SUFFIX(csskewsymv,CSSKEWSYMV)
#define F77_sger 		F77_GLOBAL_SUFFIX(csger,CSGER)
#define F77_ssyr 		F77_GLOBAL_SUFFIX(cssyr,CSSYR)
#define F77_sspr 		F77_GLOBAL_SUFFIX(csspr,CSSPR)
#define F77_ssyr2 		F77_GLOBAL_SUFFIX(cssyr2,CSSYR2)
#define F77_sspr2 		F77_GLOBAL_SUFFIX(csspr2,CSSPR2)
#define F77_sskewsyr2 F77_GLOBAL_SUFFIX(csskewsyr2,CSSKEWSYR2)
#define F77_dsymv 		F77_GLOBAL_SUFFIX(cdsymv,CDSYMV)
#define F77_dsbmv 		F77_GLOBAL_SUFFIX(cdsbmv,CDSBMV)
#define F77_dspmv 		F77_GLOBAL_SUFFIX(cdspmv,CDSPMV)
#define F77_dskewsymv F77_GLOBAL_SUFFIX(cdskewsymv,CDSKEWSYMV)
#define F77_dger 		F77_GLOBAL_SUFFIX(cdger,CDGER)
#define F77_dsyr 		F77_GLOBAL_SUFFIX(cdsyr,CDSYR)
#define F77_dspr 		F77_GLOBAL_SUFFIX(cdspr,CDSPR)
#define F77_dsyr2 		F77_GLOBAL_SUFFIX(cdsyr2,CDSYR2)
#define F77_dspr2 		F77_GLOBAL_SUFFIX(cdspr2,CDSPR2)
#define F77_dskewsyr2 F77_GLOBAL_SUFFIX(cdskewsyr2,CDSKEWSYR2)
#define F77_chemv 		F77_GLOBAL_SUFFIX(cchemv,CCHEMV)
#define F77_chbmv 		F77_GLOBAL_SUFFIX(cchbmv,CCHBMV)
#define F77_chpmv 		F77_GLOBAL_SUFFIX(cchpmv,CCHPMV)
#define F77_cgeru 		F77_GLOBAL_SUFFIX(ccgeru,CCGERU)
#define F77_cgerc 		F77_GLOBAL_SUFFIX(ccgerc,CCGERC)
#define F77_cher 		F77_GLOBAL_SUFFIX(ccher,CCHER)
#define F77_chpr 		F77_GLOBAL_SUFFIX(cchpr,CCHPR)
#define F77_cher2 		F77_GLOBAL_SUFFIX(ccher2,CCHER2)
#define F77_chpr2 		F77_GLOBAL_SUFFIX(cchpr2,CCHPR2)
#define F77_zhemv 		F77_GLOBAL_SUFFIX(czhemv,CZHEMV)
#define F77_zhbmv 		F77_GLOBAL_SUFFIX(czhbmv,CZHBMV)
#define F77_zhpmv 		F77_GLOBAL_SUFFIX(czhpmv,CZHPMV)
#define F77_zgeru 		F77_GLOBAL_SUFFIX(czgeru,CZGERU)
#define F77_zgerc 		F77_GLOBAL_SUFFIX(czgerc,CZGERC)
#define F77_zher 		F77_GLOBAL_SUFFIX(czher,CZHER)
#define F77_zhpr 		F77_GLOBAL_SUFFIX(czhpr,CZHPR)
#define F77_zher2 		F77_GLOBAL_SUFFIX(czher2,CZHER2)
#define F77_zhpr2 		F77_GLOBAL_SUFFIX(czhpr2,CZHPR2)
#define F77_sgemv 		F77_GLOBAL_SUFFIX(csgemv,CSGEMV)
#define F77_sgbmv 		F77_GLOBAL_SUFFIX(csgbmv,CSGBMV)
#define F77_strmv 		F77_GLOBAL_SUFFIX(cstrmv,CSTRMV)
#define F77_stbmv 		F77_GLOBAL_SUFFIX(cstbmv,CSTBMV)
#define F77_stpmv 		F77_GLOBAL_SUFFIX(cstpmv,CSTPMV)
#define F77_strsv 		F77_GLOBAL_SUFFIX(cstrsv,CSTRSV)
#define F77_stbsv 		F77_GLOBAL_SUFFIX(cstbsv,CSTBSV)
#define F77_stpsv 		F77_GLOBAL_SUFFIX(cstpsv,CSTPSV)
#define F77_dgemv 		F77_GLOBAL_SUFFIX(cdgemv,CDGEMV)
#define F77_dgbmv 		F77_GLOBAL_SUFFIX(cdgbmv,CDGBMV)
#define F77_dtrmv 		F77_GLOBAL_SUFFIX(cdtrmv,CDTRMV)
#define F77_dtbmv 		F77_GLOBAL_SUFFIX(cdtbmv,CDTBMV)
#define F77_dtpmv 		F77_GLOBAL_SUFFIX(cdtpmv,CDTPMV)
#define F77_dtrsv 		F77_GLOBAL_SUFFIX(cdtrsv,CDTRSV)
#define F77_dtbsv 		F77_GLOBAL_SUFFIX(cdtbsv,CDTBSV)
#define F77_dtpsv 		F77_GLOBAL_SUFFIX(cdtpsv,CDTPSV)
#define F77_cgemv 		F77_GLOBAL_SUFFIX(ccgemv,CCGEMV)
#define F77_cgbmv 		F77_GLOBAL_SUFFIX(ccgbmv,CCGBMV)
#define F77_ctrmv 		F77_GLOBAL_SUFFIX(cctrmv,CCTRMV)
#define F77_ctbmv 		F77_GLOBAL_SUFFIX(cctbmv,CCTBMV)
#define F77_ctpmv 		F77_GLOBAL_SUFFIX(cctpmv,CCTPMV)
#define F77_ctrsv 		F77_GLOBAL_SUFFIX(cctrsv,CCTRSV)
#define F77_ctbsv 		F77_GLOBAL_SUFFIX(cctbsv,CCTBSV)
#define F77_ctpsv 		F77_GLOBAL_SUFFIX(cctpsv,CCTPSV)
#define F77_zgemv 		F77_GLOBAL_SUFFIX(czgemv,CZGEMV)
#define F77_zgbmv 		F77_GLOBAL_SUFFIX(czgbmv,CZGBMV)
#define F77_ztrmv 		F77_GLOBAL_SUFFIX(cztrmv,CZTRMV)
#define F77_ztbmv 		F77_GLOBAL_SUFFIX(cztbmv,CZTBMV)
#define F77_ztpmv 		F77_GLOBAL_SUFFIX(cztpmv,CZTPMV)
#define F77_ztrsv 		F77_GLOBAL_SUFFIX(cztrsv,CZTRSV)
#define F77_ztbsv 		F77_GLOBAL_SUFFIX(cztbsv,CZTBSV)
#define F77_ztpsv 		F77_GLOBAL_SUFFIX(cztpsv,CZTPSV)
/*
 * Level 3 BLAS
 */
#define F77_s3chke 		F77_GLOBAL_SUFFIX(cs3chke,CS3CHKE)
#define F77_d3chke 		F77_GLOBAL_SUFFIX(cd3chke,CD3CHKE)
#define F77_c3chke 		F77_GLOBAL_SUFFIX(cc3chke,CC3CHKE)
#define F77_z3chke 		F77_GLOBAL_SUFFIX(cz3chke,CZ3CHKE)
#define F77_chemm 		F77_GLOBAL_SUFFIX(cchemm,CCHEMM)
#define F77_cherk 		F77_GLOBAL_SUFFIX(ccherk,CCHERK)
#define F77_cher2k 		F77_GLOBAL_SUFFIX(ccher2k,CCHER2K)
#define F77_zhemm 		F77_GLOBAL_SUFFIX(czhemm,CZHEMM)
#define F77_zherk 		F77_GLOBAL_SUFFIX(czherk,CZHERK)
#define F77_zher2k 		F77_GLOBAL_SUFFIX(czher2k,CZHER2K)
#define F77_sgemm 		F77_GLOBAL_SUFFIX(csgemm,CSGEMM)
#define F77_sgemmtr 		F77_GLOBAL_SUFFIX(csgemmtr,CSGEMMTR)
#define F77_ssymm 		F77_GLOBAL_SUFFIX(cssymm,CSSYMM)
#define F77_sskewsymm F77_GLOBAL_SUFFIX(csskewsymm,CSSKEWSYMM)
#define F77_ssyrk 		F77_GLOBAL_SUFFIX(cssyrk,CSSYRK)
#define F77_ssyr2k 		F77_GLOBAL_SUFFIX(cssyr2k,CSSYR2K)
#define F77_sskewsyr2k 	F77_GLOBAL_SUFFIX(csskewsyr2k,CSSKEWSYR2K)
#define F77_strmm 		F77_GLOBAL_SUFFIX(cstrmm,CSTRMM)
#define F77_strsm 		F77_GLOBAL_SUFFIX(cstrsm,CSTRSM)
#define F77_dgemm 		F77_GLOBAL_SUFFIX(cdgemm,CDGEMM)
#define F77_dgemmtr 		F77_GLOBAL_SUFFIX(cdgemmtr,CDGEMMTR)
#define F77_dsymm 		F77_GLOBAL_SUFFIX(cdsymm,CDSYMM)
#define F77_dskewsymm F77_GLOBAL_SUFFIX(cdskewsymm,CDSKEWSYMM)
#define F77_dsyrk 		F77_GLOBAL_SUFFIX(cdsyrk,CDSYRK)
#define F77_dsyr2k 		F77_GLOBAL_SUFFIX(cdsyr2k,CDSYR2K)
#define F77_dskewsyr2k 	F77_GLOBAL_SUFFIX(cdskewsyr2k,CDSKEWSYR2K)
#define F77_dtrmm 		F77_GLOBAL_SUFFIX(cdtrmm,CDTRMM)
#define F77_dtrsm 		F77_GLOBAL_SUFFIX(cdtrsm,CDTRSM)
#define F77_cgemm 		F77_GLOBAL_SUFFIX(ccgemm,CCGEMM)
#define F77_cgemmtr 		F77_GLOBAL_SUFFIX(ccgemmtr,CCGEMMTR)
#define F77_csymm 		F77_GLOBAL_SUFFIX(ccsymm,CCSYMM)
#define F77_csyrk 		F77_GLOBAL_SUFFIX(ccsyrk,CCSYRK)
#define F77_csyr2k 		F77_GLOBAL_SUFFIX(ccsyr2k,CCSYR2K)
#define F77_ctrmm 		F77_GLOBAL_SUFFIX(cctrmm,CCTRMM)
#define F77_ctrsm 		F77_GLOBAL_SUFFIX(cctrsm,CCTRSM)
#define F77_zgemm 		F77_GLOBAL_SUFFIX(czgemm,CZGEMM)
#define F77_zgemmtr 		F77_GLOBAL_SUFFIX(czgemmtr,CZGEMMTR)
#define F77_zsymm 		F77_GLOBAL_SUFFIX(czsymm,CZSYMM)
#define F77_zsyrk 		F77_GLOBAL_SUFFIX(czsyrk,CZSYRK)
#define F77_zsyr2k 		F77_GLOBAL_SUFFIX(czsyr2k,CZSYR2K)
#define F77_ztrmm 		F77_GLOBAL_SUFFIX(cztrmm,CZTRMM)
#define F77_ztrsm 		F77_GLOBAL_SUFFIX(cztrsm, CZTRSM)

void get_transpose_type(char *type, CBLAS_TRANSPOSE *trans);
void get_uplo_type(char *type, CBLAS_UPLO *uplo);
void get_diag_type(char *type, CBLAS_DIAG *diag);
void get_side_type(char *type, CBLAS_SIDE *side);

#endif /* CBLAS_TEST_H */
