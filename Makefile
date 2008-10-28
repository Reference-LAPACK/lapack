#
#  Top Level Makefile for LAPACK
#  Version 3.1.1
#  February 2007
#

include make.inc

all: lapack_install lib lapack_testing blas_testing

lib: lapacklib variants tmglib
#lib: blaslib lapacklib variants tmglib

clean: cleanlib cleantesting cleanblas_testing cleantiming

lapack_install:
	( cd INSTALL; $(MAKE); ./testlsame; ./testslamch; \
	  ./testdlamch; ./testsecond; ./testdsecnd; ./testversion )

blaslib:
	( cd BLAS/SRC; $(MAKE) )

xblaslib:
	( cd XBLAS;  $(MAKE) lib )

lapacklib:	lapack_install
	( cd SRC; $(MAKE) )
	
variants:
	( cd SRC/VARIANTS ; $(MAKE))

tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

lapack_testing:	lib
	( cd TESTING ; $(MAKE) )
	
variants_testing: lib
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/cholrl.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_cholrl.out ; mv dtest.out dtest_cholrl.out ; mv ctest.out ctest_cholrl.out ; mv ztest.out ztest_cholrl.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/choltop.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_choltop.out ; mv dtest.out dtest_choltop.out ; mv ctest.out ctest_choltop.out ; mv ztest.out ztest_choltop.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/lucr.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_lucr.out ; mv dtest.out dtest_lucr.out ; mv ctest.out ctest_lucr.out ; mv ztest.out ztest_lucr.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/lull.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_lull.out ; mv dtest.out dtest_lull.out ; mv ctest.out ctest_lull.out ; mv ztest.out ztest_lull.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/lurec.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_lurec.out ; mv dtest.out dtest_lurec.out ; mv ctest.out ctest_lurec.out ; mv ztest.out ztest_lurec.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/qrll.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_qrll.out ; mv dtest.out dtest_qrll.out ; mv ctest.out ctest_qrll.out ; mv ztest.out ztest_qrll.out )

blas_testing:
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 )
	( cd BLAS; ./xblat1s > sblat1.out    ; \
	           ./xblat1d > dblat1.out    ; \
	           ./xblat1c > cblat1.out    ; \
	           ./xblat1z > zblat1.out    ) 
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 )
	( cd BLAS; ./xblat2s < sblat2.in     ; \
	           ./xblat2d < dblat2.in     ; \
	           ./xblat2c < cblat2.in     ; \
	           ./xblat2z < zblat2.in     )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 )
	( cd BLAS; ./xblat3s < sblat3.in     ; \
	           ./xblat3d < dblat3.in     ; \
	           ./xblat3c < cblat3.in     ; \
	           ./xblat3z < zblat3.in     ) 

xblas_testing:
	( cd XBLAS;  $(MAKE) tests )

lapack_timing:	lib lapack_testing blas_testing
	( cd TIMING; $(MAKE) )

blas_timing:	lapack_timing
	( cd TIMING/LIN; $(MAKE) )
	( cd TIMING; ./xlintims < sblasa_small.in > sblasa_small.out ; \
	             ./xlintims < sblasb_small.in > sblasb_small.out ; \
	             ./xlintims < sblasc_small.in > sblasc_small.out )
	( cd TIMING; ./xlintimd < dblasa_small.in > dblasa_small.out ; \
	             ./xlintimd < dblasb_small.in > dblasb_small.out ; \
	             ./xlintimd < dblasc_small.in > dblasc_small.out )
	( cd TIMING; ./xlintimc < cblasa_small.in > cblasa_small.out ; \
	             ./xlintimc < cblasb_small.in > cblasb_small.out ; \
	             ./xlintimc < cblasc_small.in > cblasc_small.out )
	( cd TIMING; ./xlintimz < zblasa_small.in > zblasa_small.out ; \
	             ./xlintimz < zblasb_small.in > zblasb_small.out ; \
	             ./xlintimz < zblasc_small.in > zblasc_small.out )

cleanlib:
	( cd INSTALL; $(MAKE) clean )
	( cd BLAS/SRC; $(MAKE) clean )
	( cd SRC; $(MAKE) clean )
	( cd SRC/VARIANTS; $(MAKE) clean )
	( cd TESTING/MATGEN; $(MAKE) clean )
	( cd XBLAS; $(MAKE) clean )

cleanblas_testing:	
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 clean )
	( cd BLAS; rm -f xblat* )

cleanxblas:
	( cd XBLAS; $(MAKE) clean ) 

cleantesting:
	( cd TESTING/LIN; $(MAKE) clean )
	( cd TESTING/EIG; $(MAKE) clean )
	( cd TESTING; rm -f xlin* xeig* )

cleantiming:
	( cd TIMING/LIN; $(MAKE) clean )
	( cd TIMING/LIN/LINSRC; $(MAKE) clean )
	( cd TIMING/EIG; $(MAKE) clean )
	( cd TIMING/EIG/EIGSRC; $(MAKE) clean )
	( cd TIMING; rm -f xlin* xeig* )

cleanall: cleanlib cleanblas_testing cleantesting cleantiming
	(cd XBLAS; $(MAKE) dist-clean )
	rm -f *.a TESTING/*.out TIMING/*.out INSTALL/test* \
                 BLAS/*.out TIMING/LIN/*.a TIMING/EIG/*.a

