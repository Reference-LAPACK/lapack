#
#  Top Level Makefile for LAPACK
#  Version 3.4.1
#  April 2012
#

include make.inc

all: lapack_install lib blas_testing lapack_testing

lib: lapacklib tmglib
#lib: blaslib variants lapacklib tmglib

clean: cleanlib cleantesting cleanblas_testing cleancblas_testing

blaslib:
	( cd BLAS/SRC; $(MAKE) )

cblaslib:
	( cd CBLAS; $(MAKE) )

lapacklib:
	( cd SRC; $(MAKE) )

lapackelib:
	( cd LAPACKE; $(MAKE) )

tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

variants:
	( cd SRC/VARIANTS; $(MAKE) )

lapack_install:
	( cd INSTALL; $(MAKE) run )

blas_testing: blaslib
	( cd BLAS/TESTING; $(MAKE) run )

cblas_testing: cblaslib blaslib
	( cd CBLAS; $(MAKE) cblas_testing )

lapack_testing: tmglib lapacklib blaslib
	( cd TESTING/LIN && rm -f xlintst* )
	( cd TESTING && $(MAKE) )
	./lapack_testing.py

variants_testing: tmglib variants lapacklib blaslib
	( cd TESTING/LIN && rm -f xlintst* && $(MAKE) VARLIB='SRC/VARIANTS/cholrl.a' )
	( cd TESTING && $(MAKE) stest.out && mv stest.out stest_cholrl.out )
	( cd TESTING && $(MAKE) dtest.out && mv dtest.out dtest_cholrl.out )
	( cd TESTING && $(MAKE) ctest.out && mv ctest.out ctest_cholrl.out )
	( cd TESTING && $(MAKE) ztest.out && mv ztest.out ztest_cholrl.out )
	( cd TESTING/LIN; rm -f xlintst*; $(MAKE) VARLIB='SRC/VARIANTS/choltop.a' )
	( cd TESTING && $(MAKE) stest.out && mv stest.out stest_choltop.out )
	( cd TESTING && $(MAKE) dtest.out && mv dtest.out dtest_choltop.out )
	( cd TESTING && $(MAKE) ctest.out && mv ctest.out ctest_choltop.out )
	( cd TESTING && $(MAKE) ztest.out && mv ztest.out ztest_choltop.out )
	( cd TESTING/LIN; rm -f xlintst*; $(MAKE) VARLIB='SRC/VARIANTS/lucr.a' )
	( cd TESTING && $(MAKE) stest.out && mv stest.out stest_lucr.out )
	( cd TESTING && $(MAKE) dtest.out && mv dtest.out dtest_lucr.out )
	( cd TESTING && $(MAKE) ctest.out && mv ctest.out ctest_lucr.out )
	( cd TESTING && $(MAKE) ztest.out && mv ztest.out ztest_lucr.out )
	( cd TESTING/LIN; rm -f xlintst*; $(MAKE) VARLIB='SRC/VARIANTS/lull.a' )
	( cd TESTING && $(MAKE) stest.out && mv stest.out stest_lull.out )
	( cd TESTING && $(MAKE) dtest.out && mv dtest.out dtest_lull.out )
	( cd TESTING && $(MAKE) ctest.out && mv ctest.out ctest_lull.out )
	( cd TESTING && $(MAKE) ztest.out && mv ztest.out ztest_lull.out )
	( cd TESTING/LIN; rm -f xlintst*; $(MAKE) VARLIB='SRC/VARIANTS/lurec.a' )
	( cd TESTING && $(MAKE) stest.out && mv stest.out stest_lurec.out )
	( cd TESTING && $(MAKE) dtest.out && mv dtest.out dtest_lurec.out )
	( cd TESTING && $(MAKE) ctest.out && mv ctest.out ctest_lurec.out )
	( cd TESTING && $(MAKE) ztest.out && mv ztest.out ztest_lurec.out )
	( cd TESTING/LIN; rm -f xlintst*; $(MAKE) VARLIB='SRC/VARIANTS/qrll.a' )
	( cd TESTING && $(MAKE) stest.out && mv stest.out stest_qrll.out )
	( cd TESTING && $(MAKE) dtest.out && mv dtest.out dtest_qrll.out )
	( cd TESTING && $(MAKE) ctest.out && mv ctest.out ctest_qrll.out )
	( cd TESTING && $(MAKE) ztest.out && mv ztest.out ztest_qrll.out )

cblas_example: cblaslib blaslib
	( cd CBLAS; $(MAKE) cblas_example )

lapacke_example: lapackelib lapacklib blaslib
	( cd LAPACKE; $(MAKE) lapacke_example )

html:
	@echo "LAPACK HTML PAGES GENERATION with Doxygen"
	doxygen DOCS/Doxyfile
	@echo "=================="
	@echo "LAPACK HTML PAGES GENERATED in DOCS/explore-html"
	@echo "Usage: open DOCS/explore-html/index.html"
	@echo "Online version available at http://www.netlib.org/lapack/explore-html/"
	@echo "=================="

man:
	@echo "LAPACK MAN PAGES GENERATION with Doxygen"
	doxygen DOCS/Doxyfile_man
	@echo "=================="
	@echo "LAPACK MAN PAGES GENERATED in DOCS/MAN"
	@echo "Set your MANPATH env variable accordingly"
	@echo "Usage: man dgetrf.f"
	@echo "=================="

cleanlib:
	( cd INSTALL; $(MAKE) clean )
	( cd BLAS/SRC; $(MAKE) clean )
	( cd CBLAS; $(MAKE) clean )
	( cd SRC; $(MAKE) clean )
	( cd SRC/VARIANTS; $(MAKE) clean )
	( cd TESTING/MATGEN; $(MAKE) clean )
	( cd LAPACKE; $(MAKE) clean )


cleanblas_testing:
	( cd BLAS/TESTING; $(MAKE) clean )

cleancblas_testing:
	( cd CBLAS/testing; $(MAKE) clean )

cleantesting:
	( cd TESTING/LIN; $(MAKE) clean )
	( cd TESTING/EIG; $(MAKE) clean )

cleanall: cleanlib cleanblas_testing cleancblas_testing cleantesting
	( cd INSTALL; $(MAKE) cleanall )
	rm -f *.a TESTING/*.out
