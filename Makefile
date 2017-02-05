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
	$(MAKE) -C BLAS/SRC

cblaslib:
	$(MAKE) -C CBLAS

lapacklib:
	$(MAKE) -C SRC

lapackelib:
	$(MAKE) -C LAPACKE

tmglib:
	$(MAKE) -C TESTING/MATGEN

variants:
	$(MAKE) -C SRC/VARIANTS

lapack_install:
	$(MAKE) -C INSTALL run

blas_testing: blaslib
	$(MAKE) -C BLAS/TESTING run

cblas_testing: cblaslib blaslib
	$(MAKE) -C CBLAS cblas_testing

lapack_testing: tmglib lapacklib blaslib
	rm -f TESTING/LIN/xlintst*
	$(MAKE) -C TESTING
	./lapack_testing.py

variants_testing: tmglib variants lapacklib blaslib
	rm -f TESTING/LIN/xlintst*
	$(MAKE) -C TESTING/LIN VARLIB='SRC/VARIANTS/cholrl.a'
	$(MAKE) -C TESTING stest.out && mv TESTING/stest.out TESTING/stest_cholrl.out
	$(MAKE) -C TESTING dtest.out && mv TESTING/dtest.out TESTING/dtest_cholrl.out
	$(MAKE) -C TESTING ctest.out && mv TESTING/ctest.out TESTING/ctest_cholrl.out
	$(MAKE) -C TESTING ztest.out && mv TESTING/ztest.out TESTING/ztest_cholrl.out
	rm -f TESTING/LIN/xlintst*
	$(MAKE) -C TESTING/LIN VARLIB='SRC/VARIANTS/choltop.a'
	$(MAKE) -C TESTING stest.out && mv TESTING/stest.out TESTING/stest_choltop.out
	$(MAKE) -C TESTING dtest.out && mv TESTING/dtest.out TESTING/dtest_choltop.out
	$(MAKE) -C TESTING ctest.out && mv TESTING/ctest.out TESTING/ctest_choltop.out
	$(MAKE) -C TESTING ztest.out && mv TESTING/ztest.out TESTING/ztest_choltop.out
	rm -f TESTING/LIN/xlintst*
	$(MAKE) -C TESTING/LIN VARLIB='SRC/VARIANTS/lucr.a'
	$(MAKE) -C TESTING stest.out && mv TESTING/stest.out TESTING/stest_lucr.out
	$(MAKE) -C TESTING dtest.out && mv TESTING/dtest.out TESTING/dtest_lucr.out
	$(MAKE) -C TESTING ctest.out && mv TESTING/ctest.out TESTING/ctest_lucr.out
	$(MAKE) -C TESTING ztest.out && mv TESTING/ztest.out TESTING/ztest_lucr.out
	rm -f TESTING/LIN/xlintst*
	$(MAKE) -C TESTING/LIN VARLIB='SRC/VARIANTS/lull.a'
	$(MAKE) -C TESTING stest.out && mv TESTING/stest.out TESTING/stest_lull.out
	$(MAKE) -C TESTING dtest.out && mv TESTING/dtest.out TESTING/dtest_lull.out
	$(MAKE) -C TESTING ctest.out && mv TESTING/ctest.out TESTING/ctest_lull.out
	$(MAKE) -C TESTING ztest.out && mv TESTING/ztest.out TESTING/ztest_lull.out
	rm -f TESTING/LIN/xlintst*
	$(MAKE) -C TESTING/LIN VARLIB='SRC/VARIANTS/lurec.a'
	$(MAKE) -C TESTING stest.out && mv TESTING/stest.out TESTING/stest_lurec.out
	$(MAKE) -C TESTING dtest.out && mv TESTING/dtest.out TESTING/dtest_lurec.out
	$(MAKE) -C TESTING ctest.out && mv TESTING/ctest.out TESTING/ctest_lurec.out
	$(MAKE) -C TESTING ztest.out && mv TESTING/ztest.out TESTING/ztest_lurec.out
	rm -f TESTING/LIN/xlintst*
	$(MAKE) -C TESTING/LIN VARLIB='SRC/VARIANTS/qrll.a'
	$(MAKE) -C TESTING stest.out && mv TESTING/stest.out TESTING/stest_qrll.out
	$(MAKE) -C TESTING dtest.out && mv TESTING/dtest.out TESTING/dtest_qrll.out
	$(MAKE) -C TESTING ctest.out && mv TESTING/ctest.out TESTING/ctest_qrll.out
	$(MAKE) -C TESTING ztest.out && mv TESTING/ztest.out TESTING/ztest_qrll.out

cblas_example: cblaslib blaslib
	$(MAKE) -C CBLAS cblas_example

lapacke_example: lapackelib lapacklib blaslib
	$(MAKE) -C LAPACKE lapacke_example

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
	$(MAKE) -C INSTALL clean
	$(MAKE) -C BLAS/SRC clean
	$(MAKE) -C CBLAS clean
	$(MAKE) -C SRC clean
	$(MAKE) -C SRC/VARIANTS clean
	$(MAKE) -C TESTING/MATGEN clean
	$(MAKE) -C LAPACKE clean


cleanblas_testing:
	$(MAKE) -C BLAS/TESTING clean

cleancblas_testing:
	$(MAKE) -C CBLAS/testing clean

cleantesting:
	$(MAKE) -C TESTING/LIN clean
	$(MAKE) -C TESTING/EIG clean

cleanall: cleanlib cleanblas_testing cleancblas_testing cleantesting
	$(MAKE) -C INSTALL cleanall
	rm -f *.a TESTING/*.out
