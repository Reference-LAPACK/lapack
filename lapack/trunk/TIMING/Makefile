#######################################################################
#  This makefile runs the timing programs for the linear equation routines
#  and the eigenvalue routines in LAPACK.  The timing output files
#  are grouped as follows:
#
#       SLINTIM,SEIGTIM  -- Single precision real time routines
#       CLINTIM,CEIGTIM  -- Single precision complex time routines
#       DLINTIM,DEIGTIM  -- Double precision real time routines
#       ZLINTIM,ZEIGTIM  -- Double precision complex time routines
#
#  Timing programs can be executed for all or some of the four different
#  precisions.  Enter 'make' followed by one or more of the data types
#  desired.
#
#  By default, 'make' alone executes the timings with SMALL data sets.
#  (the data files with _small in the name)
#  
#  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments runs all eight test programs with the SMALL
#  timing data files.
#
#  If you wish to use the LARGE data sets for your timingsi ( the data
#  files with _large in their names), you
#  need to type
#       make large
#  for all tests, which would be equivalent to typing
#       make single_large double_large complex_large complex16_large
#  which could also be executed separately, for example, as
#       make single_large
#
#  The executable files are called:
#       xlintims, xlintimd, xlintimc, and xlintimz for LIN
#       xeigtims, xeigtimd, xeigtimc, and xeigtimz for EIG
#  and exist in the current directory level.
#
#  To remove the output files after the tests have been run, enter
#       make clean
#
#  To re-run specific tests after a make, enter (for example):
#       'rm stime.out; make'  or:
#       'make stime.out' or:
#       'touch stime.in; make' (to re-run the timings.)
#
#       'rm *time.out; make'  (to re-run all the timings.)
#
#######################################################################

include ../make.inc

all:
	$(MAKE) single
	$(MAKE) complex
	$(MAKE) double
	$(MAKE) complex16
large:
	$(MAKE) single_large
	$(MAKE) complex_large
	$(MAKE) double_large
	$(MAKE) complex16_large
small:
	$(MAKE) single
	$(MAKE) complex
	$(MAKE) double
	$(MAKE) complex16

SEIGTIM= sgeptim_small.out \
         sneptim_small.out \
         sseptim_small.out \
         ssvdtim_small.out

SEIGTM2= sgeptim_large.out \
         sneptim_large.out \
         sseptim_large.out \
         ssvdtim_large.out

CEIGTIM= cgeptim_small.out \
         cneptim_small.out \
         cseptim_small.out \
         csvdtim_small.out

CEIGTM2= cgeptim_large.out \
         cneptim_large.out \
         cseptim_large.out \
         csvdtim_large.out

DEIGTIM= dgeptim_small.out \
         dneptim_small.out \
         dseptim_small.out \
         dsvdtim_small.out

DEIGTM2= dgeptim_large.out \
         dneptim_large.out \
         dseptim_large.out \
         dsvdtim_large.out

ZEIGTIM= zgeptim_small.out \
         zneptim_small.out \
         zseptim_small.out \
         zsvdtim_small.out

ZEIGTM2= zgeptim_large.out \
         zneptim_large.out \
         zseptim_large.out \
         zsvdtim_large.out

SLINTIM= stime_small.out \
	 sband_small.out \
	 stime2_small.out

SLINTM2= stime_large.out \
	 sband_large.out \
	 stime2_large.out

CLINTIM= ctime_small.out \
	 cband_small.out \
	 ctime2_small.out

CLINTM2= ctime_large.out \
	 cband_large.out \
	 ctime2_large.out

DLINTIM= dtime_small.out \
	 dband_small.out \
	 dtime2_small.out

DLINTM2= dtime_large.out \
	 dband_large.out \
	 dtime2_large.out

ZLINTIM= ztime_small.out \
	 zband_small.out \
	 ztime2_small.out

ZLINTM2= ztime_large.out \
	 zband_large.out \
	 ztime2_large.out

single:	xlintims xeigtims
	for outfile in $(SLINTIM) $(SEIGTIM) ; do $(MAKE) $$outfile; done
complex:	xlintimc xeigtimc
	for outfile in $(CLINTIM) $(CEIGTIM) ; do $(MAKE) $$outfile; done
double:	xlintimd xeigtimd
	for outfile in $(DLINTIM) $(DEIGTIM) ; do $(MAKE) $$outfile; done
complex16:	xlintimz xeigtimz
	for outfile in $(ZLINTIM) $(ZEIGTIM) ; do $(MAKE) $$outfile; done

single_large:	xlintims xeigtims
	for outfile in $(SLINTM2) $(SEIGTM2) ; do $(MAKE) $$outfile; done
complex_large:	xlintimc xeigtimc
	for outfile in $(CLINTM2) $(CEIGTM2) ; do $(MAKE) $$outfile; done
double_large:	xlintimd xeigtimd
	for outfile in $(DLINTM2) $(DEIGTM2) ; do $(MAKE) $$outfile; done
complex16_large:	xlintimz xeigtimz
	for outfile in $(ZLINTM2) $(ZEIGTM2) ; do $(MAKE) $$outfile; done

#
# ======== SINGLE LIN TIMINGS ===========================

stime_small.out: stime_small.in xlintims
	@echo Timing square REAL LAPACK linear equation routines
	./xlintims < stime_small.in > $@ 2>&1

stime_large.out: stime_large.in xlintims
	@echo Timing square REAL LAPACK linear equation routines
	./xlintims < stime_large.in > $@ 2>&1

sband_small.out: sband_small.in xlintims
	@echo Timing banded REAL LAPACK linear equation routines
	./xlintims < sband_small.in > $@ 2>&1

sband_large.out: sband_large.in xlintims
	@echo Timing banded REAL LAPACK linear equation routines
	./xlintims < sband_large.in > $@ 2>&1

stime2_small.out: stime2_small.in xlintims
	@echo Timing rectangular REAL LAPACK linear equation routines
	./xlintims < stime2_small.in > $@ 2>&1

stime2_large.out: stime2_large.in xlintims
	@echo Timing rectangular REAL LAPACK linear equation routines
	./xlintims < stime2_large.in > $@ 2>&1
#
# ======== COMPLEX LIN TIMINGS ==========================

ctime_small.out: ctime_small.in xlintimc
	@echo Timing square COMPLEX LAPACK linear equation routines
	./xlintimc < ctime_small.in > $@ 2>&1

ctime_large.out: ctime_large.in xlintimc
	@echo Timing square COMPLEX LAPACK linear equation routines
	./xlintimc < ctime_large.in > $@ 2>&1

cband_small.out: cband_small.in xlintimc
	@echo Timing banded COMPLEX LAPACK linear equation routines
	./xlintimc < cband_small.in > $@ 2>&1

cband_large.out: cband_large.in xlintimc
	@echo Timing banded COMPLEX LAPACK linear equation routines
	./xlintimc < cband_large.in > $@ 2>&1

ctime2_small.out: ctime2_small.in xlintimc
	@echo Timing rectangular COMPLEX LAPACK linear equation routines
	./xlintimc < ctime2_small.in > $@ 2>&1

ctime2_large.out: CTIME2_large.in xlintimc
	@echo Timing rectangular COMPLEX LAPACK linear equation routines
	./xlintimc < ctime2_large.in > $@ 2>&1
#
# ======== DOUBLE LIN TIMINGS ===========================

dtime_small.out: dtime_small.in xlintimd
	@echo Timing square DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < dtime_small.in > $@ 2>&1

dtime_large.out: dtime_large.in xlintimd
	@echo Timing square DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < dtime_large.in > $@ 2>&1

dband_small.out: dband_small.in xlintimd
	@echo Timing banded DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < dband_small.in > $@ 2>&1

dband_large.out: dband_large.in xlintimd
	@echo Timing banded DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < dband_large.in > $@ 2>&1

dtime2_small.out: dtime2_small.in xlintimd
	@echo Timing rectangular DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < dtime2_small.in > $@ 2>&1

dtime2_large.out: dtime2_large.in xlintimd
	@echo Timing rectangular DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < dtime2_large.in > $@ 2>&1
#
# ======== COMPLEX16 LIN TIMINGS ========================

ztime_small.out: ztime_small.in xlintimz
	@echo Timing square COMPLEX16 LAPACK linear equation routines
	./xlintimz < ztime_small.in > $@ 2>&1

ztime_large.out: ztime_large.in xlintimz
	@echo Timing square COMPLEX16 LAPACK linear equation routines
	./xlintimz < ztime_large.in > $@ 2>&1

zband_small.out: zband_small.in xlintimz
	@echo Timing banded COMPLEX16 LAPACK linear equation routines
	./xlintimz < zband_small.in > $@ 2>&1

zband_large.out: zband_large.in xlintimz
	@echo Timing banded COMPLEX16 LAPACK linear equation routines
	./xlintimz < zband_large.in > $@ 2>&1

ztime2_small.out: ztime2_small.in xlintimz
	@echo Timing rectangular COMPLEX16 LAPACK linear equation routines
	./xlintimz < ztime2_small.in > $@ 2>&1

ztime2_large.out: ztime2_large.in xlintimz
	@echo Timing rectangular COMPLEX16 LAPACK linear equation routines
	./xlintimz < ztime2_large.in > $@ 2>&1
#
#
# ======== SINGLE EIG TIMINGS ===========================
#
sgeptim_small.out: sgeptim_small.in xeigtims
	@echo GEP: Timing REAL Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtims < sgeptim_small.in > $@ 2>&1

sgeptim_large.out: sgeptim_large.in xeigtims
	@echo GEP: Timing REAL Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtims < sgeptim_large.in > $@ 2>&1

sneptim_small.out: sneptim_small.in xeigtims
	@echo NEP: Timing REAL Nonsymmetric Eigenvalue Problem routines
	./xeigtims < sneptim_small.in > $@ 2>&1

sneptim_large.out: sneptim_large.in xeigtims
	@echo NEP: Timing REAL Nonsymmetric Eigenvalue Problem routines
	./xeigtims < sneptim_large.in > $@ 2>&1

sseptim_small.out: sseptim_small.in xeigtims
	@echo SEP: Timing REAL Symmetric Eigenvalue Problem routines
	./xeigtims < sseptim_small.in > $@ 2>&1

sseptim_large.out: sseptim_large.in xeigtims
	@echo SEP: Timing REAL Symmetric Eigenvalue Problem routines
	./xeigtims < sseptim_large.in > $@ 2>&1

ssvdtim_small.out: ssvdtim_small.in xeigtims
	@echo SVD: Timing REAL Singular Value Decomposition routines
	./xeigtims < ssvdtim_small.in > $@ 2>&1

ssvdtim_large.out: ssvdtim_large.in xeigtims
	@echo SVD: Timing REAL Singular Value Decomposition routines
	./xeigtims < ssvdtim_large.in > $@ 2>&1
#
# ======== COMPLEX EIG TIMINGS ===========================
#
cgeptim_small.out: cgeptim_small.in xeigtimc
	@echo GEP: Timing COMPLEX Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimc < cgeptim_small.in > $@ 2>&1

cgeptim_large.out: cgeptim_large.in xeigtimc
	@echo GEP: Timing COMPLEX Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimc < cgeptim_large_large.in > $@ 2>&1

cneptim_small.out: cneptim_small.in xeigtimc
	@echo NEP: Timing COMPLEX Nonsymmetric Eigenvalue Problem routines
	./xeigtimc < cneptim_small.in > $@ 2>&1

cneptim_large.out: cneptim_large.in xeigtimc
	@echo NEP: Timing COMPLEX Nonsymmetric Eigenvalue Problem routines
	./xeigtimc < cneptim_large.in > $@ 2>&1

cseptim_small.out: cseptim_small.in xeigtimc
	@echo SEP: Timing COMPLEX Symmetric Eigenvalue Problem routines
	./xeigtimc < cseptim_small.in > $@ 2>&1

cseptim_large.out: cseptim_large.in xeigtimc
	@echo SEP: Timing COMPLEX Symmetric Eigenvalue Problem routines
	./xeigtimc < cseptim_large.in > $@ 2>&1

csvdtim_small.out: csvdtim_small.in xeigtimc
	@echo SVD: Timing COMPLEX Singular Value Decomposition routines
	./xeigtimc < csvdtim_small.in > $@ 2>&1

csvdtim_large.out: csvdtim_large.in xeigtimc
	@echo SVD: Timing COMPLEX Singular Value Decomposition routines
	./xeigtimc < csvdtim_large.in > $@ 2>&1
#
# ======== DOUBLE EIG TIMINGS ===========================
#
dgeptim_small.out: dgeptim_small.in xeigtimd
	@echo GEP: Timing DOUBLE PRECISION Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimd < dgeptim_small.in > $@ 2>&1

dgeptim_large.out: dgeptim_large.in xeigtimd
	@echo GEP: Timing DOUBLE PRECISION Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimd < dgeptim_large.in > $@ 2>&1

dneptim_small.out: dneptim_small.in xeigtimd
	@echo NEP: Timing DOUBLE PRECISION Nonsymmetric Eigenvalue Problem routines
	./xeigtimd < dneptim_small.in > $@ 2>&1

dneptim_large.out: dneptim_large.in xeigtimd
	@echo NEP: Timing DOUBLE PRECISION Nonsymmetric Eigenvalue Problem routines
	./xeigtimd < dneptim_large.in > $@ 2>&1

dseptim_small.out: dseptim_small.in xeigtimd
	@echo SEP: Timing DOUBLE PRECISION Symmetric Eigenvalue Problem routines
	./xeigtimd < dseptim_small.in > $@ 2>&1

dseptim_large.out: dseptim_large.in xeigtimd
	@echo SEP: Timing DOUBLE PRECISION Symmetric Eigenvalue Problem routines
	./xeigtimd < dseptim_large.in > $@ 2>&1

dsvdtim_small.out: dsvdtim_small.in xeigtimd
	@echo SVD: Timing DOUBLE PRECISION Singular Value Decomposition routines
	./xeigtimd < dsvdtim_small.in > $@ 2>&1

dsvdtim_large.out: dsvdtim_large.in xeigtimd
	@echo SVD: Timing DOUBLE PRECISION Singular Value Decomposition routines
	./xeigtimd < dsvdtim_large.in > $@ 2>&1
#
# ======== COMPLEX16 EIG TIMINGS ===========================
#
zgeptim_small.out: zgeptim_small.in xeigtimz
	@echo GEP: Timing COMPLEX16 Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimz < zgeptim_small.in > $@ 2>&1

zgeptim_large.out: zgeptim_large.in xeigtimz
	@echo GEP: Timing COMPLEX16 Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimz < zgeptim_large.in > $@ 2>&1

zneptim_small.out: zneptim_small.in xeigtimz
	@echo NEP: Timing COMPLEX16 Nonsymmetric Eigenvalue Problem routines
	./xeigtimz < zneptim_small.in > $@ 2>&1

zneptim_large.out: zneptim_large.in xeigtimz
	@echo NEP: Timing COMPLEX16 Nonsymmetric Eigenvalue Problem routines
	./xeigtimz < zneptim_large.in > $@ 2>&1

zseptim_small.out: zseptim_small.in xeigtimz
	@echo SEP: Timing COMPLEX16 Symmetric Eigenvalue Problem routines
	./xeigtimz < zseptim_small.in > $@ 2>&1

zseptim_large.out: zseptim_large.in xeigtimz
	@echo SEP: Timing COMPLEX16 Symmetric Eigenvalue Problem routines
	./xeigtimz < zseptim_large.in > $@ 2>&1

zsvdtim_small.out: zsvdtim_small.in xeigtimz
	@echo SVD: Timing COMPLEX16 Singular Value Decomposition routines
	./xeigtimz < zsvdtim_small.in > $@ 2>&1

zsvdtim_large.out: zsvdtim_large.in xeigtimz
	@echo SVD: Timing COMPLEX16 Singular Value Decomposition routines
	./xeigtimz < zsvdtim_large.in > $@ 2>&1
# ==============================================================================

xlintims:
	cd LIN ; $(MAKE) single

xlintimc:
	cd LIN ; $(MAKE) complex

xlintimd:
	cd LIN ; $(MAKE) double

xlintimz:
	cd LIN ; $(MAKE) complex16


xeigtims:
	cd EIG ; $(MAKE) single

xeigtimc:
	cd EIG ; $(MAKE) complex

xeigtimd:
	cd EIG ; $(MAKE) double

xeigtimz:
	cd EIG ; $(MAKE) complex16

clean:
	rm -f *.out core

cleanup:
	rm -f x* *.out core
