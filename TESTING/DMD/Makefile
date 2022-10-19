#######################################################################
#  This makefile creates the test programs for the linear equation
#  routines in LAPACK.  The test files are grouped as follows:
#
#       ALINTST -- Auxiliary test routines
#       SLINTST -- Single precision real test routines
#       CLINTST -- Single precision complex test routines
#       SCLNTST -- Single precision real and complex routines in common
#       DLINTST -- Double precision real test routines
#       ZLINTST -- Double precision complex test routines
#       DZLNTST -- Double precision real and complex routines in common
#
#  Test programs can be generated for all or some of the four different
#  precisions.  Enter make followed by one or more of the data types
#  desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates all four test programs.
#  The executable files are called
#       xlintsts, xlintstd, xlintstc, and xlintstz
#
#  To remove the object files after the executable files have been
#  created, enter
#       make cleanobj
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################

TOPSRCDIR = ../..
include $(TOPSRCDIR)/make.inc

#ALINTST = \
#  aladhd.o alaerh.o alaesm.o alahd.o alareq.o \
#  alasum.o alasvm.o chkxer.o icopy.o ilaenv.o xlaenv.o xerbla.o

SLINTST = schkdmd.o

DZLNTST = dchkdmd.o

CLINTST = cchkdmd.o

ZLINTST = zchkdmd.o

.PHONY: all
all: single double complex complex16 proto-single proto-double proto-complex proto-complex16

.PHONY: single double complex complex16
single: xlintsts
double: xlintstd
complex: xlintstc
complex16: xlintstz

.PHONY: proto-single proto-double proto-complex proto-complex16
proto-single: xlintstrfs
proto-double: xlintstds xlintstrfd
proto-complex: xlintstrfc
proto-complex16: xlintstzc xlintstrfz

xlintsts: $(ALINTST) $(SLINTST) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(XBLASLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

xlintstc: $(ALINTST) $(CLINTST) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(XBLASLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

xlintstd: $(ALINTST) $(DLINTST) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(XBLASLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

xlintstz: $(ALINTST) $(ZLINTST) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(XBLASLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

xlintstds: $(DSLINTST) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

xlintstzc: $(ZCLINTST) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

xlintstrfs: $(SLINTSTRFP) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

xlintstrfd: $(DLINTSTRFP) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

xlintstrfc: $(CLINTSTRFP) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

xlintstrfz: $(ZLINTSTRFP) $(TMGLIB) $(VARLIB) $(LAPACKLIB) $(BLASLIB)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

$(ALINTST): $(FRC)
$(SCLNTST): $(FRC)
$(DZLNTST): $(FRC)
$(SLINTST): $(FRC)
$(CLINTST): $(FRC)
$(DLINTST): $(FRC)
$(ZLINTST): $(FRC)

FRC:
	@FRC=$(FRC)

.PHONY: clean cleanobj cleanexe
clean: cleanobj cleanexe
cleanobj:
	rm -f *.o
cleanexe:
	rm -f xlintst*

schkaa.o: schkaa.F
	$(FC) $(FFLAGS_DRV) -c -o $@ $<
dchkaa.o: dchkaa.F
	$(FC) $(FFLAGS_DRV) -c -o $@ $<
cchkaa.o: cchkaa.F
	$(FC) $(FFLAGS_DRV) -c -o $@ $<
zchkaa.o: zchkaa.F
	$(FC) $(FFLAGS_DRV) -c -o $@ $<

.NOTPARALLEL:
