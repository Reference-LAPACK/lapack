include ../../make.inc

#######################################################################
#  This makefile creates the timing programs for the linear equation
#  routines in LAPACK.  The timing files are grouped as follows:
#       ALINTIM -- Timing routines for all data types
#       SCINTIM -- Single precision real and complex timing routines
#       DZINTIM -- Double precision real and complex timing routines
#
#       SLINTIM -- Single precision real timing routines
#       CLINTIM -- Single precision complex timing routines
#       DLINTIM -- Double precision real timing routines
#       ZLINTIM -- Double precision complex timing routines
#
#  Timing programs can be generated for all or some of the four
#  data types.  Enter make followed by one or more of the data types
#  desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates all four timing programs.
#  The executable files are called
#       xlintims, xlintimd, xlintimc, and xlintimz
#  and are created in the next higher directory level.
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################
 
ALINTIM = atimin.o atimck.o icopy.o ilaenv.o xlaenv.o
 
SCINTIM = slaord.o sprtbl.o sprtb2.o sprtb3.o sprtb4.o sprtb5.o \
        sprtls.o smflop.o sopaux.o sopbl2.o sopbl3.o sopgb.o sopla.o
 
DZINTIM = dlaord.o dprtbl.o dprtb2.o dprtb3.o dprtb4.o dprtb5.o \
        dprtls.o dmflop.o dopaux.o dopbl2.o dopbl3.o dopgb.o dopla.o
 
SLINTIM = stimaa.o stimmg.o \
        stimmv.o stimmm.o stimb2.o stimb3.o stimge.o stimgb.o stimgt.o \
        stimls.o stimpo.o stimpp.o stimpb.o stimpt.o stimsy.o stimsp.o \
        stimtr.o stimtp.o stimtb.o stimhr.o stimtd.o stimbr.o \
        stimq3.o stimqr.o stimlq.o stimql.o stimrq.o stimqp.o \
        sqrt13.o sqrt15.o slinpk.o seispk.o
 
CLINTIM = ctimaa.o ctimmg.o \
        ctimmv.o ctimmm.o ctimb2.o ctimb3.o ctimge.o ctimgb.o ctimgt.o \
        ctimls.o ctimpo.o ctimpp.o ctimpb.o ctimpt.o \
        ctimsy.o ctimsp.o ctimhe.o ctimhp.o \
        ctimtr.o ctimtp.o ctimtb.o ctimhr.o ctimtd.o ctimbr.o \
        ctimq3.o ctimqr.o ctimlq.o ctimql.o ctimrq.o ctimqp.o \
        cqrt13.o cqrt15.o clinpk.o
 
DLINTIM = dtimaa.o dtimmg.o \
        dtimmv.o dtimmm.o dtimb2.o dtimb3.o dtimge.o dtimgb.o dtimgt.o \
        dtimls.o dtimpo.o dtimpp.o dtimpb.o dtimpt.o dtimsy.o dtimsp.o \
        dtimtr.o dtimtp.o dtimtb.o dtimhr.o dtimtd.o dtimbr.o \
        dtimq3.o dtimqr.o dtimlq.o dtimql.o dtimrq.o dtimqp.o \
        dqrt13.o dqrt15.o dlinpk.o deispk.o
 
ZLINTIM = ztimaa.o ztimmg.o \
        ztimmv.o ztimmm.o ztimb2.o ztimb3.o ztimge.o ztimgb.o ztimgt.o \
        ztimls.o ztimpo.o ztimpp.o ztimpb.o ztimpt.o \
        ztimsy.o ztimsp.o ztimhe.o ztimhp.o \
        ztimtr.o ztimtp.o ztimtb.o ztimhr.o ztimtd.o ztimbr.o \
        ztimq3.o ztimqr.o ztimlq.o ztimql.o ztimrq.o ztimqp.o \
        zqrt13.o zqrt15.o zlinpk.o
 
all:  single double complex complex16

single: ../xlintims
double: ../xlintimd
complex: ../xlintimc
complex16: ../xlintimz

../xlintims : $(ALINTIM) $(SCINTIM) $(SLINTIM) $(LINSRCLIB)
	$(LOADER) $(LOADOPTS)  $(ALINTIM) $(SCINTIM) $(SLINTIM) \
	../../$(TMGLIB) $(LINSRCLIB) ../../$(LAPACKLIB) $(BLASLIB) -o $@

../xlintimc : $(ALINTIM) $(SCINTIM) $(CLINTIM) $(LINSRCLIB)
	$(LOADER) $(LOADOPTS)  $(ALINTIM) $(SCINTIM) $(CLINTIM) \
	../../$(TMGLIB) $(LINSRCLIB) ../../$(LAPACKLIB) $(BLASLIB) -o $@

../xlintimd : $(ALINTIM) $(DZINTIM) $(DLINTIM) $(LINSRCLIB)
	$(LOADER) $(LOADOPTS)  $(ALINTIM) $(DZINTIM) $(DLINTIM) \
	../../$(TMGLIB) $(LINSRCLIB) ../../$(LAPACKLIB) $(BLASLIB)  -o $@

../xlintimz : $(ALINTIM) $(DZINTIM) $(ZLINTIM) $(LINSRCLIB)
	$(LOADER) $(LOADOPTS)  $(ALINTIM) $(DZINTIM) $(ZLINTIM) \
	../../$(TMGLIB) $(LINSRCLIB) ../../$(LAPACKLIB) $(BLASLIB)  -o $@

$(LINSRCLIB):
	cd LINSRC; $(MAKE) ../$(LINSRCLIB)

$(ALINTIM): $(FRC)
$(SCINTIM): $(FRC)
$(DZINTIM): $(FRC)
$(SLINTIM): $(FRC)
$(CLINTIM): $(FRC)
$(DLINTIM): $(FRC)
$(ZLINTIM): $(FRC)
 
FRC:
	@FRC=$(FRC)
 
clean:
	rm -f *.o

stimaa.o: stimaa.f
	$(FORTRAN) $(DRVOPTS) -c $<
dtimaa.o: dtimaa.f
	$(FORTRAN) $(DRVOPTS) -c $<
ctimaa.o: ctimaa.f
	$(FORTRAN) $(DRVOPTS) -c $<
ztimaa.o: ztimaa.f
	$(FORTRAN) $(DRVOPTS) -c $<
 
.f.o : ; $(FORTRAN) $(OPTS) -c $*.f
