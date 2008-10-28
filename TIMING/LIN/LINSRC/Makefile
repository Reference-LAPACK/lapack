include ../../../make.inc

#######################################################################
#  This is the makefile to create a library of the least squares
#  routines from LAPACK that have been instrumented to time specific
#  parts of the code and count operations. 
#  The files are organized as follows:
#
#     SLINSRC -- Single precision real instrumented LAPACK routines
#     CLINSRC -- Single precision complex instrumented LAPACK routines
#     DLINSRC -- Double precision real instrumented LAPACK routines
#     ZLINSRC -- Double precision complex instrumented LAPACK routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       linsrc.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################

SCAUX = sopla.o sopla2.o

DZAUX = dopla.o dopla2.o
 
SLINSRC = \
   sgels.o  sgelsd.o sgelss.o sgelsx.o sgelsy.o slaic1.o slals0.o \
   slalsa.o slalsd.o

CLINSRC = \
   cgels.o  cgelsd.o cgelss.o cgelsx.o cgelsy.o claic1.o clals0.o \
   clalsa.o clalsd.o

DLINSRC = \
   dgels.o  dgelsd.o dgelss.o dgelsx.o dgelsy.o dlaic1.o dlals0.o \
   dlalsa.o dlalsd.o

ZLINSRC = \
   zgels.o  zgelsd.o zgelss.o zgelsx.o zgelsy.o zlaic1.o zlals0.o \
   zlalsa.o zlalsd.o

all: ../$(LINSRCLIB)
ALLOBJ=$(SLINSRC) $(CLINSRC) $(SCAUX) $(DLINSRC) $(ZLINSRC) $(DZAUX)
../$(LINSRCLIB): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

single: $(SLINSRC) $(SCAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LINSRCLIB) $(SLINSRC) $(SCAUX)
	$(RANLIB) ../$(LINSRCLIB)

complex: $(CLINSRC) $(SCAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LINSRCLIB) $(CLINSRC) $(SCAUX)
	$(RANLIB) ../$(LINSRCLIB)

double: $(DLINSRC) $(DZAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LINSRCLIB) $(DLINSRC) $(DZAUX)
	$(RANLIB) ../$(LINSRCLIB)

complex16: $(ZLINSRC) $(DZAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LINSRCLIB) $(ZLINSRC) $(DZAUX)
	$(RANLIB) ../$(LINSRCLIB)

$(SCAUX): $(FRC)
$(DZAUX): $(FRC)
$(SLINSRC): $(FRC)
$(CLINSRC): $(FRC)
$(DLINSRC): $(FRC)
$(ZLINSRC): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

.f.o: 
	$(FORTRAN) -c $(OPTS) $<
