dlvl = ./.
include $(dlvl)/Makefile.in

all: alllib alltst 

help:
	@ echo "Make sure you are using correct Makefile.in for your system."
	@ echo "At this level, assuming you have downloded all necessary    "
	@ echo "files and made an archive file of BLAS routines for your    "
	@ echo "system."
	@ echo " "
	@ echo "The Makefile compiles the routines of CBLAS (C interface of "
	@ echo "BLAS) and testers for all the precisions.                   "
	@ echo "If there is no directory for archives in CBLAS/lib, it      "
	@ echo "creates new directory with the name of the platform of your "
	@ echo "machine." 
	@ echo " "
	@ echo "To compile, you have to type as follows"
	@ echo "make <target>"
	@ echo " where <target> is one of:"
	@ echo "slib1 --- make an archive of level 1 REAL."
	@ echo "dlib1 --- make an archive of level 1 DOUBLE PRECISION."
	@ echo "clib1 --- make an archive of level 1 COMPLEX."
	@ echo "zlib1 --- make an archive of level 1 COMPLEX*16."
	@ echo "alllib1 - make an archive of level 1 all precisions."
	@ echo " "
	@ echo "slib2 --- make an archive of level 2 REAL."
	@ echo "dlib2 --- make an archive of level 2 DOUBLE PRECSION."
	@ echo "clib2 --- make an archive of level 2 COMPLEX."
	@ echo "zlib2 --- make an archive of level 2 COMPLEX*16."
	@ echo "alllib2 - make an archive of level 2 all precisions."
	@ echo " "
	@ echo "slib3 --- make an archive of level 3 REAL."
	@ echo "dlib3 --- make an archive of level 3 DOUBLE PRECISION ."
	@ echo "clib3 --- make an archive of level 3 COMPLEX."
	@ echo "zlib3 --- make an archive of level 3 COMPLEX*16."
	@ echo "alllib3 - make an archive of level 3 all precisions."
	@ echo " "
	@ echo "alllib -- make an archive for all precisions."
	@ echo " "
	@ echo "stest1 -- Compiles the tester for level 1 REAL."
	@ echo "dtest1 -- Compiles the tester for level 1 DOUBLE PRECISION. "
	@ echo "ctest1 -- Compiles the tester for level 1 COMPLEX."
	@ echo "ztest1 -- Compiles the tester for level 1 COMPLEX*16."
	@ echo "alltst1 - Compiles testers for all precisions of level 1." 
	@ echo " "
	@ echo "stest2 -- Compiles the tester for level 2 REAL."
	@ echo "dtest2 -- Compiles the tester for level 2 DOUBLE PRECISION. "
	@ echo "ctest2 -- Compiles the tester for level 2 COMPLEX."
	@ echo "ztest2 -- Compiles the tester for level 2 COMPLEX*16."
	@ echo "alltst2 - Compiles testers for all precisions of level 2." 
	@ echo " "
	@ echo "stest3 -- Compiles the tester for level 3 REAL."
	@ echo "dtest3 -- Compiles the tester for level 3 DOUBLE PRECISON. "
	@ echo "ctest3 -- Compiles the tester for level 3 COMPLEX."
	@ echo "ztest3 -- Compiles the tester for level 3 COMPLEX*16."
	@ echo "alltst3 - Compiles testers for all precisions of level 3." 
	@ echo " "
	@ echo "alltst -- Compiles testers for all CBLAS routines." 
	@ echo "runtst -- Execute testers for all CBLAS routines." 
	@ echo " "
	@ echo "all ----- Creates a library and testers for ALL." 
	@ echo " "
	@ echo "clean --- Erase all the .o and excutable files" 
	@ echo "cleanlib -- Erase all the .o  files" 
	@ echo "cleanexe -- Erase all the excutable files" 
	@ echo "rmlib --- Remove a library file." 
	@ echo " "
	@ echo "example -- Creates example1 and example2"
	@ echo "example1 -- A small example to exercise the interface "
	@ echo "example2 -- Test that cblas_xerbla() is working correctly"
	@ echo " "
	@ echo " ------- Warning ------- "
	@ echo "If you want just to make a tester, make sure you have"
	@ echo "already made an archive file out of CBLAS routines."
	@ echo " "
	@ echo "Written by Keita Teranishi"
	@ echo "3/4/98 "


# In general, the Makefile call other Makefiles in the sub-directories.


clean:
	( cd testing && make clean )
	( cd src && make clean )
	rm -f *.o cblas_ex1 cblas_ex2

cleanobj:
	( cd testing && make cleanobj )
	( cd src && make clean )

cleanexe:
	( cd testing && make cleanexe )

rmlib:
	( rm -f $(CBLIB) )
slib1:  sreal1
dlib1:  dreal1
clib1:  scplx1
zlib1:  dcplx1
slib2:  sreal2
dlib2:  dreal2
clib2:  scplx2
zlib2:  dcplx2
slib3:  sreal3
dlib3:  dreal3
clib3:  scplx3 
zlib3:  dcplx3 
alllib1: allprecision1
alllib2: allprecision2
alllib3: allprecision3
alllib:  allprecision


sreal1:
	( cd src && make slib1)
dreal1:
	( cd src && make dlib1)
scplx1:
	( cd src && make clib1)
dcplx1:
	( cd src && make zlib1)
allprecision1:
	( cd src && make all1)
sreal2:
	( cd src && make slib2)
dreal2:
	( cd src && make dlib2)
scplx2:
	( cd src && make clib2)
dcplx2:
	( cd src && make zlib2)
allprecision2:
	( cd src && make all2)
sreal3:
	( cd src && make slib3)
dreal3:
	( cd src && make dlib3)
scplx3:
	( cd src && make clib3)
dcplx3:
	( cd src && make zlib3)
allprecision3:
	( cd src && make all3)
allprecision:
	( cd src && make all)

stest1: 
	( cd testing && make stest1 )
dtest1: 
	( cd testing && make dtest1 )
ctest1: 
	( cd testing && make ctest1 )
ztest1: 
	( cd testing && make ztest1 )
alltst1:
	( cd testing && make all1 )
stest2:
	( cd testing && make stest2 )
dtest2:
	( cd testing && make dtest2 )
ctest2:
	( cd testing && make ctest2 )
ztest2:
	( cd testing && make ztest2 )
alltst2:
	( cd testing && make all2 )
stest3:
	( cd testing && make stest3 )
dtest3:
	( cd testing && make dtest3 )
ctest3:
	( cd testing && make ctest3 )
ztest3:
	( cd testing && make ztest3 )
alltst3:
	( cd testing && make all3 )
alltst:
	( cd testing && make all )
runtst:
	( cd testing && make run )
	
example: alllib
	( cd examples && make all )
example1: alllib
	( cd examples && make example1 )
example2: alllib
	( cd examples && make example1 )

   
cleanall:
	( cd src && rm -f a.out core *.o $(CBLIB) )
	( cd testing && rm -f *.out core *.o x[sdcz]cblat[123] )
	( cd examples && rm -f *.o cblas_ex1 cblas_ex2 )
