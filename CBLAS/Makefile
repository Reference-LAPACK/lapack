include ../make.inc

all:
	cd include && cp cblas_mangling_with_flags.h.in cblas_mangling.h
	cd src && $(MAKE) all


clean: cleanlib

cleanlib:
	cd src && $(MAKE) clean

cleanexe:
	cd testing && $(MAKE) cleanexe

cleanall: clean cleanexe
	cd src && $(MAKE) cleanall
	cd examples && $(MAKE) cleanall

cblas_testing:
	cd testing && $(MAKE) all

runtst:
	cd testing && $(MAKE) run

example: all
	cd examples && $(MAKE) all
