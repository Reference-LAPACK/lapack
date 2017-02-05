include ../make.inc

all: cblas

cblas: include/cblas_mangling.h
	cd src && $(MAKE)

include/cblas_mangling.h: include/cblas_mangling_with_flags.h.in
	cp $< $@

cblas_testing: cblas
	cd testing && $(MAKE) run

cblas_example: cblas
	cd examples && $(MAKE)

clean: cleanlib

cleanlib:
	cd src && $(MAKE) clean

cleanexe:
	cd testing && $(MAKE) cleanexe

cleanall: clean cleanexe
	cd src && $(MAKE) cleanall
	cd examples && $(MAKE) cleanall
