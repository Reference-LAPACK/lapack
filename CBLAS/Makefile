include ../make.inc

all: cblas

cblas: include/cblas_mangling.h
	$(MAKE) -C src

include/cblas_mangling.h: include/cblas_mangling_with_flags.h.in
	cp $< $@

cblas_testing: cblas
	$(MAKE) -C testing run

cblas_example: cblas
	$(MAKE) -C examples

clean: cleanlib

cleanlib:
	$(MAKE) -C src clean

cleanexe:
	$(MAKE) -C testing cleanexe

cleanall: clean cleanexe
	$(MAKE) -C src cleanall
	$(MAKE) -C examples cleanall
