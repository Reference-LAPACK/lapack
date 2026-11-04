TOPSRCDIR = ..
include $(TOPSRCDIR)/make.inc

.PHONY: all
all: cblas

CBLAS_INCLUDES = include/cblas.h include/cblas_mangling.h include/cblas_f77.h

.PHONY: cblas
cblas: $(CBLAS_INCLUDES)
	$(MAKE) -C src

include/cblas.h: include/cblas.h.in
	sed 's/@FORTRAN_STRLEN_TYPE@/size_t/' $< > $@.tmp && mv $@.tmp $@

include/cblas_mangling.h: include/cblas_mangling_with_flags.h.in
	cp $< $@

include/cblas_f77.h: include/cblas_f77.h.in
	sed 's/@FORTRAN_STRLEN_TYPE@/size_t/' $< > $@.tmp && mv $@.tmp $@

.PHONY: cblas_testing
cblas_testing: cblas
	$(MAKE) -C testing run

.PHONY: cblas_example
cblas_example: cblas
	$(MAKE) -C examples

.PHONY: clean cleanobj cleanlib cleanexe cleantest
clean:
	$(MAKE) -C src clean
	$(MAKE) -C testing clean
	$(MAKE) -C examples clean
	rm -f $(CBLAS_INCLUDES)
cleanobj:
	$(MAKE) -C src cleanobj
	$(MAKE) -C testing cleanobj
	$(MAKE) -C examples cleanobj
cleanlib:
	$(MAKE) -C src cleanlib
cleanexe:
	$(MAKE) -C testing cleanexe
	$(MAKE) -C examples cleanexe
cleantest:
	$(MAKE) -C testing cleantest
