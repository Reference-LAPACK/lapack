include ../../make.inc

all: test_nancheck_DGELS

LIBRARIES = ../../$(LAPACKELIB) ../../$(LAPACKLIB) $(BLASLIB)

# Double Precision Examples

test_nancheck_DGELS: test_nancheck_DGELS.o lapacke_example_aux.o $(LIBRARIES)
	$(LOADER) $(LOADOPTS) -o $@ $^
	./$@

clean: cleanobj cleanexe
cleanobj:
	rm -f *.o
cleanexe:
	rm -f x*

.c.o:
	$(CC) $(CFLAGS) -I. -I../include -c -o $@ $<
