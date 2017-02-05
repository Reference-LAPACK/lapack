include ../../make.inc

all: cblas_ex1 cblas_ex2

cblas_ex1: cblas_example1.o $(CBLASLIB) $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^
cblas_ex2: cblas_example2.o $(CBLASLIB) $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^

cleanall:
	rm -f *.o cblas_ex1 cblas_ex2

.c.o:
	$(CC) $(CFLAGS) -I../include -c -o $@ $<
