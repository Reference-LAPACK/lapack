include ../../make.inc

all: cblas_ex1 cblas_ex2

cblas_ex1: cblas_example1.o
	$(LOADER) $(LOADOPTS) -o $@ cblas_example1.o $(CBLASLIB) $(BLASLIB)
cblas_ex2: cblas_example2.o
	$(LOADER) $(LOADOPTS) -o $@ cblas_example2.o $(CBLASLIB) $(BLASLIB)

cleanall:
	rm -f *.o cblas_ex1 cblas_ex2

.c.o:
	$(CC) $(CFLAGS) -I../include -c -o $@ $<
