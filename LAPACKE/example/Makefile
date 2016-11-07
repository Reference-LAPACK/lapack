include ../../make.inc

all: xexample_DGESV_rowmajor \
     xexample_DGESV_colmajor \
     xexample_DGELS_rowmajor \
     xexample_DGELS_colmajor

LIBRARIES = ../../$(LAPACKELIB) ../../$(LAPACKLIB) $(BLASLIB)

# Double Precision Examples
xexample_DGESV_rowmajor: example_DGESV_rowmajor.o lapacke_example_aux.o $(LIBRARIES)
	$(LOADER) $(LOADOPTS) example_DGESV_rowmajor.o lapacke_example_aux.o \
	$(LIBRARIES) -o $@
	./$@

xexample_DGESV_colmajor: example_DGESV_colmajor.o lapacke_example_aux.o $(LIBRARIES)
	$(LOADER) $(LOADOPTS) example_DGESV_colmajor.o lapacke_example_aux.o \
	$(LIBRARIES) -o $@
	./$@

xexample_DGELS_rowmajor: example_DGELS_rowmajor.o lapacke_example_aux.o $(LIBRARIES)
	$(LOADER) $(LOADOPTS) example_DGELS_rowmajor.o lapacke_example_aux.o \
	$(LIBRARIES) -o $@
	./$@

xexample_DGELS_colmajor: example_DGELS_colmajor.o lapacke_example_aux.o $(LIBRARIES)
	$(LOADER) $(LOADOPTS) example_DGELS_colmajor.o lapacke_example_aux.o \
	$(LIBRARIES) -o $@
	./$@

.c.o:
	$(CC) -c $(CFLAGS) -I. -I ../include -o $@ $<

clean:
	rm -f *.o x*
