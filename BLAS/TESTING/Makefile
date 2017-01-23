include ../../make.inc

all: single double complex complex16
single:    xblat1s xblat2s xblat3s
double:    xblat1d xblat2d xblat3d
complex:   xblat1c xblat2c xblat3c
complex16: xblat1z xblat2z xblat3z

xblat1s: sblat1.o
	$(LOADER) $(LOADOPTS) -o $@ sblat1.o $(BLASLIB)
xblat1d: dblat1.o
	$(LOADER) $(LOADOPTS) -o $@ dblat1.o $(BLASLIB)
xblat1c: cblat1.o
	$(LOADER) $(LOADOPTS) -o $@ cblat1.o $(BLASLIB)
xblat1z: zblat1.o
	$(LOADER) $(LOADOPTS) -o $@ zblat1.o $(BLASLIB)

xblat2s: sblat2.o
	$(LOADER) $(LOADOPTS) -o $@ sblat2.o $(BLASLIB)
xblat2d: dblat2.o
	$(LOADER) $(LOADOPTS) -o $@ dblat2.o $(BLASLIB)
xblat2c: cblat2.o
	$(LOADER) $(LOADOPTS) -o $@ cblat2.o $(BLASLIB)
xblat2z: zblat2.o
	$(LOADER) $(LOADOPTS) -o $@ zblat2.o $(BLASLIB)

xblat3s: sblat3.o
	$(LOADER) $(LOADOPTS) -o $@ sblat3.o $(BLASLIB)
xblat3d: dblat3.o
	$(LOADER) $(LOADOPTS) -o $@ dblat3.o $(BLASLIB)
xblat3c: cblat3.o
	$(LOADER) $(LOADOPTS) -o $@ cblat3.o $(BLASLIB)
xblat3z: zblat3.o
	$(LOADER) $(LOADOPTS) -o $@ zblat3.o $(BLASLIB)

run: all
	./xblat1s > sblat1.out
	./xblat1d > dblat1.out
	./xblat1c > cblat1.out
	./xblat1z > zblat1.out
	./xblat2s < sblat2.in
	./xblat2d < dblat2.in
	./xblat2c < cblat2.in
	./xblat2z < zblat2.in
	./xblat3s < sblat3.in
	./xblat3d < dblat3.in
	./xblat3c < cblat3.in
	./xblat3z < zblat3.in

clean:
	rm -f *.o xblat* *.out

.f.o:
	$(FORTRAN) $(OPTS) -c -o $@ $<
