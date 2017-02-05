include ../make.inc

all: testlsame testslamch testdlamch testsecond testdsecnd testieee testversion

clean:
	rm -f *.o
cleanall:
	rm -f test*

testlsame: lsame.o lsametst.o
	$(LOADER) $(LOADOPTS) -o $@ $^

testslamch: slamch.o lsame.o slamchtst.o
	$(LOADER) $(LOADOPTS) -o $@ $^

testdlamch: dlamch.o lsame.o dlamchtst.o
	$(LOADER) $(LOADOPTS) -o $@ $^

testsecond: second_$(TIMER).o secondtst.o
	@echo "[INFO] : TIMER value: $(TIMER) (given by make.inc)"
	$(LOADER) $(LOADOPTS) -o $@ $^

testdsecnd: dsecnd_$(TIMER).o dsecndtst.o
	@echo "[INFO] : TIMER value: $(TIMER) (given by make.inc)"
	$(LOADER) $(LOADOPTS) -o $@ $^

testieee: tstiee.o
	$(LOADER) $(LOADOPTS) -o $@ $^

testversion: ilaver.o LAPACK_version.o
	$(LOADER) $(LOADOPTS) -o $@ $^

run: all
	./testlsame
	./testslamch
	./testdlamch
	./testsecond
	./testdsecnd
	./testieee
	./testversion

.SUFFIXES: .o .f
.f.o:
	$(FORTRAN) $(OPTS) -c -o $@ $<

slamch.o: slamch.f ; $(FORTRAN) $(NOOPT) -c -o $@ $<
dlamch.o: dlamch.f ; $(FORTRAN) $(NOOPT) -c -o $@ $<
