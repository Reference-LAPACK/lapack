TOPSRCDIR = ..
include $(TOPSRCDIR)/make.inc

all: testlsame testslamch testdlamch testsecond testdsecnd testieee testversion

testlsame: lsame.o lsametst.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

testslamch: slamch.o lsame.o slamchtst.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

testdlamch: dlamch.o lsame.o dlamchtst.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

testsecond: second_$(TIMER).o secondtst.o
	@echo "[INFO] : TIMER value: $(TIMER) (given by make.inc)"
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

testdsecnd: dsecnd_$(TIMER).o dsecndtst.o
	@echo "[INFO] : TIMER value: $(TIMER) (given by make.inc)"
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

testieee: tstiee.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

testversion: ilaver.o LAPACK_version.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

run: all
	./testlsame
	./testslamch
	./testdlamch
	./testsecond
	./testdsecnd
	./testieee
	./testversion

clean: cleanobj cleanexe cleantest
cleanobj:
	rm -f *.o
cleanexe:
	rm -f test*
cleantest:
	rm -f core

slamch.o: slamch.f ; $(FC) $(FFLAGS_NOOPT) -c -o $@ $<
dlamch.o: dlamch.f ; $(FC) $(FFLAGS_NOOPT) -c -o $@ $<
