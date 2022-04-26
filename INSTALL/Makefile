TOPSRCDIR = ..
include $(TOPSRCDIR)/make.inc

.PHONY: all testlsame testslamch testdlamch testsecond testdsecnd testieee testversion test_zcomplexabs test_zcomplexdiv test_zcomplexmult test_zminMax
all: testlsame testslamch testdlamch testsecond testdsecnd testieee testversion test_zcomplexabs test_zcomplexdiv test_zcomplexmult test_zminMax

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

testieee: tstiee.o $(TOPSRCDIR)/SRC/ieeeck.o $(TOPSRCDIR)/SRC/ilaenv.o $(TOPSRCDIR)/SRC/iparmq.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

testversion: ilaver.o LAPACK_version.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

test_zcomplexabs: test_zcomplexabs.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

test_zcomplexdiv: test_zcomplexdiv.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

test_zcomplexmult: test_zcomplexmult.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

test_zminMax: test_zminMax.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

.PHONY: run
run: all
	./testlsame
	./testslamch
	./testdlamch
	./testsecond
	./testdsecnd
	./testieee
	./testversion
	./test_zcomplexabs 2> test_zcomplexabs.err
	./test_zcomplexdiv 2> test_zcomplexdiv.err
	./test_zcomplexmult 2> test_zcomplexmult.err
	./test_zminMax 2> test_zminMax.err

.PHONY: clean cleanobj cleanexe cleantest
clean: cleanobj cleanexe cleantest
cleanobj:
	rm -f *.o
cleanexe:
	rm -f testlsame testslamch testdlamch testsecond testdsecnd testieee testversion test_zcomplexabs test_zcomplexdiv test_zcomplexmult test_zminMax
cleantest:
	rm -f core

slamch.o: slamch.f ; $(FC) $(FFLAGS_NOOPT) -c -o $@ $<
dlamch.o: dlamch.f ; $(FC) $(FFLAGS_NOOPT) -c -o $@ $<

sroundup_lwork.o: sroundup_lwork.f ; $(FC) $(FFLAGS_NOOPT) -c -o $@ $<
droundup_lwork.o: droundup_lwork.f ; $(FC) $(FFLAGS_NOOPT) -c -o $@ $<
