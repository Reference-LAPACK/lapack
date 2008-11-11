# Top-level makefile for XBLAS.
#
# To generate all libraries and test results, just type make
#
# To generate the library, type make lib
#
# To generate all the library, and the test library, type make test-lib
#
# To clean out all the object files, type make clean
#
include make.conf
include $(MAKEINC)

LIB = dot-lib sum-lib axpby-lib waxpby-lib gemv-lib ge_sum_mv-lib gbmv-lib symv-lib spmv-lib sbmv-lib hemv-lib hpmv-lib hbmv-lib trmv-lib tpmv-lib trsv-lib tbsv-lib gemm-lib symm-lib hemm-lib gemv2-lib symv2-lib
TESTS = dot-test sum-test axpby-test waxpby-test gemv-test ge_sum_mv-test gbmv-test symv-test spmv-test sbmv-test hemv-test hpmv-test hbmv-test trmv-test tpmv-test trsv-test tbsv-test gemm-test symm-test hemm-test gemv2-test symv2-test
TEST_LIB = dot-test-lib sum-test-lib axpby-test-lib waxpby-test-lib gemv-test-lib ge_sum_mv-test-lib gbmv-test-lib symv-test-lib spmv-test-lib sbmv-test-lib hemv-test-lib hpmv-test-lib hbmv-test-lib trmv-test-lib tpmv-test-lib trsv-test-lib tbsv-test-lib gemm-test-lib symm-test-lib hemm-test-lib gemv2-test-lib symv2-test-lib
CLEAN = dot-clean sum-clean axpby-clean waxpby-clean gemv-clean ge_sum_mv-clean gbmv-clean symv-clean spmv-clean sbmv-clean hemv-clean hpmv-clean hbmv-clean trmv-clean tpmv-clean trsv-clean tbsv-clean gemm-clean symm-clean hemm-clean gemv2-clean symv2-clean dot2-clean

SRC_DIR = src
TEST_DIR = testing

all: tests

$(TEST_LIB): lib

test-lib: lib $(TEST_LIB)

$(TESTS): test-lib

tests: test-lib $(TESTS)
	rm -f testall.result testall.summary
	cat $(TEST_DIR)/test-dot/dot.results >> testall.result
	cat $(TEST_DIR)/test-sum/sum.results >> testall.result
	cat $(TEST_DIR)/test-axpby/axpby.results >> testall.result
	cat $(TEST_DIR)/test-waxpby/waxpby.results >> testall.result
	cat $(TEST_DIR)/test-gemv/gemv.results >> testall.result
	cat $(TEST_DIR)/test-ge_sum_mv/ge_sum_mv.results >> testall.result
	cat $(TEST_DIR)/test-gbmv/gbmv.results >> testall.result
	cat $(TEST_DIR)/test-symv/symv.results >> testall.result
	cat $(TEST_DIR)/test-spmv/spmv.results >> testall.result
	cat $(TEST_DIR)/test-sbmv/sbmv.results >> testall.result
	cat $(TEST_DIR)/test-hemv/hemv.results >> testall.result
	cat $(TEST_DIR)/test-hpmv/hpmv.results >> testall.result
	cat $(TEST_DIR)/test-hbmv/hbmv.results >> testall.result
	cat $(TEST_DIR)/test-trmv/trmv.results >> testall.result
	cat $(TEST_DIR)/test-tpmv/tpmv.results >> testall.result
	cat $(TEST_DIR)/test-trsv/trsv.results >> testall.result
	cat $(TEST_DIR)/test-tbsv/tbsv.results >> testall.result
	cat $(TEST_DIR)/test-gemm/gemm.results >> testall.result
	cat $(TEST_DIR)/test-symm/symm.results >> testall.result
	cat $(TEST_DIR)/test-hemm/hemm.results >> testall.result
	cat $(TEST_DIR)/test-gemv2/gemv2.results >> testall.result
	cat $(TEST_DIR)/test-symv2/symv2.results >> testall.result
	grep 'FAIL/TOTAL' testall.result >testall.summary
	cat testall.summary

common-lib:
	cd $(SRC_DIR)/common && $(MAKE)

common-test-lib:
	cd $(TEST_DIR)/common && $(MAKE)

lib: $(LIB)
	cd $(SRC_DIR)/common && $(MAKE) lib
	cd $(SRC_DIR)/dot && $(MAKE) lib
	cd $(SRC_DIR)/sum && $(MAKE) lib
	cd $(SRC_DIR)/axpby && $(MAKE) lib
	cd $(SRC_DIR)/waxpby && $(MAKE) lib
	cd $(SRC_DIR)/gemv && $(MAKE) lib
	cd $(SRC_DIR)/ge_sum_mv && $(MAKE) lib
	cd $(SRC_DIR)/gbmv && $(MAKE) lib
	cd $(SRC_DIR)/symv && $(MAKE) lib
	cd $(SRC_DIR)/spmv && $(MAKE) lib
	cd $(SRC_DIR)/sbmv && $(MAKE) lib
	cd $(SRC_DIR)/hemv && $(MAKE) lib
	cd $(SRC_DIR)/hpmv && $(MAKE) lib
	cd $(SRC_DIR)/hbmv && $(MAKE) lib
	cd $(SRC_DIR)/trmv && $(MAKE) lib
	cd $(SRC_DIR)/tpmv && $(MAKE) lib
	cd $(SRC_DIR)/trsv && $(MAKE) lib
	cd $(SRC_DIR)/tbsv && $(MAKE) lib
	cd $(SRC_DIR)/gemm && $(MAKE) lib
	cd $(SRC_DIR)/symm && $(MAKE) lib
	cd $(SRC_DIR)/hemm && $(MAKE) lib
	cd $(SRC_DIR)/gemv2 && $(MAKE) lib
	cd $(SRC_DIR)/symv2 && $(MAKE) lib

# custom test-dot2 stuff

dot2-test-lib: lib
	cd $(TEST_DIR)/test-dot2 && $(MAKE) do_test_dot2 

dot2-test: dot2-test-lib
	cd $(TEST_DIR)/test-dot2 && $(MAKE)

dot2-clean:
	cd $(TEST_DIR)/test-dot2 && $(MAKE) clean


# dot stuff

dot: dot-test


dot-lib: common-lib
	cd $(SRC_DIR)/dot && $(MAKE)

dot-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-dot && $(MAKE) do_test_dot 

dot-test: dot-test-lib
	cd $(TEST_DIR)/test-dot && $(MAKE)

dot-clean:
	cd $(SRC_DIR)/dot && $(MAKE) clean
	cd $(TEST_DIR)/test-dot && $(MAKE) clean


# sum stuff

sum: sum-test


sum-lib: common-lib
	cd $(SRC_DIR)/sum && $(MAKE)

sum-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-sum && $(MAKE) do_test_sum 

sum-test: sum-test-lib
	cd $(TEST_DIR)/test-sum && $(MAKE)

sum-clean:
	cd $(SRC_DIR)/sum && $(MAKE) clean
	cd $(TEST_DIR)/test-sum && $(MAKE) clean


# axpby stuff

axpby: axpby-test


axpby-lib: common-lib
	cd $(SRC_DIR)/axpby && $(MAKE)

axpby-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-axpby && $(MAKE) do_test_axpby 

axpby-test: axpby-test-lib
	cd $(TEST_DIR)/test-axpby && $(MAKE)

axpby-clean:
	cd $(SRC_DIR)/axpby && $(MAKE) clean
	cd $(TEST_DIR)/test-axpby && $(MAKE) clean


# waxpby stuff

waxpby: waxpby-test


waxpby-lib: common-lib
	cd $(SRC_DIR)/waxpby && $(MAKE)

waxpby-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-waxpby && $(MAKE) do_test_waxpby 

waxpby-test: waxpby-test-lib
	cd $(TEST_DIR)/test-waxpby && $(MAKE)

waxpby-clean:
	cd $(SRC_DIR)/waxpby && $(MAKE) clean
	cd $(TEST_DIR)/test-waxpby && $(MAKE) clean


# gemv stuff

gemv: gemv-test


gemv-lib: common-lib
	cd $(SRC_DIR)/gemv && $(MAKE)

gemv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-gemv && $(MAKE) do_test_gemv 

gemv-test: gemv-test-lib
	cd $(TEST_DIR)/test-gemv && $(MAKE)

gemv-clean:
	cd $(SRC_DIR)/gemv && $(MAKE) clean
	cd $(TEST_DIR)/test-gemv && $(MAKE) clean


# ge_sum_mv stuff

ge_sum_mv: ge_sum_mv-test


ge_sum_mv-lib: common-lib
	cd $(SRC_DIR)/ge_sum_mv && $(MAKE)

ge_sum_mv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-ge_sum_mv && $(MAKE) do_test_ge_sum_mv 

ge_sum_mv-test: ge_sum_mv-test-lib
	cd $(TEST_DIR)/test-ge_sum_mv && $(MAKE)

ge_sum_mv-clean:
	cd $(SRC_DIR)/ge_sum_mv && $(MAKE) clean
	cd $(TEST_DIR)/test-ge_sum_mv && $(MAKE) clean


# gbmv stuff

gbmv: gbmv-test


gbmv-lib: common-lib
	cd $(SRC_DIR)/gbmv && $(MAKE)

gbmv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-gbmv && $(MAKE) do_test_gbmv 

gbmv-test: gbmv-test-lib
	cd $(TEST_DIR)/test-gbmv && $(MAKE)

gbmv-clean:
	cd $(SRC_DIR)/gbmv && $(MAKE) clean
	cd $(TEST_DIR)/test-gbmv && $(MAKE) clean


# symv stuff

symv: symv-test


symv-lib: common-lib
	cd $(SRC_DIR)/symv && $(MAKE)

symv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-symv && $(MAKE) do_test_symv 

symv-test: symv-test-lib
	cd $(TEST_DIR)/test-symv && $(MAKE)

symv-clean:
	cd $(SRC_DIR)/symv && $(MAKE) clean
	cd $(TEST_DIR)/test-symv && $(MAKE) clean


# spmv stuff

spmv: spmv-test


spmv-lib: common-lib
	cd $(SRC_DIR)/spmv && $(MAKE)

spmv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-spmv && $(MAKE) do_test_spmv 

spmv-test: spmv-test-lib
	cd $(TEST_DIR)/test-spmv && $(MAKE)

spmv-clean:
	cd $(SRC_DIR)/spmv && $(MAKE) clean
	cd $(TEST_DIR)/test-spmv && $(MAKE) clean


# sbmv stuff

sbmv: sbmv-test


sbmv-lib: common-lib
	cd $(SRC_DIR)/sbmv && $(MAKE)

sbmv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-sbmv && $(MAKE) do_test_sbmv 

sbmv-test: sbmv-test-lib
	cd $(TEST_DIR)/test-sbmv && $(MAKE)

sbmv-clean:
	cd $(SRC_DIR)/sbmv && $(MAKE) clean
	cd $(TEST_DIR)/test-sbmv && $(MAKE) clean


# hemv stuff

hemv: hemv-test


hemv-lib: common-lib
	cd $(SRC_DIR)/hemv && $(MAKE)

hemv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-hemv && $(MAKE) do_test_hemv 

hemv-test: hemv-test-lib
	cd $(TEST_DIR)/test-hemv && $(MAKE)

hemv-clean:
	cd $(SRC_DIR)/hemv && $(MAKE) clean
	cd $(TEST_DIR)/test-hemv && $(MAKE) clean


# hpmv stuff

hpmv: hpmv-test


hpmv-lib: common-lib
	cd $(SRC_DIR)/hpmv && $(MAKE)

hpmv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-hpmv && $(MAKE) do_test_hpmv 

hpmv-test: hpmv-test-lib
	cd $(TEST_DIR)/test-hpmv && $(MAKE)

hpmv-clean:
	cd $(SRC_DIR)/hpmv && $(MAKE) clean
	cd $(TEST_DIR)/test-hpmv && $(MAKE) clean


# hbmv stuff

hbmv: hbmv-test


hbmv-lib: common-lib
	cd $(SRC_DIR)/hbmv && $(MAKE)

hbmv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-hbmv && $(MAKE) do_test_hbmv 

hbmv-test: hbmv-test-lib
	cd $(TEST_DIR)/test-hbmv && $(MAKE)

hbmv-clean:
	cd $(SRC_DIR)/hbmv && $(MAKE) clean
	cd $(TEST_DIR)/test-hbmv && $(MAKE) clean


# trmv stuff

trmv: trmv-test


trmv-lib: common-lib
	cd $(SRC_DIR)/trmv && $(MAKE)

trmv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-trmv && $(MAKE) do_test_trmv 

trmv-test: trmv-test-lib
	cd $(TEST_DIR)/test-trmv && $(MAKE)

trmv-clean:
	cd $(SRC_DIR)/trmv && $(MAKE) clean
	cd $(TEST_DIR)/test-trmv && $(MAKE) clean


# tpmv stuff

tpmv: tpmv-test


tpmv-lib: common-lib
	cd $(SRC_DIR)/tpmv && $(MAKE)

tpmv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-tpmv && $(MAKE) do_test_tpmv 

tpmv-test: tpmv-test-lib
	cd $(TEST_DIR)/test-tpmv && $(MAKE)

tpmv-clean:
	cd $(SRC_DIR)/tpmv && $(MAKE) clean
	cd $(TEST_DIR)/test-tpmv && $(MAKE) clean


# trsv stuff

trsv: trsv-test


trsv-lib: common-lib
	cd $(SRC_DIR)/trsv && $(MAKE)

trsv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-trsv && $(MAKE) do_test_trsv 

trsv-test: trsv-test-lib
	cd $(TEST_DIR)/test-trsv && $(MAKE)

trsv-clean:
	cd $(SRC_DIR)/trsv && $(MAKE) clean
	cd $(TEST_DIR)/test-trsv && $(MAKE) clean


# tbsv stuff

tbsv: tbsv-test


tbsv-lib: common-lib
	cd $(SRC_DIR)/tbsv && $(MAKE)

tbsv-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-tbsv && $(MAKE) do_test_tbsv 

tbsv-test: tbsv-test-lib
	cd $(TEST_DIR)/test-tbsv && $(MAKE)

tbsv-clean:
	cd $(SRC_DIR)/tbsv && $(MAKE) clean
	cd $(TEST_DIR)/test-tbsv && $(MAKE) clean


# gemm stuff

gemm: gemm-test


gemm-lib: common-lib
	cd $(SRC_DIR)/gemm && $(MAKE)

gemm-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-gemm && $(MAKE) do_test_gemm 

gemm-test: gemm-test-lib
	cd $(TEST_DIR)/test-gemm && $(MAKE)

gemm-clean:
	cd $(SRC_DIR)/gemm && $(MAKE) clean
	cd $(TEST_DIR)/test-gemm && $(MAKE) clean


# symm stuff

symm: symm-test


symm-lib: common-lib
	cd $(SRC_DIR)/symm && $(MAKE)

symm-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-symm && $(MAKE) do_test_symm 

symm-test: symm-test-lib
	cd $(TEST_DIR)/test-symm && $(MAKE)

symm-clean:
	cd $(SRC_DIR)/symm && $(MAKE) clean
	cd $(TEST_DIR)/test-symm && $(MAKE) clean


# hemm stuff

hemm: hemm-test


hemm-lib: common-lib
	cd $(SRC_DIR)/hemm && $(MAKE)

hemm-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-hemm && $(MAKE) do_test_hemm 

hemm-test: hemm-test-lib
	cd $(TEST_DIR)/test-hemm && $(MAKE)

hemm-clean:
	cd $(SRC_DIR)/hemm && $(MAKE) clean
	cd $(TEST_DIR)/test-hemm && $(MAKE) clean


# gemv2 stuff

gemv2: gemv2-test


gemv2-lib: common-lib
	cd $(SRC_DIR)/gemv2 && $(MAKE)

gemv2-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-gemv2 && $(MAKE) do_test_gemv2 

gemv2-test: gemv2-test-lib
	cd $(TEST_DIR)/test-gemv2 && $(MAKE)

gemv2-clean:
	cd $(SRC_DIR)/gemv2 && $(MAKE) clean
	cd $(TEST_DIR)/test-gemv2 && $(MAKE) clean


# symv2 stuff

symv2: symv2-test


symv2-lib: common-lib
	cd $(SRC_DIR)/symv2 && $(MAKE)

symv2-test-lib: lib common-test-lib
	cd $(TEST_DIR)/test-symv2 && $(MAKE) do_test_symv2 

symv2-test: symv2-test-lib
	cd $(TEST_DIR)/test-symv2 && $(MAKE)

symv2-clean:
	cd $(SRC_DIR)/symv2 && $(MAKE) clean
	cd $(TEST_DIR)/test-symv2 && $(MAKE) clean



# Test library dependencies
sum-test-lib: dot-test-lib
axpby-test-lib: dot-test-lib
waxpby-test-lib: dot-test-lib
gemv-test-lib: dot-test-lib
ge_sum_mv-test-lib: dot-test-lib gemv-test-lib symv-test-lib gemm-test-lib
gbmv-test-lib: dot-test-lib
symv-test-lib: dot-test-lib gemm-test-lib
spmv-test-lib: dot-test-lib symv-test-lib
sbmv-test-lib: dot-test-lib symv-test-lib
hemv-test-lib: dot-test-lib symv-test-lib
hbmv-test-lib: dot-test-lib symv-test-lib hemv-test-lib sbmv-test-lib
hpmv-test-lib: dot-test-lib symv-test-lib hemv-test-lib
trmv-test-lib: dot-test-lib
tpmv-test-lib: dot-test-lib
trsv-test-lib: dot-test-lib trmv-test-lib
tbsv-test-lib: dot-test-lib gbmv-test-lib trsv-test-lib
gemm-test-lib: dot-test-lib gemv-test-lib
symm-test-lib: dot-test-lib gemv-test-lib symv-test-lib
hemm-test-lib: dot-test-lib gemv-test-lib hemv-test-lib symv-test-lib symm-test-lib
gemv2-test-lib: dot2-test-lib gemv-test-lib
dot2-test-lib: dot-test-lib


# Cleaning stuff

clean: $(CLEAN)
	cd $(SRC_DIR)/common && $(MAKE) clean
	cd $(TEST_DIR)/common && $(MAKE) clean

dist-clean: clean
	rm -f $(XBLASLIB) testall.result testall.summary

.PHONY: $(LIBS) $(TEST_LIB) $(TESTS) $(CLEAN) \
        all test-lib tests common-lib common-test-lib lib clean dist-clean

