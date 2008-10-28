include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
TEST_GEMM_PATH = ../test-gemm
TEST_SYMV_PATH = ../test-symv
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

SYMV_OBJS = $(TEST_SYMV_PATH)/symv-support.o \
		$(TEST_SYMV_PATH)/BLAS_symv_testgen.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o \
		$(TEST_DOT_PATH)/print_vector.o \
		$(TEST_DOT_PATH)/copy_vector.o

TEST_SRCS = do_test_sbmv.c BLAS_sbmv_testgen.c sbmv-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(GEMM_TEST_OBJS) \
		$(SYMV_OBJS) $(COMMON_OBJS)

test: do_test_sbmv
	@echo Testing SBMV - default quick test
	./do_test_sbmv 7 1 1.0 3 0.001 > sbmv.results

test-few: do_test_sbmv
	@echo Testing SBMV - preforming a small set of tests
	./do_test_sbmv 7 1 1.0 3 $(DO_FEW_TESTS) > sbmv.results

test-some: do_test_sbmv
	@echo Testing SBMV - preforming some tests : not all tests 
	./do_test_sbmv 7 1 1.0 3 $(DO_SOME_TESTS) > sbmv.results

test-all: do_test_sbmv
	@echo Testing SBMV Very Thoroughly
	./do_test_sbmv 7 1 1.0 3 $(DO_ALL_TESTS) > sbmv.results

test-none: do_test_sbmv
	@echo NOT Testing SBMV - PREFORMING NO TESTS. NOT TESTING
	./do_test_sbmv 7 1 1.0 3 $(DO_NO_TESTS) > sbmv.results

do_test_sbmv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_sbmv $(EXTRA_LIBS)

test-lib: $(TEST_OBJS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK sbmv.results do_test_sbmv core

