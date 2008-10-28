include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
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

TEST_SRCS = do_test_spmv.c BLAS_spmv_testgen.c spmv-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS =  $(TEST_OBJS) $(TEST_DOT_OBJS) $(SYMV_OBJS) $(COMMON_OBJS)

test: do_test_spmv
	@echo Testing SPMV - default quick test
	./do_test_spmv 7 1 1.0 3 0.01 > spmv.results

test-few: do_test_spmv
	@echo Testing SPMV - preforming a small set of tests
	./do_test_spmv 7 1 1.0 3 $(DO_FEW_TESTS) > spmv.results

test-some: do_test_spmv
	@echo Testing SPMV - preforming some tests : not all tests 
	./do_test_spmv 7 1 1.0 3 $(DO_SOME_TESTS) > spmv.results

test-all: do_test_spmv
	@echo Testing SPMV Very Thoroughly
	./do_test_spmv 7 1 1.0 3 $(DO_ALL_TESTS) > spmv.results

test-none: do_test_spmv
	@echo NOT Testing SPMV - PREFORMING NO TESTS. NOT TESTING
	./do_test_spmv 7 1 1.0 3 $(DO_NO_TESTS) > spmv.results

do_test_spmv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_spmv $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK spmv.results do_test_spmv core

