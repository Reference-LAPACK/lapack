include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o \
		$(TEST_DOT_PATH)/print_vector.o \
		$(TEST_DOT_PATH)/copy_vector.o

TEST_SRCS = do_test_tpmv.c BLAS_tpmv_testgen.c tpmv-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(COMMON_OBJS)

test: do_test_tpmv
	@echo Testing TPMV - default quick test
	./do_test_tpmv 12 1 1.0 0 0.01 > tpmv.results

test-few: do_test_tpmv
	@echo Testing TPMV - preforming a small set of tests
	./do_test_tpmv 12 1 1.0 0 $(DO_FEW_TESTS)

test-some: do_test_tpmv
	@echo Testing TPMV - preforming some tests : not all tests 
	./do_test_tpmv 12 1 1.0 3 $(DO_SOME_TESTS)

test-all: do_test_tpmv
	@echo Testing TPMV Very Thoroughly
	./do_test_tpmv 12 1 1.0 0 $(DO_ALL_TESTS)

test-none: do_test_tpmv
	@echo NOT Testing TPMV - PREFORMING NO TESTS. NOT TESTING
	./do_test_tpmv 12 1 1.0 0 $(DO_NO_TESTS)

do_test_tpmv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_tpmv $(EXTRA_LIBS)

test-lib: $(TEST_OBJS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK tpmv.results do_test_tpmv core

