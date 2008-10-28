include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
TEST_GEMV_PATH = ../test-gemv
TEST_SYMV_PATH = ../test-symv
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

GEMV_TEST_OBJS  = $(TEST_GEMV_PATH)/gemv-support.o
SYMV_TEST_OBJS  = $(TEST_SYMV_PATH)/symv-support.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o \
		$(TEST_DOT_PATH)/print_vector.o \
		$(TEST_DOT_PATH)/copy_vector.o

TEST_SRCS = do_test_symm.c BLAS_symm_testgen.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(GEMV_TEST_OBJS) $(SYMV_TEST_OBJS) $(COMMON_OBJS)

test: do_test_symm
	@echo Testing SYMM - default quick test
	./do_test_symm 7 1 1.0 3 0.001 > symm.results

test-few: do_test_symm
	@echo Testing SYMM - preforming a small set of tests
	./do_test_symm 7 1 1.0 3 $(DO_FEW_TESTS) > symm.results

test-some: do_test_symm
	@echo Testing SYMM - preforming some tests : not all tests 
	./do_test_symm 7 1 1.0 3 $(DO_SOME_TESTS) > symm.results

test-all: do_test_symm
	@echo Testing SYMM Very Thoroughly
	./do_test_symm 7 1 1.0 3 $(DO_ALL_TESTS) > symm.results

test-none: do_test_symm
	@echo NOT Testing SYMM - PREFORMING NO TESTS. NOT TESTING
	./do_test_symm 7 1 1.0 3 $(DO_NO_TESTS) > symm.results

do_test_symm: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_symm $(EXTRA_LIBS)

test-lib: $(TEST_OBJS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK symm.results do_test_symm core

