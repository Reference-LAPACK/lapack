include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
TEST_GEMV_PATH = ../test-gemv
TEST_HEMV_PATH = ../test-hemv
TEST_SYMM_PATH = ../test-symm
TEST_SYMV_PATH = ../test-symv
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

GEMV_TEST_OBJS  = $(TEST_GEMV_PATH)/gemv-support.o
HEMV_TEST_OBJS  = $(TEST_HEMV_PATH)/hemv-support.o
SYMM_TEST_OBJS  = $(TEST_SYMM_PATH)/BLAS_symm_testgen.o
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

TEST_SRCS = do_test_hemm.c BLAS_hemm_testgen.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(GEMV_TEST_OBJS) $(HEMV_TEST_OBJS) \
	   $(SYMM_TEST_OBJS) $(SYMV_TEST_OBJS) $(COMMON_OBJS)

test: do_test_hemm
	@echo Testing HEMM - default quick test
	./do_test_hemm 7 1 1.0 0 0.01 > hemm.results

test-few: do_test_hemm
	@echo Testing HEMM - preforming a small set of tests
	./do_test_hemm 7 1 1.0 0 $(DO_FEW_TESTS) > hemm.results

test-some: do_test_hemm
	@echo Testing HEMM - preforming some tests : not all tests 
	./do_test_hemm 7 1 1.0 0 $(DO_SOME_TESTS) > hemm.results

test-all: do_test_hemm
	@echo Testing HEMM Very Thoroughly
	./do_test_hemm 7 1 1.0 0 $(DO_ALL_TESTS) > hemm.results

test-none: do_test_hemm
	@echo NOT Testing HEMM - PREFORMING NO TESTS. NOT TESTING
	./do_test_hemm 7 1 1.0 0 $(DO_NO_TESTS) > hemm.results

do_test_hemm: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_hemm $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK hemm.results do_test_hemm core

