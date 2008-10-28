include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
TEST_HEMV_PATH = ../test-hemv
TEST_SYMV_PATH = ../test-symv
TEST_SBMV_PATH = ../test-sbmv
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

HEMV_TEST_OBJS	= $(TEST_HEMV_PATH)/hemv-support.o
SYMV_TEST_OBJS  = $(TEST_SYMV_PATH)/symv-support.o \
		$(TEST_SYMV_PATH)/BLAS_symv_testgen.o
SBMV_TEST_OBJS  = $(TEST_SBMV_PATH)/sbmv-support.o \
		$(TEST_SBMV_PATH)/BLAS_sbmv_testgen.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o \
		$(TEST_DOT_PATH)/print_vector.o \
		$(TEST_DOT_PATH)/copy_vector.o

TEST_SRCS = do_test_hbmv.c BLAS_hbmv_testgen.c hbmv-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(HEMV_TEST_OBJS) \
		$(SYMV_TEST_OBJS) $(SBMV_TEST_OBJS) $(COMMON_OBJS)

test: do_test_hbmv
	@echo Testing HBMV - default quick test
	./do_test_hbmv 7 1 1.0 3 0.01 > hbmv.results

test-few: do_test_hbmv
	@echo Testing HBMV - preforming a small set of tests
	./do_test_hbmv 7 1 1.0 3 $(DO_FEW_TESTS) > hbmv.results

test-some: do_test_hbmv
	@echo Testing HBMV - preforming some tests : not all tests 
	./do_test_hbmv 7 1 1.0 3 $(DO_SOME_TESTS) > hbmv.results

test-all: do_test_hbmv
	@echo Testing HBMV Very Thoroughly
	./do_test_hbmv 7 1 1.0 3 $(DO_ALL_TESTS) > hbmv.results

test-none: do_test_hbmv
	@echo NOT Testing HBMV - PREFORMING NO TESTS. NOT TESTING
	./do_test_hbmv 7 1 1.0 3 $(DO_NO_TESTS) > hbmv.results

do_test_hbmv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_hbmv $(EXTRA_LIBS)

test-lib: $(TEST_OBJS)

.c.o: 
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK hbmv.results do_test_hbmv core

