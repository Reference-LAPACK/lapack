include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
TEST_GEMV_PATH = ../test-gemv
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

TEST_GEMV_OBJS = $(TEST_GEMV_PATH)/gemv-support.o

TEST_SRCS = do_test_gemm.c BLAS_gemm_testgen.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(TEST_GEMV_OBJS) $(COMMON_OBJS)

test: do_test_gemm
	@echo Testing GEMM - default quick test
	./do_test_gemm 9 1 1.0 0 0.001 > gemm.results

test-lib: $(TEST_OBJS)

do_test_gemm: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_gemm $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK gemm.results do_test_gemm core

