include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
TEST_DOT2_PATH = ../test-dot2
TEST_GEMV_PATH = ../test-gemv
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/print_vector.o \
		$(TEST_DOT_PATH)/copy_vector.o
TEST_DOT2_OBJS = $(TEST_DOT2_PATH)/dot2.o \
		 $(TEST_DOT2_PATH)/testgen_BLAS_sdot2.o \
		 $(TEST_DOT2_PATH)/testgen_BLAS_ddot2.o \
		 $(TEST_DOT2_PATH)/testgen_BLAS_cdot2.o \
		 $(TEST_DOT2_PATH)/testgen_BLAS_zdot2.o \
		 $(TEST_DOT2_PATH)/test_dot2.o \
		 $(TEST_DOT2_PATH)/r_truth2.o \
		 $(TEST_DOT2_PATH)/BLAS_dot2_testgen.o
TEST_GEMV_OBJS = $(TEST_GEMV_PATH)/gemv-support.o

TEST_SRCS = do_test_gemv2.c BLAS_gemv2_testgen.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(TEST_DOT2_OBJS) \
	   $(TEST_GEMV_OBJS) $(COMMON_OBJS)

test: do_test_gemv2
	@echo Testing GEMV2 - default quick test
	./do_test_gemv2 12 1 1.0 0 0.01 > gemv2.results

do_test_gemv2: $(LIB_PATH)/$(XBLASLIB) $(ALL_OBJS)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_gemv2 $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK gemv2.results do_test_gemv2 core

