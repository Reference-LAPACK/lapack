include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
TEST_DOT2_PATH = ../test-dot2
TEST_HEMV_PATH = ../test-hemv
TEST_SYMV_PATH = ../test-symv
TEST_SYMV2_PATH = ../test-symv2
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

TEST_DOT_OBJS =   $(TEST_DOT_PATH)/testgen_aux.o \
                  $(TEST_DOT_PATH)/print_vector.o \
                  $(TEST_DOT_PATH)/copy_vector.o
TEST_DOT2_OBJS =  $(TEST_DOT2_PATH)/dot2.o \
                  $(TEST_DOT2_PATH)/testgen_BLAS_sdot2.o \
                  $(TEST_DOT2_PATH)/testgen_BLAS_ddot2.o \
                  $(TEST_DOT2_PATH)/testgen_BLAS_cdot2.o \
                  $(TEST_DOT2_PATH)/testgen_BLAS_zdot2.o \
                  $(TEST_DOT2_PATH)/test_dot2.o \
                  $(TEST_DOT2_PATH)/r_truth2.o \
                  $(TEST_DOT2_PATH)/BLAS_dot2_testgen.o
TEST_HEMV_OBJS =  $(TEST_HEMV_PATH)/hemv-support.o
TEST_SYMV_OBJS =  $(TEST_SYMV_PATH)/symv-support.o
TEST_SYMV2_OBJS = $(TEST_SYMV2_PATH)/BLAS_symv2_testgen.o

TEST_SRCS = do_test_hemv2.c BLAS_hemv2_testgen.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(TEST_DOT2_OBJS) \
           $(TEST_HEMV_OBJS) $(TEST_SYMV_OBJS) $(TEST_SYMV2_OBJS) \
           $(COMMON_OBJS)

test: do_test_hemv2
	@echo Testing HEMV2 - default test
	./do_test_hemv2 7 1 1.0 0 0.1 > hemv2.results

do_test_hemv2: $(LIB_PATH)/$(XBLASLIB) $(ALL_OBJS)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_hemv2 $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK hemv2.results do_test_hemv2 core

