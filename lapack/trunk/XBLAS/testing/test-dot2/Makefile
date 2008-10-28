include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

TEST_SRCS = do_test_dot2.c test_dot2.c BLAS_dot2_testgen.c \
	    dot2.c r_truth2.c \
	    testgen_BLAS_sdot2.c testgen_BLAS_ddot2.c \
	    testgen_BLAS_zdot2.c testgen_BLAS_cdot2.c

TEST_OBJS = $(TEST_SRCS:.c=.o)
TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/print_vector.o \
		$(TEST_DOT_PATH)/copy_vector.o

COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(COMMON_OBJS)

test: do_test_dot2
	@echo Testing DOT2
	./do_test_dot2 20 1 1.0 0 1 > dot2.results

test-lib: $(TEST_OBJS)

do_test_dot2: $(TEST_OBJS) $(COMMON_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_dot2 $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK dot2.results do_test_dot2 core

