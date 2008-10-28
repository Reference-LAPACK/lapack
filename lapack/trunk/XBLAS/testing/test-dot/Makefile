include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

TEST_SRCS = do_test_dot.c test_dot.c BLAS_dot_testgen.c testgen_aux.c \
	    testgen_BLAS_sdot.c testgen_BLAS_ddot.c \
	    testgen_BLAS_cdot.c testgen_BLAS_zdot.c \
	    print_vector.c copy_vector.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

test: do_test_dot
	@echo Testing DOT
	./do_test_dot 20 1 1.0 0 1 > dot.results

test-lib: $(TEST_OBJS)

do_test_dot: $(TEST_OBJS) $(COMMON_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(TEST_OBJS) $(COMMON_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_dot $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK dot.results do_test_dot core

