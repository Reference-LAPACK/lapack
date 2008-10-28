include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..
TEST_DOT_PATH = ../test-dot

COMMON_OBJS = ../common/dummy_main.o

TEST_SRCS = do_test_axpby.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)

TEST_DOT_SRCS = $(TEST_DOT_PATH)/test_dot.c \
           $(TEST_DOT_PATH)/testgen_BLAS_sdot.c \
           $(TEST_DOT_PATH)/testgen_BLAS_ddot.c \
           $(TEST_DOT_PATH)/testgen_BLAS_cdot.c \
           $(TEST_DOT_PATH)/testgen_BLAS_zdot.c \
	   $(TEST_DOT_PATH)/testgen_aux.c \
           $(TEST_DOT_PATH)/BLAS_dot_testgen.c \
	   $(TEST_DOT_PATH)/print_vector.o \
	   $(TEST_DOT_PATH)/copy_vector.o

TEST_DOT_OBJS = $(TEST_DOT_SRCS:.c=.o)

ALL_OBJS = $(TEST_DOT_OBJS) $(TEST_OBJS) $(COMMON_OBJS)

test: do_test_axpby
	@echo Testing AXPBY
	./do_test_axpby 10 1 1.0 1 1 > axpby.results

do_test_axpby: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) -o do_test_axpby $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK axpby.results do_test_axpby core

