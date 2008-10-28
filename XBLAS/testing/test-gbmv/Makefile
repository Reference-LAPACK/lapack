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

TEST_SRCS = do_test_gbmv.c BLAS_gbmv_testgen.c gbmv-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(COMMON_OBJS)

test: do_test_gbmv
	@echo Testing GBMV
	./do_test_gbmv 12 1 1.0 0 0.01 > gbmv.results

do_test_gbmv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_gbmv $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK gbmv.results do_test_gbmv core

