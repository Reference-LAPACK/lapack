include ../../make.conf
include ../../$(MAKEINC)

TEST_DOT_PATH = ../test-dot
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..
TEST_GBMV_PATH = ../test-gbmv
TEST_TRSV_PATH = ../test-trsv

TEST_GBMV_OBJS = $(TEST_GBMV_PATH)/gbmv-support.o
TEST_TRSV_OBJS = $(TEST_TRSV_PATH)/BLAS_dot_testgen_trsv.o \
	$(TEST_TRSV_PATH)/testgen_BLAS_sdot_x.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o \
		$(TEST_DOT_PATH)/print_vector.o \
		$(TEST_DOT_PATH)/copy_vector.o

TEST_SRCS = do_test_tbsv.c BLAS_tbsv_testgen.c tbsv-support.c

TEST_OBJS = $(TEST_SRCS:.c=.o) 
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) \
	$(TEST_GBMV_OBJS) $(TEST_TRSV_OBJS) $(COMMON_OBJS)

test: do_test_tbsv
	@echo Testing TBSV - default quick test
	./do_test_tbsv 9 1 1.0 3 0.01 > tbsv.results

test-few: do_test_tbsv
	@echo Testing TBSV - preforming a small set of tests
	./do_test_tbsv 9 1 1.0 3 $(DO_FEW_TESTS) > tbsv.results

test-some: do_test_tbsv
	@echo Testing TBSV - preforming some tests : not all tests 
	./do_test_tbsv 9 1 1.0 3 $(DO_SOME_TESTS) > tbsv.results

test-all: do_test_tbsv
	@echo Testing TBSV Very Thoroughly
	./do_test_tbsv 9 1 1.0 3 $(DO_ALL_TESTS) > tbsv.results

test-none: do_test_tbsv
	@echo NOT Testing TBSV - PREFORMING NO TESTS. NOT TESTING
	./do_test_tbsv 9 1 1.0 3 $(DO_NO_TESTS) > tbsv.results

do_test_tbsv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_tbsv $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK tbsv.results do_test_tbsv core

