include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

AXPBY_SRCS =\
	$(SRC_PATH)/BLAS_daxpby_s.c\
	$(SRC_PATH)/BLAS_caxpby_s.c\
	$(SRC_PATH)/BLAS_zaxpby_c.c\
	$(SRC_PATH)/BLAS_zaxpby_d.c\
	$(SRC_PATH)/BLAS_saxpby_x.c\
	$(SRC_PATH)/BLAS_daxpby_x.c\
	$(SRC_PATH)/BLAS_caxpby_x.c\
	$(SRC_PATH)/BLAS_zaxpby_x.c\
	$(SRC_PATH)/BLAS_daxpby_s_x.c\
	$(SRC_PATH)/BLAS_zaxpby_c_x.c\
	$(SRC_PATH)/BLAS_caxpby_s_x.c\
	$(SRC_PATH)/BLAS_zaxpby_d_x.c \
	\
	$(SRC_PATH)/BLAS_daxpby_s-f2c.c\
	$(SRC_PATH)/BLAS_caxpby_s-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_c-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_d-f2c.c\
	$(SRC_PATH)/BLAS_saxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_daxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_caxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_daxpby_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_c_x-f2c.c\
	$(SRC_PATH)/BLAS_caxpby_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_d_x-f2c.c 

AXPBY_OBJS = $(AXPBY_SRCS:.c=.o)

all: $(AXPBY_OBJS)

lib: $(AXPBY_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(AXPBY_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

