include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

TRMV_SRCS =\
	$(SRC_PATH)/BLAS_dtrmv_s.c\
	$(SRC_PATH)/BLAS_ztrmv_c.c\
	$(SRC_PATH)/BLAS_ctrmv_s.c\
	$(SRC_PATH)/BLAS_ztrmv_d.c\
	$(SRC_PATH)/BLAS_strmv_x.c\
	$(SRC_PATH)/BLAS_dtrmv_x.c\
	$(SRC_PATH)/BLAS_ctrmv_x.c\
	$(SRC_PATH)/BLAS_ztrmv_x.c\
	$(SRC_PATH)/BLAS_dtrmv_s_x.c\
	$(SRC_PATH)/BLAS_ztrmv_c_x.c\
	$(SRC_PATH)/BLAS_ctrmv_s_x.c\
	$(SRC_PATH)/BLAS_ztrmv_d_x.c \
	\
	$(SRC_PATH)/BLAS_dtrmv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_c-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_d-f2c.c\
	$(SRC_PATH)/BLAS_strmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_c_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_d_x-f2c.c 

TRMV_OBJS = $(TRMV_SRCS:.c=.o)

all: $(TRMV_OBJS)

lib: $(TRMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(TRMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

