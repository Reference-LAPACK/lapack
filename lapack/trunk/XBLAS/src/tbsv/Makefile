include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

TBSV_SRCS =\
	$(SRC_PATH)/BLAS_dtbsv_s.c\
	$(SRC_PATH)/BLAS_ztbsv_c.c\
	$(SRC_PATH)/BLAS_ctbsv_s.c\
	$(SRC_PATH)/BLAS_ztbsv_d.c\
	$(SRC_PATH)/BLAS_stbsv_x.c\
	$(SRC_PATH)/BLAS_dtbsv_x.c\
	$(SRC_PATH)/BLAS_dtbsv_s_x.c\
	$(SRC_PATH)/BLAS_ctbsv_x.c\
	$(SRC_PATH)/BLAS_ztbsv_x.c\
	$(SRC_PATH)/BLAS_ztbsv_c_x.c\
	$(SRC_PATH)/BLAS_ctbsv_s_x.c\
	$(SRC_PATH)/BLAS_ztbsv_d_x.c \
	\
	$(SRC_PATH)/BLAS_dtbsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_c-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_d-f2c.c\
	$(SRC_PATH)/BLAS_stbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_c_x-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_d_x-f2c.c 

TBSV_OBJS = $(TBSV_SRCS:.c=.o)

all: $(TBSV_OBJS)

lib: $(TBSV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(TBSV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

