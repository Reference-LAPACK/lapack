include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

SPMV_SRCS =\
	$(SRC_PATH)/BLAS_dspmv_d_s.c\
	$(SRC_PATH)/BLAS_dspmv_s_d.c\
	$(SRC_PATH)/BLAS_dspmv_s_s.c\
	$(SRC_PATH)/BLAS_zspmv_z_c.c\
	$(SRC_PATH)/BLAS_zspmv_c_z.c\
	$(SRC_PATH)/BLAS_zspmv_c_c.c\
	$(SRC_PATH)/BLAS_cspmv_c_s.c\
	$(SRC_PATH)/BLAS_cspmv_s_c.c\
	$(SRC_PATH)/BLAS_cspmv_s_s.c\
	$(SRC_PATH)/BLAS_zspmv_z_d.c\
	$(SRC_PATH)/BLAS_zspmv_d_z.c\
	$(SRC_PATH)/BLAS_zspmv_d_d.c\
	$(SRC_PATH)/BLAS_sspmv_x.c\
	$(SRC_PATH)/BLAS_dspmv_x.c\
	$(SRC_PATH)/BLAS_cspmv_x.c\
	$(SRC_PATH)/BLAS_zspmv_x.c\
	$(SRC_PATH)/BLAS_dspmv_d_s_x.c\
	$(SRC_PATH)/BLAS_dspmv_s_d_x.c\
	$(SRC_PATH)/BLAS_dspmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zspmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zspmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zspmv_c_c_x.c\
	$(SRC_PATH)/BLAS_cspmv_c_s_x.c\
	$(SRC_PATH)/BLAS_cspmv_s_c_x.c\
	$(SRC_PATH)/BLAS_cspmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zspmv_z_d_x.c\
	$(SRC_PATH)/BLAS_zspmv_d_z_x.c\
	$(SRC_PATH)/BLAS_zspmv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_dspmv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sspmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_x-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_d_d_x-f2c.c 

SPMV_OBJS = $(SPMV_SRCS:.c=.o)

all: $(SPMV_OBJS)

lib: $(SPMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(SPMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

