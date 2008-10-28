include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

HPMV_SRCS =\
	$(SRC_PATH)/BLAS_zhpmv_z_c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_z.c\
	$(SRC_PATH)/BLAS_zhpmv_c_c.c\
	$(SRC_PATH)/BLAS_chpmv_c_s.c\
	$(SRC_PATH)/BLAS_zhpmv_z_d.c\
	$(SRC_PATH)/BLAS_chpmv_x.c\
	$(SRC_PATH)/BLAS_zhpmv_x.c\
	$(SRC_PATH)/BLAS_zhpmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zhpmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zhpmv_c_c_x.c\
	$(SRC_PATH)/BLAS_chpmv_c_s_x.c\
	$(SRC_PATH)/BLAS_zhpmv_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_zhpmv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chpmv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chpmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chpmv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_z_d_x-f2c.c 

HPMV_OBJS = $(HPMV_SRCS:.c=.o)

all: $(HPMV_OBJS)

lib: $(HPMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(HPMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

