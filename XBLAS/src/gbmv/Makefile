include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

GBMV_SRCS =\
	$(SRC_PATH)/BLAS_dgbmv_d_s.c\
	$(SRC_PATH)/BLAS_dgbmv_s_d.c\
	$(SRC_PATH)/BLAS_dgbmv_s_s.c\
	$(SRC_PATH)/BLAS_zgbmv_z_c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_z.c\
	$(SRC_PATH)/BLAS_zgbmv_c_c.c\
	$(SRC_PATH)/BLAS_cgbmv_c_s.c\
	$(SRC_PATH)/BLAS_cgbmv_s_c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_s.c\
	$(SRC_PATH)/BLAS_zgbmv_z_d.c\
	$(SRC_PATH)/BLAS_zgbmv_d_z.c\
	$(SRC_PATH)/BLAS_zgbmv_d_d.c\
	$(SRC_PATH)/BLAS_sgbmv_x.c\
	$(SRC_PATH)/BLAS_dgbmv_x.c\
	$(SRC_PATH)/BLAS_cgbmv_x.c\
	$(SRC_PATH)/BLAS_zgbmv_x.c\
	$(SRC_PATH)/BLAS_dgbmv_d_s_x.c\
	$(SRC_PATH)/BLAS_dgbmv_s_d_x.c\
	$(SRC_PATH)/BLAS_dgbmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zgbmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zgbmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zgbmv_c_c_x.c\
	$(SRC_PATH)/BLAS_cgbmv_c_s_x.c\
	$(SRC_PATH)/BLAS_cgbmv_s_c_x.c\
	$(SRC_PATH)/BLAS_cgbmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zgbmv_z_d_x.c\
	$(SRC_PATH)/BLAS_zgbmv_d_z_x.c\
	$(SRC_PATH)/BLAS_zgbmv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_dgbmv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sgbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_d_d_x-f2c.c 

GBMV_OBJS = $(GBMV_SRCS:.c=.o)

all: $(GBMV_OBJS)

lib: $(GBMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(GBMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

