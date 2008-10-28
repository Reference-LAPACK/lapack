include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

SYMV_SRCS =\
	$(SRC_PATH)/BLAS_dsymv_d_s.c\
	$(SRC_PATH)/BLAS_dsymv_s_d.c\
	$(SRC_PATH)/BLAS_dsymv_s_s.c\
	$(SRC_PATH)/BLAS_zsymv_z_c.c\
	$(SRC_PATH)/BLAS_zsymv_c_z.c\
	$(SRC_PATH)/BLAS_zsymv_c_c.c\
	$(SRC_PATH)/BLAS_csymv_c_s.c\
	$(SRC_PATH)/BLAS_csymv_s_c.c\
	$(SRC_PATH)/BLAS_csymv_s_s.c\
	$(SRC_PATH)/BLAS_zsymv_z_d.c\
	$(SRC_PATH)/BLAS_zsymv_d_z.c\
	$(SRC_PATH)/BLAS_zsymv_d_d.c\
	$(SRC_PATH)/BLAS_ssymv_x.c\
	$(SRC_PATH)/BLAS_dsymv_x.c\
	$(SRC_PATH)/BLAS_csymv_x.c\
	$(SRC_PATH)/BLAS_zsymv_x.c\
	$(SRC_PATH)/BLAS_dsymv_d_s_x.c\
	$(SRC_PATH)/BLAS_dsymv_s_d_x.c\
	$(SRC_PATH)/BLAS_dsymv_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymv_z_c_x.c\
	$(SRC_PATH)/BLAS_zsymv_c_z_x.c\
	$(SRC_PATH)/BLAS_zsymv_c_c_x.c\
	$(SRC_PATH)/BLAS_csymv_c_s_x.c\
	$(SRC_PATH)/BLAS_csymv_s_c_x.c\
	$(SRC_PATH)/BLAS_csymv_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymv_z_d_x.c\
	$(SRC_PATH)/BLAS_zsymv_d_z_x.c\
	$(SRC_PATH)/BLAS_zsymv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_dsymv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_csymv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_csymv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_csymv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_ssymv_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_d_d_x-f2c.c 

SYMV_OBJS = $(SYMV_SRCS:.c=.o)

all: $(SYMV_OBJS)

lib: $(SYMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(SYMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

