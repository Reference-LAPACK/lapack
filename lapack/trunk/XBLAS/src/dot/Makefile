include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

DOT_SRCS =\
	$(SRC_PATH)/BLAS_ddot_d_s.c\
	$(SRC_PATH)/BLAS_ddot_s_d.c\
	$(SRC_PATH)/BLAS_ddot_s_s.c\
	$(SRC_PATH)/BLAS_zdot_z_c.c\
	$(SRC_PATH)/BLAS_zdot_c_z.c\
	$(SRC_PATH)/BLAS_zdot_c_c.c\
	$(SRC_PATH)/BLAS_cdot_c_s.c\
	$(SRC_PATH)/BLAS_cdot_s_c.c\
	$(SRC_PATH)/BLAS_cdot_s_s.c\
	$(SRC_PATH)/BLAS_zdot_z_d.c\
	$(SRC_PATH)/BLAS_zdot_d_z.c\
	$(SRC_PATH)/BLAS_zdot_d_d.c\
	$(SRC_PATH)/BLAS_sdot_x.c\
	$(SRC_PATH)/BLAS_ddot_x.c\
	$(SRC_PATH)/BLAS_cdot_x.c\
	$(SRC_PATH)/BLAS_zdot_x.c\
	$(SRC_PATH)/BLAS_ddot_d_s_x.c\
	$(SRC_PATH)/BLAS_ddot_s_d_x.c\
	$(SRC_PATH)/BLAS_ddot_s_s_x.c\
	$(SRC_PATH)/BLAS_zdot_z_c_x.c\
	$(SRC_PATH)/BLAS_zdot_c_z_x.c\
	$(SRC_PATH)/BLAS_zdot_c_c_x.c\
	$(SRC_PATH)/BLAS_cdot_c_s_x.c\
	$(SRC_PATH)/BLAS_cdot_s_c_x.c\
	$(SRC_PATH)/BLAS_cdot_s_s_x.c\
	$(SRC_PATH)/BLAS_zdot_z_d_x.c\
	$(SRC_PATH)/BLAS_zdot_d_z_x.c\
	$(SRC_PATH)/BLAS_zdot_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_ddot_d_s-f2c.c\
	$(SRC_PATH)/BLAS_ddot_s_d-f2c.c\
	$(SRC_PATH)/BLAS_ddot_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zdot_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zdot_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zdot_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cdot_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cdot_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cdot_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zdot_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zdot_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zdot_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sdot_x-f2c.c\
	$(SRC_PATH)/BLAS_ddot_x-f2c.c\
	$(SRC_PATH)/BLAS_cdot_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_x-f2c.c\
	$(SRC_PATH)/BLAS_ddot_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ddot_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_ddot_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cdot_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cdot_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cdot_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_d_d_x-f2c.c 

DOT_OBJS = $(DOT_SRCS:.c=.o)

all: $(DOT_OBJS)

lib: $(DOT_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(DOT_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

