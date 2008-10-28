include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

WAXPBY_SRCS =\
	$(SRC_PATH)/BLAS_dwaxpby_d_s.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_d.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_s.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_z.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_c.c\
	$(SRC_PATH)/BLAS_cwaxpby_c_s.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_s.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_d.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_z.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_d.c\
	$(SRC_PATH)/BLAS_swaxpby_x.c\
	$(SRC_PATH)/BLAS_dwaxpby_x.c\
	$(SRC_PATH)/BLAS_cwaxpby_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_x.c\
	$(SRC_PATH)/BLAS_dwaxpby_d_s_x.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_d_x.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_s_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_c_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_z_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_c_x.c\
	$(SRC_PATH)/BLAS_cwaxpby_c_s_x.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_c_x.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_s_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_d_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_z_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_dwaxpby_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_d-f2c.c\
	$(SRC_PATH)/BLAS_swaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_d_x-f2c.c 

WAXPBY_OBJS = $(WAXPBY_SRCS:.c=.o)

all: $(WAXPBY_OBJS)

lib: $(WAXPBY_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(WAXPBY_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

