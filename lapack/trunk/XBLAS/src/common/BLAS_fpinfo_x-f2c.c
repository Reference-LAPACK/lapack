#include "blas_enum.h"
#include "f2c-bridge.h"

extern int BLAS_fpinfo_x(enum blas_cmach_type, enum blas_prec_type);

int FC_FUNC_(blas_fpinfo_x,BLAS_FPINFO_X)(int *cmach, int *prec)
{
  return BLAS_fpinfo_x((enum blas_cmach_type)*cmach,
		       (enum blas_prec_type)*prec);
}
