#ifndef BLAS_EXTENDED_PRIVATE_H
#define BLAS_EXTENDED_PRIVATE_H

/* constants */
#define BITS_S  24
#define BITS_D  53
#define BITS_E  106

/* Split a double into 2 parts with at most 26 bits each. (2^27 + 1) */
#define split 	(134217729.0)

/* macros */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#include "blas_fpu.h"

#endif /* BLAS_EXTENDED_PRIVATE_H */
