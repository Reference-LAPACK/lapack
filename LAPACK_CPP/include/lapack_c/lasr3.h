#ifndef LAPACK_C_ROT_HPP
#define LAPACK_C_ROT_HPP

#include "lapack_c/util.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void lapack_c_slasr3(char side,
                     char direct,
                     lapack_idx m,
                     lapack_idx n,
                     lapack_idx k,
                     const float* C,
                     lapack_idx ldc,
                     const float* S,
                     lapack_idx lds,
                     float* A,
                     lapack_idx lda,
                     float* work,
                     lapack_idx lwork);

void lapack_c_dlasr3(char side,
                     char direct,
                     lapack_idx m,
                     lapack_idx n,
                     lapack_idx k,
                     const double* C,
                     lapack_idx ldc,
                     const double* S,
                     lapack_idx lds,
                     double* A,
                     lapack_idx lda,
                     double* work,
                     lapack_idx lwork);

void lapack_c_clasr3(char side,
                     char direct,
                     lapack_idx m,
                     lapack_idx n,
                     lapack_idx k,
                     const float* C,
                     lapack_idx ldc,
                     const lapack_float_complex* S,
                     lapack_idx lds,
                     lapack_float_complex* A,
                     lapack_idx lda,
                     lapack_float_complex* work,
                     lapack_idx lwork);

void lapack_c_sclasr3(char side,
                      char direct,
                      lapack_idx m,
                      lapack_idx n,
                      lapack_idx k,
                      const float* C,
                      lapack_idx ldc,
                      const float* S,
                      lapack_idx lds,
                      lapack_float_complex* A,
                      lapack_idx lda,
                      lapack_float_complex* work,
                      lapack_idx lwork);

void lapack_c_zlasr3(lapack_idx m,
                     lapack_idx n,
                     lapack_idx k,
                     const double* C,
                     lapack_idx ldc,
                     const lapack_double_complex* S,
                     lapack_idx lds,
                     lapack_double_complex* A,
                     lapack_idx lda,
                     lapack_double_complex* work,
                     lapack_idx lwork);

void lapack_c_dzlasr3(lapack_idx m,
                      lapack_idx n,
                      lapack_idx k,
                      const double* C,
                      lapack_idx ldc,
                      const double* S,
                      lapack_idx lds,
                      lapack_double_complex* A,
                      lapack_idx lda,
                      lapack_double_complex* work,
                      lapack_idx lwork);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif