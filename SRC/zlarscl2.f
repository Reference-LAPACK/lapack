      SUBROUTINE ZLARSCL2 ( M, N, D, X, LDX )
*
*     -- LAPACK routine (version 3.2)                                 --
*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and --
*     -- Jason Riedy of Univ. of California Berkeley.                 --
*     -- November 2008                                                --
*
*     -- LAPACK is a software package provided by Univ. of Tennessee, --
*     -- Univ. of California Berkeley and NAG Ltd.                    --
*
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      INTEGER            M, N, LDX
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( LDX, * )
      DOUBLE PRECISION   D( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARSCL2 performs a reciprocal diagonal scaling on an vector:
*    x <-- inv(D) * x
*  where the diagonal matrix D is stored as a vector.
*  Eventually to be replaced by BLAS_sge_diag_scale in the new BLAS
*  standard.
*
*  Arguments
*  =========
*  N      (input) INTEGER
*         The size of the vectors X and D.
*
*  D      (input) DOUBLE PRECISION array, length N
*         Diagonal matrix D, stored as a vector of length N.
*  X      (input/output) COMPLEX*16 array, length N
*         On entry, the vector X to be scaled by D.
*         On exit, the scaled vector.
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Executable Statements ..
*
      DO J = 1, N
         DO I = 1, M
            X(I,J) = X(I,J) / D(I)
         END DO
      END DO
*
      RETURN
      END
*
