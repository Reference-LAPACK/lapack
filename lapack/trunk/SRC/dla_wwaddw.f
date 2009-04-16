      SUBROUTINE DLA_WWADDW( N, X, Y, W )
*
*     -- LAPACK routine (version 3.2.1)                                 --
*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and --
*     -- Jason Riedy of Univ. of California Berkeley.                 --
*     -- April 2009                                                   --
*
*     -- LAPACK is a software package provided by Univ. of Tennessee, --
*     -- Univ. of California Berkeley and NAG Ltd.                    --
*
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * ), W( * )
*     ..
*
*     Purpose
*     =======
*
*     DLA_WWADDW adds a vector W into a doubled-single vector (X, Y).
*
*     This works for all extant IBM's hex and binary floating point
*     arithmetics, but not for decimal.
*
*     Arguments
*     =========
*
*     N      (input) INTEGER
*            The length of vectors X, Y, and W.
*
*     X, Y   (input/output) DOUBLE PRECISION array, length N
*            The doubled-single accumulation vector.
*
*     W      (input) DOUBLE PRECISION array, length N
*            The vector to be added.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   S
      INTEGER            I
*     ..
*     .. Executable Statements ..
*
      DO 10 I = 1, N
        S = X(I) + W(I)
        S = (S + S) - S
        Y(I) = ((X(I) - S) + W(I)) + Y(I)
        X(I) = S
 10   CONTINUE
      RETURN
      END
