      DOUBLE PRECISION FUNCTION ZLA_RPVGRW( N, NCOLS, A, LDA, AF, LDAF )
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
      INTEGER            N, NCOLS, LDA, LDAF
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), AF( LDAF, * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   AMAX, UMAX, RPVGRW
      COMPLEX*16         ZDUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, ABS, REAL, DIMAG
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      RPVGRW = 1.0D+0
*
      DO J = 1, NCOLS
         AMAX = 0.0D+0
         UMAX = 0.0D+0
         DO I = 1, N
            AMAX = MAX( CABS1( A( I, J ) ), AMAX )
         END DO
         DO I = 1, J
            UMAX = MAX( CABS1( AF( I, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0D+0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      ZLA_RPVGRW = RPVGRW
      END FUNCTION
