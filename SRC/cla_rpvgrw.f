      REAL FUNCTION CLA_RPVGRW( N, NCOLS, A, LDA, AF, LDAF )
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
      COMPLEX            A( LDA, * ), AF( LDAF, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               AMAX, UMAX, RPVGRW
      COMPLEX            ZDUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, ABS, REAL, AIMAG
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      RPVGRW = 1.0
*
      DO J = 1, NCOLS
         AMAX = 0.0
         UMAX = 0.0
         DO I = 1, N
            AMAX = MAX( CABS1( A( I, J ) ), AMAX )
         END DO
         DO I = 1, J
            UMAX = MAX( CABS1( AF( I, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      CLA_RPVGRW = RPVGRW
      END FUNCTION
