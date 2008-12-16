      REAL FUNCTION CLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB,
     $                            LDAFB )
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
      INTEGER            N, KL, KU, NCOLS, LDAB, LDAFB
*     ..
*     .. Array Arguments ..
      COMPLEX            AB( LDAB, * ), AFB( LDAFB, * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, KD
      REAL               AMAX, UMAX, RPVGRW
      COMPLEX            ZDUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, REAL, AIMAG
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
      KD = KU + 1
      DO J = 1, NCOLS
         AMAX = 0.0
         UMAX = 0.0
         DO I = MAX( J-KU, 1 ), MIN( J+KL, N )
            AMAX = MAX( CABS1( AB( KD+I-J, J ) ), AMAX )
         END DO
         DO I = MAX( J-KU, 1 ), J
            UMAX = MAX( CABS1( AFB( KD+I-J, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      CLA_GBRPVGRW = RPVGRW
      END FUNCTION
