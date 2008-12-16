      REAL FUNCTION SLA_RPVGRW( N, NCOLS, A, LDA, AF, LDAF )
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
      REAL               A( LDA, * ), AF( LDAF, * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               AMAX, UMAX, RPVGRW
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      RPVGRW = 1.0
*
      DO J = 1, NCOLS
         AMAX = 0.0
         UMAX = 0.0
         DO I = 1, N
            AMAX = MAX( ABS( A( I, J ) ), AMAX )
         END DO
         DO I = 1, J
            UMAX = MAX( ABS( AF( I, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      SLA_RPVGRW = RPVGRW
      END FUNCTION
