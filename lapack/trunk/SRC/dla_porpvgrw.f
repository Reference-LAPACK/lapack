      DOUBLE PRECISION FUNCTION DLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF, 
     $     LDAF, WORK )
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
      CHARACTER*1        UPLO
      INTEGER            NCOLS, LDA, LDAF
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), WORK( * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   AMAX, UMAX, RPVGRW
      LOGICAL            UPPER
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Functions ..
      EXTERNAL           LSAME, DLASET
      LOGICAL            LSAME
*     ..
*     .. Executable Statements ..
*
      UPPER = LSAME( 'Upper', UPLO )
*
*     DPOTRF will have factored only the NCOLSxNCOLS leading minor, so
*     we restrict the growth search to that minor and use only the first
*     2*NCOLS workspace entries.
*
      RPVGRW = 1.0D+0
      DO I = 1, 2*NCOLS
         WORK( I ) = 0.0D+0
      END DO
*
*     Find the max magnitude entry of each column.
*
      IF ( UPPER ) THEN
         DO J = 1, NCOLS
            DO I = 1, J
               WORK( NCOLS+J ) =
     $              MAX( ABS( A( I, J ) ), WORK( NCOLS+J ) )
            END DO
         END DO
      ELSE
         DO J = 1, NCOLS
            DO I = J, NCOLS
               WORK( NCOLS+J ) =
     $              MAX( ABS( A( I, J ) ), WORK( NCOLS+J ) )
            END DO
         END DO
      END IF
*
*     Now find the max magnitude entry of each column of the factor in
*     AF.  No pivoting, so no permutations.
*
      IF ( LSAME( 'Upper', UPLO ) ) THEN
         DO J = 1, NCOLS
            DO I = 1, J
               WORK( J ) = MAX( ABS( AF( I, J ) ), WORK( J ) )
            END DO
         END DO
      ELSE
         DO J = 1, NCOLS
            DO I = J, NCOLS
               WORK( J ) = MAX( ABS( AF( I, J ) ), WORK( J ) )
            END DO
         END DO
      END IF
*
*     Compute the *inverse* of the max element growth factor.  Dividing
*     by zero would imply the largest entry of the factor's column is
*     zero.  Than can happen when either the column of A is zero or
*     massive pivots made the factor underflow to zero.  Neither counts
*     as growth in itself, so simply ignore terms with zero
*     denominators.
*
      IF ( LSAME( 'Upper', UPLO ) ) THEN
         DO I = 1, NCOLS
            UMAX = WORK( I )
            AMAX = WORK( NCOLS+I )
            IF ( UMAX /= 0.0D+0 ) THEN
               RPVGRW = MIN( AMAX / UMAX, RPVGRW )
            END IF
         END DO
      ELSE
         DO I = 1, NCOLS
            UMAX = WORK( I )
            AMAX = WORK( NCOLS+I )
            IF ( UMAX /= 0.0D+0 ) THEN
               RPVGRW = MIN( AMAX / UMAX, RPVGRW )
            END IF
         END DO
      END IF

      DLA_PORPVGRW = RPVGRW
      END FUNCTION
