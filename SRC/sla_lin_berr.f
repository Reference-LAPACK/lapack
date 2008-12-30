      SUBROUTINE SLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR )
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
      INTEGER            N, NZ, NRHS
*     ..
*     .. Array Arguments ..
      REAL               AYB( N, NRHS ), BERR( NRHS )
      REAL               RES( N, NRHS )
*     ..
*
*  Purpose
*  =======
*
*     SLA_LIN_BERR computes component-wise relative backward error from
*     the formula
*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
*     where abs(Z) is the component-wise absolute value of the matrix
*     or vector Z.
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               TMP
      INTEGER            I, J
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. External Functions ..
      EXTERNAL           SLAMCH
      REAL               SLAMCH
      REAL               SAFE1
*     ..
*     .. Executable Statements ..
*
*     Adding SAFE1 to the numerator guards against spuriously zero
*     residuals.  A similar safeguard is in the SLA_yyAMV routine used
*     to compute AYB.
*
      SAFE1 = SLAMCH( 'Safe minimum' )
      SAFE1 = (NZ+1)*SAFE1

      DO J = 1, NRHS
         BERR(J) = 0.0
         DO I = 1, N
            IF (AYB(I,J) .NE. 0.0) THEN
               TMP = (SAFE1+ABS(RES(I,J)))/AYB(I,J)
               BERR(J) = MAX( BERR(J), TMP )
            END IF
*
*     If AYB is exactly 0.0 (and if computed by SLA_yyAMV), then we know
*     the true residual also must be exactly 0.0.
*
         END DO
      END DO
      END SUBROUTINE
