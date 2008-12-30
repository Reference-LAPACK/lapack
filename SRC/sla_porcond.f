      REAL FUNCTION SLA_PORCOND( UPLO, N, A, LDA, AF, LDAF, CMODE, C,
     $                           INFO, WORK, IWORK )
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
      CHARACTER          UPLO
      INTEGER            N, LDA, LDAF, INFO, CMODE
      REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ),
     $                   C( * )
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
*     ..
*
*  Purpose
*  =======
*
*     SLA_PORCOND Estimates the Skeel condition number of  op(A) * op2(C)
*     where op2 is determined by CMODE as follows
*     CMODE =  1    op2(C) = C
*     CMODE =  0    op2(C) = I
*     CMODE = -1    op2(C) = inv(C)
*     The Skeel condition number  cond(A) = norminf( |inv(A)||A| )
*     is computed by computing scaling factors R such that
*     diag(R)*A*op2(C) is row equilibrated and computing the standard
*     infinity-norm condition number.
*
*  Arguments
*  ==========
*
*     WORK    real workspace of size 3*N, and
*
*     IWORK   integer workspace of size N.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            KASE, I, J
      REAL               AINVNM, TMP
      LOGICAL            UP
*     ..
*     .. Array Arguments ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ISAMAX
      EXTERNAL           LSAME, ISAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLACN2, SPOTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      SLA_PORCOND = 0.0
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLA_PORCOND', -INFO )
         RETURN
      END IF

      IF( N.EQ.0 ) THEN
         SLA_PORCOND = 1.0
         RETURN
      END IF
      UP = .FALSE.
      IF ( LSAME( UPLO, 'U' ) ) UP = .TRUE.
*
*     Compute the equilibration matrix R such that
*     inv(R)*A*C has unit 1-norm.
*
      IF ( UP ) THEN
         DO I = 1, N
            TMP = 0.0
            IF ( CMODE .EQ. 1 ) THEN
               DO J = 1, I
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               END DO
            ELSE IF ( CMODE .EQ. 0 ) THEN
               DO J = 1, I
                  TMP = TMP + ABS( A( J, I ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( I, J ) )
               END DO
            ELSE
               DO J = 1, I
                  TMP = TMP + ABS( A( J ,I ) / C( J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               END DO
            END IF
            WORK( 2*N+I ) = TMP
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0
            IF ( CMODE .EQ. 1 ) THEN
               DO J = 1, I
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               END DO
            ELSE IF ( CMODE .EQ. 0 ) THEN
               DO J = 1, I
                  TMP = TMP + ABS( A( I, J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( J, I ) )
               END DO
            ELSE
               DO J = 1, I
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( J, I ) / C( J ) )
               END DO
            END IF
            WORK( 2*N+I ) = TMP
         END DO
      ENDIF
*
*     Estimate the norm of inv(op(A)).
*
      AINVNM = 0.0

      KASE = 0
   10 CONTINUE
      CALL SLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.2 ) THEN
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO

            IF (UP) THEN
               CALL SPOTRS( 'Upper', N, 1, AF, LDAF, WORK, N, INFO )
            ELSE
               CALL SPOTRS( 'Lower', N, 1, AF, LDAF, WORK, N, INFO )
            ENDIF
*
*           Multiply by inv(C).
*
            IF ( CMODE .EQ. 1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            ELSE IF ( CMODE .EQ. -1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF
         ELSE
*
*           Multiply by inv(C').
*
            IF ( CMODE .EQ. 1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            ELSE IF ( CMODE .EQ. -1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF

            IF ( UP ) THEN
               CALL SPOTRS( 'Upper', N, 1, AF, LDAF, WORK, N, INFO )
            ELSE
               CALL SPOTRS( 'Lower', N, 1, AF, LDAF, WORK, N, INFO )
            ENDIF
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM .NE. 0.0 )
     $   SLA_PORCOND = ( 1.0 / AINVNM )
*
      RETURN
*
      END
