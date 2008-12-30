      REAL FUNCTION SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB,
     $                           IPIV, CMODE, C, INFO, WORK, IWORK )
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
      CHARACTER          TRANS
      INTEGER            N, LDAB, LDAFB, INFO, KL, KU, CMODE
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), IPIV( * )
      REAL               AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ),
     $                   C( * )
*    ..
*
*  Purpose
*  =======
*
*     SLA_GERCOND Estimates the Skeel condition number of  op(A) * op2(C)
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
*  WORK     real workspace of size 5*N.
*
*  IWORK    integer workspace of size N.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            NOTRANS
      INTEGER            KASE, I, J, KD
      REAL               AINVNM, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLACN2, SGBTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      SLA_GBRCOND = 0.0
*
      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      IF ( .NOT. NOTRANS .AND. .NOT. LSAME(TRANS, 'T')
     $     .AND. .NOT. LSAME(TRANS, 'C') ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -4
      ELSE IF( KU.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KL+KU+1 ) THEN
         INFO = -8
      ELSE IF( LDAFB.LT.2*KL+KU+1 ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLA_GBRCOND', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 ) THEN
         SLA_GBRCOND = 1.0
         RETURN
      END IF
*
*     Compute the equilibration matrix R such that
*     inv(R)*A*C has unit 1-norm.
*
      KD = KU + 1
      IF ( NOTRANS ) THEN
         DO I = 1, N
            TMP = 0.0
               IF ( CMODE .EQ. 1 ) THEN
                  DO J = 1, N
                     IF ( I.GE.MAX( 1, J-KU )
     $                    .AND. I.LE.MIN( N, J+KL ) ) THEN
                        TMP = TMP + ABS( AB( KD+I-J, J ) * C( J ) )
                     END IF
                  END DO
               ELSE IF ( CMODE .EQ. 0 ) THEN
                  DO J = 1, N
                     IF ( I.GE.MAX( 1, J-KU )
     $                    .AND. I.LE.MIN( N, J+KL ) ) THEN
                        TMP = TMP + ABS( AB( KD+I-J, J ) )
                     END IF
                  END DO
               ELSE
                  DO J = 1, N
                     IF ( I.GE.MAX( 1, J-KU )
     $                    .AND. I.LE.MIN( N, J+KL ) ) THEN
                        TMP = TMP + ABS( AB( KD+I-J, J ) / C( J ) )
                     END IF
                  END DO
               END IF
            WORK( 2*N+I ) = TMP
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0
            IF ( CMODE .EQ. 1 ) THEN
               DO J = 1, N
                  IF ( I.GE.MAX( 1, J-KU )
     $                 .AND. I.LE.MIN( N, J+KL ) ) THEN
                     TMP = TMP + ABS( AB( J, KD+I-J ) * C( J ) )
                  END IF
               END DO
            ELSE IF ( CMODE .EQ. 0 ) THEN
               DO J = 1, N
                  IF ( I.GE.MAX( 1, J-KU )
     $                 .AND. I.LE.MIN( N, J+KL ) ) THEN
                     TMP = TMP + ABS(AB(J,KD+I-J))
                  END IF
               END DO
            ELSE
               DO J = 1, N
                  IF ( I.GE.MAX( 1, J-KU )
     $                 .AND. I.LE.MIN( N, J+KL ) ) THEN
                     TMP = TMP + ABS( AB( J, KD+I-J ) / C( J ) )
                  END IF
               END DO
            END IF
            WORK( 2*N+I ) = TMP
         END DO
      END IF
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

            IF ( NOTRANS ) THEN
               CALL SGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            ELSE
               CALL SGBTRS( 'Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV,
     $              WORK, N, INFO )
            END IF
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

            IF ( NOTRANS ) THEN
               CALL SGBTRS( 'Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV,
     $              WORK, N, INFO )
            ELSE
               CALL SGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            END IF
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
     $   SLA_GBRCOND = ( 1.0 / AINVNM )
*
      RETURN
*
      END
