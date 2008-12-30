      REAL FUNCTION CLA_GBRCOND_C( TRANS, N, KL, KU, AB, LDAB, AFB,
     $                             LDAFB, IPIV, C, CAPPLY, INFO, WORK,
     $                             RWORK )
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
      LOGICAL            CAPPLY
      INTEGER            N, KL, KU, KD, LDAB, LDAFB, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), WORK( * )
      REAL               C( * ), RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*     CLA_GBRCOND_C Computes the infinity norm condition number of
*     op(A) * inv(diag(C)) where C is a REAL vector.
*
*  Arguments
*  =========
*
*  C       REAL vector.
*
*  WORK    COMPLEX workspace of size 2*N.
*
*  RWORK   REAL workspace of size 3*N.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            NOTRANS
      INTEGER            KASE, I, J
      REAL               AINVNM, ANORM, TMP
      COMPLEX            ZDUM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACN2, CGBTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
      CLA_GBRCOND_C = 0.0E+0
*
      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      IF ( .NOT. NOTRANS .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT.
     $     LSAME( TRANS, 'C' ) ) THEN
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLA_GBRCOND_C', -INFO )
         RETURN
      END IF
*
*     Compute norm of op(A)*op2(C).
*
      ANORM = 0.0E+0
      KD = KU + 1
      IF ( NOTRANS ) THEN
         DO I = 1, N
            TMP = 0.0E+0
            IF ( CAPPLY ) THEN
               DO J = 1, N
                  IF ( I.GE.MAX( 1, J-KU )
     $                 .AND. I.LE.MIN( N, J+KL ) ) THEN
                     TMP = TMP + CABS1(AB( KD+I-J, J ) ) / C( J )
                  END IF
               END DO
            ELSE
               DO J = 1, N
                  IF ( I.GE.MAX( 1, J-KU )
     $                 .AND. I.LE.MIN( N, J+KL ) ) THEN
                     TMP = TMP + CABS1( AB( KD+I-J, J ) )
                  END IF
               END DO
            END IF
            RWORK( 2*N+I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0E+0
            IF ( CAPPLY ) THEN
               DO J = 1, N
                  IF ( I.GE.MAX( 1, J-KU )
     $                 .AND. I.LE.MIN( N, J+KL ) ) THEN
                     TMP = TMP + CABS1( AB( J, KD+I-J ) ) / C( J )
                  END IF
               END DO
            ELSE
               DO J = 1, N
                  IF ( I.GE.MAX( 1, J-KU )
     $                 .AND. I.LE.MIN( N, J+KL ) ) THEN
                     TMP = TMP + CABS1( AB( J, KD+I-J ) )
                  END IF
               END DO
            END IF
            RWORK( 2*N+I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         CLA_GBRCOND_C = 1.0E+0
         RETURN
      ELSE IF( ANORM .EQ. 0.0E+0 ) THEN
         RETURN
      END IF
*
*     Estimate the norm of inv(op(A)).
*
      AINVNM = 0.0E+0
*
      KASE = 0
   10 CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.2 ) THEN
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( 2*N+I )
            END DO
*
            IF ( NOTRANS ) THEN
               CALL CGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            ELSE
               CALL CGBTRS( 'Conjugate transpose', N, KL, KU, 1, AFB,
     $              LDAFB, IPIV, WORK, N, INFO )
            ENDIF
*
*           Multiply by inv(C).
*
            IF ( CAPPLY ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF
         ELSE
*
*           Multiply by inv(C').
*
            IF ( CAPPLY ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF
*
            IF ( NOTRANS ) THEN
               CALL CGBTRS( 'Conjugate transpose', N, KL, KU, 1, AFB,
     $              LDAFB, IPIV,  WORK, N, INFO )
            ELSE
               CALL CGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            END IF
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( 2*N+I )
            END DO
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM .NE. 0.0E+0 )
     $   CLA_GBRCOND_C = 1.0E+0 / AINVNM
*
      RETURN
*
      END
