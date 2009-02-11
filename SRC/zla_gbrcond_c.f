      DOUBLE PRECISION FUNCTION ZLA_GBRCOND_C( TRANS, N, KL, KU, AB,
     $                             LDAB, AFB, LDAFB, IPIV, C, CAPPLY,
     $                             INFO, WORK, RWORK )
*
*     -- LAPACK routine (version 3.2.1)                               --
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
      INTEGER            N, KL, KU, KD, KE, LDAB, LDAFB, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), WORK( * )
      DOUBLE PRECISION   C( * ), RWORK( * )
*
*
*  Purpose
*  =======
*
*     ZLA_GBRCOND_C Computes the infinity norm condition number of
*     op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector.
*
*  Arguments
*  =========
*
*   C      DOUBLE PRECISION vector.
*
*   WORK   COMPLEX*16 workspace of size 2*N.
*
*   RWORK  DOUBLE PRECISION workspace of size N.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            NOTRANS
      INTEGER            KASE, I, J
      DOUBLE PRECISION   AINVNM, ANORM, TMP
      COMPLEX*16         ZDUM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLACN2, ZGBTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
      ZLA_GBRCOND_C = 0.0D+0
*
      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      IF ( .NOT. NOTRANS .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT.
     $     LSAME( TRANS, 'C' ) ) THEN
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 .OR. KL.GT.N-1 ) THEN
         INFO = -4
      ELSE IF( KU.LT.0 .OR. KU.GT.N-1 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KL+KU+1 ) THEN
         INFO = -8
      ELSE IF( LDAFB.LT.2*KL+KU+1 ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLA_GBRCOND_C', -INFO )
         RETURN
      END IF
*
*     Compute norm of op(A)*op2(C).
*
      ANORM = 0.0D+0
      KD = KU + 1
      KE = KL + 1
      IF ( NOTRANS ) THEN
         DO I = 1, N
            TMP = 0.0D+0
            IF ( CAPPLY ) THEN
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + CABS1( AB( KD+I-J, J ) ) / C( J )
               END DO
            ELSE
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + CABS1( AB( KD+I-J, J ) )
               END DO
            END IF
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0D+0
            IF ( CAPPLY ) THEN
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + CABS1( AB( KE-I+J, I ) ) / C( J )
               END DO
            ELSE
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + CABS1( AB( KE-I+J, I ) )
               END DO
            END IF
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         ZLA_GBRCOND_C = 1.0D+0
         RETURN
      ELSE IF( ANORM .EQ. 0.0D+0 ) THEN
         RETURN
      END IF
*
*     Estimate the norm of inv(op(A)).
*
      AINVNM = 0.0D+0
*
      KASE = 0
   10 CONTINUE
      CALL ZLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.2 ) THEN
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
*
            IF ( NOTRANS ) THEN
               CALL ZGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            ELSE
               CALL ZGBTRS( 'Conjugate transpose', N, KL, KU, 1, AFB,
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
               CALL ZGBTRS( 'Conjugate transpose', N, KL, KU, 1, AFB,
     $              LDAFB, IPIV,  WORK, N, INFO )
            ELSE
               CALL ZGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            END IF
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM .NE. 0.0D+0 )
     $   ZLA_GBRCOND_C = 1.0D+0 / AINVNM
*
      RETURN
*
      END
