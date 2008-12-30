      DOUBLE PRECISION FUNCTION ZLA_GBRCOND_X( TRANS, N, KL, KU, AB, 
     $                             LDAB, AFB, LDAFB, IPIV, X, INFO, 
     $     WORK, RWORK )
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
      INTEGER            N, KL, KU, KD, LDAB, LDAFB, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ),
     $                   X( * )
      DOUBLE PRECISION   RWORK( * )
*
*
*  Purpose
*  =======
*
*     ZLA_GBRCOND_X Computes the infinity norm condition number of
*     op(A) * diag(X) where X is a COMPLEX*16 vector.
*
*  Arguments
*  =========
*
*  X     COMPLEX*16 vector.
*
*  WORK  COMPLEX*16 workspace of size 2*N.
*
*  RWORK DOUBLE PRECISION workspace of size 3*N.
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
*
      ZLA_GBRCOND_X = 0.0D+0
*
      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      IF ( .NOT. NOTRANS .AND. .NOT. LSAME(TRANS, 'T') .AND. .NOT.
     $     LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLA_GBRCOND_X', -INFO )
         RETURN
      END IF
*
*     Compute norm of op(A)*op2(C).
*
      KD = KU + 1
      ANORM = 0.0D+0
      IF ( NOTRANS ) THEN
         DO I = 1, N
            TMP = 0.0D+0
            DO J = 1, N
               IF ( I.GE.MAX( 1, J-KU ) .AND. I.LE.MIN( N, J+KL ) ) THEN
                  TMP = TMP + CABS1( AB( KD+I-J, J) * X( J ) )
               END IF
            END DO
            RWORK( 2*N+I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0D+0
            DO J = 1, N
               IF ( I.GE.MAX( 1, J-KU ) .AND. I.LE.MIN( N, J+KL ) ) THEN
                  TMP = TMP + CABS1( AB( J, KD+I-J ) * X( J ) )
               END IF
            END DO
            RWORK( 2*N+I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         ZLA_GBRCOND_X = 1.0D+0
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
               WORK( I ) = WORK( I ) * RWORK( 2*N+I )
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
*           Multiply by inv(X).
*
            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO
         ELSE
*
*           Multiply by inv(X').
*
            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO
*
            IF ( NOTRANS ) THEN
               CALL ZGBTRS( 'Conjugate transpose', N, KL, KU, 1, AFB,
     $              LDAFB, IPIV, WORK, N, INFO )
            ELSE
               CALL ZGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
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
      IF( AINVNM .NE. 0.0D+0 )
     $   ZLA_GBRCOND_X = 1.0D+0 / AINVNM
*
      RETURN
*
      END
