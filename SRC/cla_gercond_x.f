      REAL FUNCTION CLA_GERCOND_X( TRANS, N, A, LDA, AF, LDAF, IPIV, X,
     $                             INFO, WORK, RWORK )
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
      INTEGER            N, LDA, LDAF, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )
      REAL               RWORK( * )
*     ..
*
*  Purpose
*  =======
* 
*     CLA_GERCOND_X computes the infinity norm condition number of
*     op(A) * diag(X) where X is a COMPLEX vector.
*
*  Arguments
*  =========
*
*  X       COMPLEX vector.
*
*  WORK    COMPLEX workspace of size 2*N.
*
*  RWORK   REAL workspace of size 3*N.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            NOTRANS
      INTEGER            KASE
      REAL               AINVNM, ANORM, TMP
      INTEGER            I, J
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
      EXTERNAL           CLACN2, CGETRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, REAL, AIMAG
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      CLA_GERCOND_X = 0.0E+0
*
      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      IF ( .NOT. NOTRANS .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT.
     $     LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLA_GERCOND_X', -INFO )
         RETURN
      END IF
*
*     Compute norm of op(A)*op2(C).
*
      ANORM = 0.0
      IF ( NOTRANS ) THEN
         DO I = 1, N
            TMP = 0.0E+0
            DO J = 1, N
               TMP = TMP + CABS1( A( I, J ) * X( J ) )
            END DO
            RWORK( 2*N+I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0E+0
            DO J = 1, N
               TMP = TMP + CABS1( A( J, I ) * X( J ) )
            END DO
            RWORK( 2*N+I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         CLA_GERCOND_X = 1.0E+0
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
*           Multiply by R.
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( 2*N+I )
            END DO
*
            IF ( NOTRANS ) THEN
               CALL CGETRS( 'No transpose', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
            ELSE
               CALL CGETRS( 'Conjugate transpose', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
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
               CALL CGETRS( 'Conjugate transpose', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
            ELSE
               CALL CGETRS( 'No transpose', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
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
     $   CLA_GERCOND_X = 1.0E+0 / AINVNM
*
      RETURN
*
      END
