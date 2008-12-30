      REAL FUNCTION CLA_HERCOND_C( UPLO, N, A, LDA, AF, LDAF, IPIV, C,
     $                             CAPPLY, INFO, WORK, RWORK )
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
      LOGICAL            CAPPLY
      INTEGER            N, LDA, LDAF, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * )
      REAL               C ( * ), RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*     CLA_HERCOND_C computes the infinity norm condition number of
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
      INTEGER            KASE, I, J
      REAL               AINVNM, ANORM, TMP
      LOGICAL            UP
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
      EXTERNAL           CLACN2, CHETRS, XERBLA
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
*
      CLA_HERCOND_C = 0.0E+0
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLA_HERCOND_C', -INFO )
         RETURN
      END IF
      UP = .FALSE.
      IF ( LSAME( UPLO, 'U' ) ) UP = .TRUE.
*
*     Compute norm of op(A)*op2(C).
*
      ANORM = 0.0E+0
      IF ( UP ) THEN
         DO I = 1, N
            TMP = 0.0E+0
            IF ( CAPPLY ) THEN
               DO J = 1, N
                  IF ( I.GT.J ) THEN
                     TMP = TMP + CABS1( A( J, I ) ) / C( J )
                  ELSE
                     TMP = TMP + CABS1( A( I, J ) ) / C( J )
                  END IF
               END DO
            ELSE
               DO J = 1, N
                  IF ( I.GT.J ) THEN
                     TMP = TMP + CABS1( A( J, I ) )
                  ELSE
                     TMP = TMP + CABS1( A( I, J ) )
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
                  IF ( I.LT.J ) THEN
                     TMP = TMP + CABS1( A( J, I ) ) / C( J )
                  ELSE
                     TMP = TMP + CABS1( A( I, J ) ) / C( J )
                  END IF
               END DO
            ELSE
               DO J = 1, N
                  IF ( I.LT.J ) THEN
                     TMP = TMP + CABS1( A( J, I ) )
                  ELSE
                     TMP = TMP + CABS1( A( I, J ) )
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
         CLA_HERCOND_C = 1.0E+0
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
            IF ( UP ) THEN
               CALL CHETRS( 'U', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
            ELSE
               CALL CHETRS( 'L', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
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
            IF ( UP ) THEN
               CALL CHETRS( 'U', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
            ELSE
               CALL CHETRS( 'L', N, 1, AF, LDAF, IPIV,
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
     $   CLA_HERCOND_C = 1.0E+0 / AINVNM
*
      RETURN
*
      END
