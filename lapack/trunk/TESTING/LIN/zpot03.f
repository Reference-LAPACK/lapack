      SUBROUTINE ZPOT03( UPLO, N, A, LDA, AINV, LDAINV, WORK, LDWORK,
     $                   RWORK, RCOND, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDAINV, LDWORK, N
      DOUBLE PRECISION   RCOND, RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), AINV( LDAINV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  ZPOT03 computes the residual for a Hermitian matrix times its
*  inverse:
*     norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS ),
*  where EPS is the machine epsilon.
*
*  Arguments
*  ==========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The original Hermitian matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N)
*
*  AINV    (input/output) COMPLEX*16 array, dimension (LDAINV,N)
*          On entry, the inverse of the matrix A, stored as a Hermitian
*          matrix in the same format as A.
*          In this version, AINV is expanded into a full matrix and
*          multiplied by A, so the opposing triangle of AINV will be
*          changed; i.e., if the upper triangular part of AINV is
*          stored, the lower triangular part will be used as work space.
*
*  LDAINV  (input) INTEGER
*          The leading dimension of the array AINV.  LDAINV >= max(1,N).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,N)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.  LDWORK >= max(1,N).
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RCOND   (output) DOUBLE PRECISION
*          The reciprocal of the condition number of A, computed as
*          ( 1/norm(A) ) / norm(AINV).
*
*  RESID   (output) DOUBLE PRECISION
*          norm(I - A*AINV) / ( N * norm(A) * norm(AINV) * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   AINVNM, ANORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANHE
      EXTERNAL           LSAME, DLAMCH, ZLANGE, ZLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZHEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCONJG
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RCOND = ONE
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
      AINVNM = ZLANHE( '1', UPLO, N, AINV, LDAINV, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
*     Expand AINV into a full matrix and call ZHEMM to multiply
*     AINV on the left by A.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, J - 1
               AINV( J, I ) = DCONJG( AINV( I, J ) )
   10       CONTINUE
   20    CONTINUE
      ELSE
         DO 40 J = 1, N
            DO 30 I = J + 1, N
               AINV( J, I ) = DCONJG( AINV( I, J ) )
   30       CONTINUE
   40    CONTINUE
      END IF
      CALL ZHEMM( 'Left', UPLO, N, N, -CONE, A, LDA, AINV, LDAINV,
     $            CZERO, WORK, LDWORK )
*
*     Add the identity matrix to WORK .
*
      DO 50 I = 1, N
         WORK( I, I ) = WORK( I, I ) + CONE
   50 CONTINUE
*
*     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = ZLANGE( '1', N, N, WORK, LDWORK, RWORK )
*
      RESID = ( ( RESID*RCOND ) / EPS ) / DBLE( N )
*
      RETURN
*
*     End of ZPOT03
*
      END
