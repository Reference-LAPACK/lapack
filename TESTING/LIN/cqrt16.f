*> \brief \b CQRT16
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CQRT16( TRANS, M, N, NRHS, A, LDA, X, LDX, B, LDB,
*                          RWORK, RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            LDA, LDB, LDX, M, N, NRHS
*       REAL               RESID
*       ..
*       .. Array Arguments ..
*       REAL               RWORK( * )
*       COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CQRT16 computes the residual for a solution of a system of linear
*> equations  A*x = b  or  A'*x = b:
*>    RESID = norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ),
*> where EPS is the machine epsilon.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the form of the system of equations:
*>          = 'N':  A *x = b
*>          = 'T':  A^T*x = b, where A^T is the transpose of A
*>          = 'C':  A^H*x = b, where A^H is the conjugate transpose of A
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of columns of B, the matrix of right hand sides.
*>          NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          The original M x N matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX array, dimension (LDX,NRHS)
*>          The computed solution vectors for the system of linear
*>          equations.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X.  If TRANS = 'N',
*>          LDX >= max(1,N); if TRANS = 'T' or 'C', LDX >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,NRHS)
*>          On entry, the right hand side vectors for the system of
*>          linear equations.
*>          On exit, B is overwritten with the difference B - A*X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  IF TRANS = 'N',
*>          LDB >= max(1,M); if TRANS = 'T' or 'C', LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (M)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is REAL
*>          The maximum over the number of right hand sides of
*>          norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex_lin
*
*  =====================================================================
      SUBROUTINE CQRT16( TRANS, M, N, NRHS, A, LDA, X, LDX, B, LDB,
     $                   RWORK, RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            LDA, LDB, LDX, M, N, NRHS
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            J, N1, N2
      REAL               ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANGE, SCASUM, SLAMCH
      EXTERNAL           LSAME, CLANGE, SCASUM, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if M = 0 or N = 0 or NRHS = 0
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.EQ.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
      IF( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) THEN
         ANORM = CLANGE( 'I', M, N, A, LDA, RWORK )
         N1 = N
         N2 = M
      ELSE
         ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
         N1 = M
         N2 = N
      END IF
*
      EPS = SLAMCH( 'Epsilon' )
*
*     Compute  B - A*X  (or  B - A'*X ) and store in B.
*
      CALL CGEMM( TRANS, 'No transpose', N1, NRHS, N2, -CONE, A, LDA, X,
     $            LDX, CONE, B, LDB )
*
*     Compute the maximum over the number of right hand sides of
*        norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = SCASUM( N1, B( 1, J ), 1 )
         XNORM = SCASUM( N2, X( 1, J ), 1 )
         IF( ANORM.EQ.ZERO .AND. BNORM.EQ.ZERO ) THEN
            RESID = ZERO
         ELSE IF( ANORM.LE.ZERO .OR. XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) /
     $              ( MAX( M, N )*EPS ) )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of CQRT16
*
      END
