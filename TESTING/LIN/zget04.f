*> \brief \b ZGET04
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGET04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
*
*       .. Scalar Arguments ..
*       INTEGER            LDX, LDXACT, N, NRHS
*       DOUBLE PRECISION   RCOND, RESID
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         X( LDX, * ), XACT( LDXACT, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGET04 computes the difference between a computed solution and the
*> true solution to a system of linear equations.
*>
*> RESID =  ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS ),
*> where RCOND is the reciprocal of the condition number and EPS is the
*> machine epsilon.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows of the matrices X and XACT.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of columns of the matrices X and XACT.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)
*>          The computed solution vectors.  Each vector is stored as a
*>          column of the matrix X.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X.  LDX >= max(1,N).
*> \endverbatim
*>
*> \param[in] XACT
*> \verbatim
*>          XACT is COMPLEX*16 array, dimension (LDX,NRHS)
*>          The exact solution vectors.  Each vector is stored as a
*>          column of the matrix XACT.
*> \endverbatim
*>
*> \param[in] LDXACT
*> \verbatim
*>          LDXACT is INTEGER
*>          The leading dimension of the array XACT.  LDXACT >= max(1,N).
*> \endverbatim
*>
*> \param[in] RCOND
*> \verbatim
*>          RCOND is DOUBLE PRECISION
*>          The reciprocal of the condition number of the coefficient
*>          matrix in the system of equations.
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is DOUBLE PRECISION
*>          The maximum over the NRHS solution vectors of
*>          ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS )
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
*> \ingroup complex16_lin
*
*  =====================================================================
      SUBROUTINE ZGET04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDX, LDXACT, N, NRHS
      DOUBLE PRECISION   RCOND, RESID
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, J
      DOUBLE PRECISION   DIFFNM, EPS, XNORM
      COMPLEX*16         ZDUM
*     ..
*     .. External Functions ..
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           IZAMAX, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if RCOND is invalid.
*
      EPS = DLAMCH( 'Epsilon' )
      IF( RCOND.LT.ZERO ) THEN
         RESID = 1.0D0 / EPS
         RETURN
      END IF
*
*     Compute the maximum of
*        norm(X - XACT) / ( norm(XACT) * EPS )
*     over all the vectors X and XACT .
*
      RESID = ZERO
      DO 20 J = 1, NRHS
         IX = IZAMAX( N, XACT( 1, J ), 1 )
         XNORM = CABS1( XACT( IX, J ) )
         DIFFNM = ZERO
         DO 10 I = 1, N
            DIFFNM = MAX( DIFFNM, CABS1( X( I, J )-XACT( I, J ) ) )
   10    CONTINUE
         IF( XNORM.LE.ZERO ) THEN
            IF( DIFFNM.GT.ZERO )
     $         RESID = 1.0D0 / EPS
         ELSE
            RESID = MAX( RESID, ( DIFFNM / XNORM )*RCOND )
         END IF
   20 CONTINUE
      IF( RESID*EPS.LT.1.0D0 )
     $   RESID = RESID / EPS
*
      RETURN
*
*     End of ZGET04
*
      END
