*> \brief \b CRSCL multiplies a vector by the reciprocal of a real scalar.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CRSCL + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/crscl.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/crscl.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/crscl.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CRSCL( N, A, X, INCX )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       COMPLEX            A
*       ..
*       .. Array Arguments ..
*       COMPLEX            X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CRSCL multiplies an n-element complex vector x by the complex scalar
*> 1/a.  This is done without overflow or underflow as long as
*> the final result x/a does not overflow or underflow.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of components of the vector x.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX
*>          The scalar a which is used to divide each component of x.
*>          A must not be 0, or the subroutine will divide by zero.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX array, dimension
*>                         (1+(N-1)*abs(INCX))
*>          The n-element vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive values of the vector X.
*>          > 0:  X(1) = X(1) and X(1+(i-1)*INCX) = x(i),     1< i<= n
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
*> \ingroup complexOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE CRSCL( N, A, X, INCX )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX            A
*     ..
*     .. Array Arguments ..
      COMPLEX            X( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      REAL               BIGNUM, SMLNUM, HUGE, AR, AI, ABSR, ABSI, UR
     %                   , UI
      COMPLEX            INVA
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      COMPLEX            CLADIV
      EXTERNAL           SLAMCH, CLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           CSCAL, CSSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      HUGE   = SLAMCH( 'O' )
*
*     Initialize constants related to A.
*
      AR = REAL( A )
      AI = AIMAG( A )
      ABSR = ABS( AR )
      ABSI = ABS( AI )
*
      IF( ABSI.EQ.ZERO ) THEN
*        If alpha is real, then we can use csrscl
         CALL CSRSCL( N, AR, X, INCX )
*
      ELSE IF( ABSR.EQ.ZERO ) THEN
*        If alpha has a zero real part, then we follow the same rules as if
*        alpha were real.
         IF( ABSI.GT.BIGNUM ) THEN
            INVA = CMPLX( ZERO, -BIGNUM / AI )
            CALL CSSCAL( N, SMLNUM, X, INCX )
            CALL CSCAL( N, INVA, X, INCX )
         ELSE IF( ABSI.LT.SMLNUM ) THEN
            INVA = CMPLX( ZERO, -SMLNUM / AI )
            CALL CSCAL( N, INVA, X, INCX )
            CALL CSSCAL( N, BIGNUM, X, INCX )
         ELSE
            INVA = CMPLX( ZERO, -ONE / AI )
            CALL CSCAL( N, INVA, X, INCX )
         END IF
*
      ELSE IF( (ABSR.GE.BIGNUM).OR.(ABSI.GE.BIGNUM) ) THEN
*        Either real or imaginary part is too large.
         INVA = CLADIV( CMPLX( BIGNUM, ZERO ), A )
         CALL CSSCAL( N, SMLNUM, X, INCX )
         CALL CSCAL( N, INVA, X, INCX )
*
      ELSE
*        The following numbers can be computed without NaNs and zeros.
*        They do not overflow simultaneously.
*        They are the inverse of the real and imaginary parts of 1/alpha.
         UR = AR + AI * ( AI / AR )
         UI = AI + AR * ( AR / AI )
*
         IF( (ABS( UR ).LT.SMLNUM).OR.(ABS( UI ).LT.SMLNUM) ) THEN
            INVA = CMPLX( SMLNUM / UR, -SMLNUM / UI )
            CALL CSCAL( N, INVA, X, INCX )
            CALL CSSCAL( N, BIGNUM, X, INCX )
         ELSE IF( ABS( UR ).GT.HUGE ) THEN
            IF( ABSR.GE.ABSI ) THEN
               UR = (SMLNUM * AR) + AI * (SMLNUM * (AI / AR))
            ELSE
               UR = (SMLNUM * AR) + AI * ((SMLNUM * AI) / AR)
            END IF
            INVA = CMPLX( ONE / UR, -BIGNUM / UI )
            CALL CSSCAL( N, SMLNUM, X, INCX )
            CALL CSCAL( N, INVA, X, INCX )
         ELSE IF( ABS( UI ).GT.HUGE ) THEN
            IF( ABSI.GE.ABSR ) THEN
               UI = (SMLNUM * AI) + AR * (SMLNUM * (AR / AI))
            ELSE
               UI = (SMLNUM * AI) + AR * ((SMLNUM * AR) / AI)
            END IF
            INVA = CMPLX( BIGNUM / UR, -ONE / UI )
            CALL CSSCAL( N, SMLNUM, X, INCX )
            CALL CSCAL( N, INVA, X, INCX )
         ELSE IF( (ABS( UR ).GT.BIGNUM).OR.(ABS( UI ).GT.BIGNUM) ) THEN
            INVA = CMPLX( BIGNUM / UR, -BIGNUM / UI )
            CALL CSSCAL( N, SMLNUM, X, INCX )
            CALL CSCAL( N, INVA, X, INCX )
         ELSE
            INVA = CMPLX( ONE / UR, -ONE / UI )
            CALL CSCAL( N, INVA, X, INCX )
         END IF
      END IF
*
      RETURN
*
*     End of CRSCL
*
      END
