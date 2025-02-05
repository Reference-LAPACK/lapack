*> \brief \b CPOEQUB
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CPOEQUB + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpoequb.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpoequb.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpoequb.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, N
*       REAL               AMAX, SCOND
*       ..
*       .. Array Arguments ..
*       COMPLEX            A( LDA, * )
*       REAL               S( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CPOEQUB computes row and column scalings intended to equilibrate a
*> Hermitian positive definite matrix A and reduce its condition number
*> (with respect to the two-norm).  S contains the scale factors,
*> S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with
*> elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This
*> choice of S puts the condition number of B within a factor N of the
*> smallest possible condition number over all possible diagonal
*> scalings.
*>
*> This routine differs from CPOEQU by restricting the scaling factors
*> to a power of the radix.  Barring over- and underflow, scaling by
*> these factors introduces no additional rounding errors.  However, the
*> scaled diagonal entries are no longer approximately 1 but lie
*> between sqrt(radix) and 1/sqrt(radix).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          The N-by-N Hermitian positive definite matrix whose scaling
*>          factors are to be computed.  Only the diagonal elements of A
*>          are referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is REAL array, dimension (N)
*>          If INFO = 0, S contains the scale factors for A.
*> \endverbatim
*>
*> \param[out] SCOND
*> \verbatim
*>          SCOND is REAL
*>          If INFO = 0, S contains the ratio of the smallest S(i) to
*>          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
*>          large nor too small, it is not worth scaling by S.
*> \endverbatim
*>
*> \param[out] AMAX
*> \verbatim
*>          AMAX is REAL
*>          Absolute value of largest matrix element.  If AMAX is very
*>          close to overflow or very close to underflow, the matrix
*>          should be scaled.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the i-th diagonal element is nonpositive.
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
*> \ingroup poequb
*
*  =====================================================================
      SUBROUTINE CPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, N
      REAL               AMAX, SCOND
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * )
      REAL               S( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      REAL               SMIN, BASE, TMP
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT, LOG, INT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
*     Positive definite only performs 1 pass of equilibration.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPOEQUB', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         SCOND = ONE
         AMAX = ZERO
         RETURN
      END IF

      BASE = SLAMCH( 'B' )
      TMP = -0.5 / LOG ( BASE )
*
*     Find the minimum and maximum diagonal elements.
*
      S( 1 ) = REAL( A( 1, 1 ) )
      SMIN = S( 1 )
      AMAX = S( 1 )
      DO 10 I = 2, N
         S( I ) = REAL( A( I, I ) )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
   10 CONTINUE
*
      IF( SMIN.LE.ZERO ) THEN
*
*        Find the first non-positive diagonal element and return.
*
         DO 20 I = 1, N
            IF( S( I ).LE.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   20    CONTINUE
      ELSE
*
*        Set the scale factors to the reciprocals
*        of the diagonal elements.
*
         DO 30 I = 1, N
            S( I ) = BASE ** INT( TMP * LOG( S( I ) ) )
   30    CONTINUE
*
*        Compute SCOND = min(S(I)) / max(S(I)).
*
         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      END IF
*
      RETURN
*
*     End of CPOEQUB
*
      END
