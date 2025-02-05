c
c
*> \brief \b DORGQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DORGQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORGQR generates an M-by-N real matrix Q with orthonormal columns,
*> which is defined as the first N columns of a product of K elementary
*> reflectors of order M
*>
*>       Q  =  H(1) H(2) . . . H(k)
*>
*> as returned by DGEQRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix Q. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix Q. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines the
*>          matrix Q. N >= K >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the i-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by DGEQRF in the first k columns of its array
*>          argument A.
*>          On exit, the M-by-N matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The first dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DGEQRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is the
*>          optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument has an illegal value
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
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I,II, IB, IINFO, IWS, J,JJ, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL             DLARFB0C2, DLARFT, DORG2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      ! Only need a workspace for dorg2r in case of bailout and 
      ! for the panel factorization
      LWKOPT = MAX( 1, N )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      ! Probably not needed anymore
      NBMIN = 2
      ! Parameter that controls when we cross from blocked to
      ! unblocked
      NX = 0
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = K - 2 * NB
         KK = K - NB
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the only block.
*
      IF( KK.EQ.0 ) THEN
            CALL DORG2R( M, N, K, A, LDA, TAU, WORK, IINFO ) 
      END IF
*
      IF( KK.GT.0 ) THEN
         I = KK + 1
         IB = NB
         ! This is a specialized form of our loop below. We could make this its
         ! own function, however this is a specialized step, so currently we
         ! don't do that.
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
         CALL DLARFT('Forward', 'Column', M-I+1, IB, A(I,I), 
     $                     LDA, TAU(I), A(I,I), LDA)
*
*           Apply H to A(i:m,i+ib:n) from the left
*
*
**        W := V2
*        C1 := V2**T
*
*         Since C1 starts as 0, we are using this instead of WORK(IB+1).
*         This helps us reduce the memory footprint by lowering WORK to
*         be of only size IB
*         CALL DLACPY('All', N-K, IB, A(I+IB,I), LDA,WORK(IB+1),LDWORK)
         DO JJ = K - NB + 1, K
           DO II = K + 1, N
              A( JJ, II ) = A( II, JJ)
           END DO
         END DO
*
*              C1 := T * C1
*
         CALL DTRMM( 'Left', 'Upper', 'No transpose', 'Non-unit', IB,
     $               N-K,ONE, A(I,I), LDA, A(I,I+IB),LDA )
*
*                 C2 := C2 - V2 * C1
*
         CALL DGEMM( 'No transpose', 'No transpose', M-IB-KK, N-K, IB,
     $               -ONE, A( I+IB, I ), LDA, A(I,I+IB),LDA, ZERO,
     $               A( I+IB, I+IB ), LDA )
         DO JJ = 1, N-K
            A(I+IB+JJ-1,I+IB+JJ-1) = 1 + A(I+IB+JJ-1,I+IB+JJ-1)
         END DO
*
*              C1 := -V1 * C1 
*
         CALL DTRMM( 'Left', 'Lower', 'No transpose', 'Unit', IB, N-K,
     $               -ONE, A(I,I), LDA, A(I,I+IB),LDA )
*
*        Apply H to rows i:m of current block
*
         CALL DORG2R(M-I+1, IB, IB, A(I,I), LDA, TAU(I), WORK, IINFO)
         DO I = KI + 1, 1, -NB
            IB = NB
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT('Forward', 'Column', M-I+1, IB, A(I,I), 
     $         LDA, TAU(I), A(I,I), LDA)
*
*           Apply H to A(i:m,i+ib:n) from the left
*
            CALL DLARFB0C2('A', 'A', 'Forward', 'Column', M-I+1, 
     $         N-(I+IB)+1, IB, A(I,I), LDA, A(I,I), LDA, 
     $         A(I,I+IB), LDA)

*
*           Apply H to rows i:m of current block
*
            CALL DORG2R(M-I+1, IB, IB, A(I,I), LDA, TAU(I), WORK, 
     $         IINFO)
         END DO
*        This checks for if K was a perfect multiple of NB
*        so that we only have a special case for the last block when
*        necessary
         IF(I.LT.1) THEN
            IB = I + NB - 1
            I = 1
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT('Forward', 'Column', M-I+1, IB, A(I,I), 
     $         LDA, TAU(I), A(I,I), LDA)
*
*           Apply H to A(i:m,i+ib:n) from the left
*
            CALL DLARFB0C2('A', 'A', 'Forward', 'Column', M-I+1, 
     $         N-(I+IB)+1, IB, A(I,I), LDA, A(I,I), LDA, 
     $         A(I,I+IB),LDA)

*
*           Apply H to rows i:m of current block
*
            CALL DORG2R(M-I+1, IB, IB, A(I,I), LDA, TAU(I), WORK, 
     $         IINFO)
         END IF
      END IF
*
*      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQR
*
      END
