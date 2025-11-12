*> \brief \b SORGRQ
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SORGRQ + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgrq.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgrq.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgrq.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SORGRQ generates an M-by-N real matrix Q with orthonormal rows,
*> which is defined as the last M rows of a product of K elementary
*> reflectors of order N
*>
*>       Q  =  H(1) H(2) . . . H(k)
*>
*> as returned by SGERQF.
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
*>          The number of columns of the matrix Q. N >= M.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines the
*>          matrix Q. M >= K >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the (m-k+i)-th row must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by SGERQF in the last k rows of its array argument
*>          A.
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
*>          TAU is REAL array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by SGERQF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,M).
*>          For optimum performance LWORK >= M*NB, where NB is the
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
*> \ingroup ungrq
*
*  =====================================================================
      SUBROUTINE SORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, II, IINFO, IWS, KK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARFB0C2, SLARFT, SORGR2,
     $                   SORGRK, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      REAL               SROUNDUP_LWORK
      EXTERNAL           ILAENV, SROUNDUP_LWORK
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'SORGRQ', ' ', M, N, K, -1 )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
*
      IF( INFO.EQ.0 ) THEN
         LWKOPT = MAX(1,M)
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
         IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORGRQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      NX = MAX( 0, ILAENV( 3, 'SORGRQ', ' ', M, N, K, -1 ) )
      IWS = M
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        We want to use the blocking method as long as our matrix is big enough
*        and it's deemed worthwhile
*
         KK = K
      ELSE
         KK = 0
      END IF
*
*     Potentially bail to the unblocked code
*
      IF( KK.EQ.0 ) THEN
         CALL SORGR2( M, N, K, A, LDA, TAU, WORK, IINFO )
      END IF
*
      IF( KK.GT.0 ) THEN
*
*        Factor the first block assuming that our first application
*        will be on the Identity matrix
*
         I = 1
         IB = NB
         II = M - K + I
*
*        Form the triangular factor of the block reflector
*        H = H(i+ib-1) . . . H(i+1) H(i)
*
         CALL SLARFT( 'Transpose', 'Rowwise', N-K+I+IB-1, IB,
     $                A( II, 1 ), LDA, TAU( I ), A( II, N-K+I ), LDA )
*
*        Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
*
         CALL SLARFB0C2(.TRUE., 'Right', 'No Transpose', 'Backward', 
     $         'Rowwise', II-1, N-K+I+IB-1, IB, A(II,1), LDA,
     $          A( II, N-K+I ), LDA, A, LDA)
*
*           Apply H to columns 1:n-k+i+ib-1 of current block
*
         CALL SORGRK( IB, N-K+I+IB-1, A( II, 1 ), LDA )

         DO I = NB + 1, K, NB
*
*           The last block may be less than size NB
*
            IB = MIN(NB, K-I+1)
            II = M - K + I
*
*           Form the triangular factor of the block reflector
*           H = H(i+ib-1) . . . H(i+1) H(i)
*
            CALL SLARFT( 'Transpose', 'Rowwise', N-K+I+IB-1, IB,
     $                A( II, 1 ), LDA, TAU( I ), A( II, N-K+I ), LDA )
*
*           Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
*
            CALL SLARFB0C2(.FALSE., 'Right', 'No Transpose',
     $            'Backward', 'Rowwise', II-1, N-K+I+IB-1, IB, A(II,1),
     $             LDA, A( II, N-K+I ), LDA, A, LDA)
*
*           Apply H to columns 1:n-k+i+ib-1 of current block
*
            CALL SORGRK( IB, N-K+I+IB-1, A( II, 1 ), LDA )
         END DO
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK(IWS)
      RETURN
*
*     End of SORGRQ
*
      END
