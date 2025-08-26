*> \brief \b CUNGLQ
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CUNGLQ + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunglq.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunglq.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunglq.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CUNGLQ generates an M-by-N complex matrix Q with orthonormal rows,
*> which is defined as the first M rows of a product of K elementary
*> reflectors of order N
*>
*>       Q  =  H(k)**H . . . H(2)**H H(1)**H
*>
*> as returned by CGELQF.
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
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the i-th row must contain the vector which defines
*>          the elementary reflector H(i), for i = 1,2,...,k, as returned
*>          by CGELQF in the first k rows of its array argument A.
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
*>          TAU is COMPLEX array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by CGELQF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,M).
*>          For optimum performance LWORK >= M*NB, where NB is
*>          the optimal blocksize.
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
*>          = 0:  successful exit;
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
*> \ingroup unglq
*
*  =====================================================================
      SUBROUTINE CUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, KI, KK, LWKOPT, LDWORK,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARFB0C2, CLARFT, CUNGL2,
     $                   CUNGLK, XERBLA
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
      NB = ILAENV( 1, 'CUNGLQ', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, M )*NB
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNGLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = MAX( 2, ILAENV( 2, 'CUNGLQ', ' ', M, N, K, -1 ) )
      NX = MAX( 0, ILAENV( 3, 'CUNGLQ', ' ', M, N, K, -1 ) )
      IWS = M
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Treat the last NB block starting at KK+1 specially then use our blocking
*        method from the block starting at KI+1 to 1 
*
         KI = K - 2 * NB
         KK = K - NB
      ELSE
         KK = 0
      END IF
*
*     Potentially bail to the unblocked version
*
      IF( KK.EQ.0 ) THEN
         CALL CUNGL2( M, N, K, A, LDA, TAU, WORK, IINFO )
      END IF
*
      IF( KK.GT.0 ) THEN
*
*        Factor the last block assuming that our first application
*        will be on the Identity matrix
*
         I = KK + 1
         IB = NB
*
*        Form the triangular factor of the block reflector
*        H = H(i) H(i+1) . . . H(i+ib-1)
*
         CALL CLARFT( 'Forward', 'Transpose', N-I+1, IB, A( I, I ),
     $               LDA, TAU( I ), A( I, I ), LDA )
*
*        Apply H to A(i+ib:m,i:n) from the right
*
         CALL CLARFB0C2(.TRUE., 'Right', 'No Transpose', 'Forward',
     $         'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I), LDA, A(I,I), 
     $         LDA, A(I+IB,I), LDA)
*
*        Apply H to columns i:n of current block

         CALL CUNGLK( IB, N-I+1, A( I, I ), LDA)
*
*        Use blocked code
*
         DO I = KI + 1, 1, -NB
            IB = NB
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL CLARFT( 'Forward', 'Transpose', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL CLARFB0C2(.FALSE., 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL CUNGLK( IB, N-I+1, A( I, I ), LDA)
         END DO
*
*        This checks for if K was a perfect multiple of NB
*        so that we only have a special case for the last block when
*        necessary
*
         IF(I.LT.1) THEN
            IB = I + NB - 1
            I = 1
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL CLARFT( 'Forward', 'Transpose', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL CLARFB0C2(.FALSE., 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL CUNGLK( IB, N-I+1, A( I, I ), LDA)
         END IF
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK(IWS)
      RETURN
*
*     End of CUNGLQ
*
      END
