*> \brief \b CUNGQL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CUNGQL + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cungql.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cungql.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cungql.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
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
*> CUNGQL generates an M-by-N complex matrix Q with orthonormal columns,
*> which is defined as the last N columns of a product of K elementary
*> reflectors of order M
*>
*>       Q  =  H(k) . . . H(2) H(1)
*>
*> as returned by CGEQLF.
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
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the (n-k+i)-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by CGEQLF in the last k columns of its array
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
*>          TAU is COMPLEX array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by CGEQLF.
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
*>          The dimension of the array WORK. LWORK >= max(1,N).
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
*> \ingroup ungql
*
*  =====================================================================
      SUBROUTINE CUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
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
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, KK, LWKOPT, NB, NBMIN,
     $                   NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARFB0C2, CLARFT, CUNG2L,
     $                   CUNGKL, XERBLA
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
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'CUNGQL', ' ', M, N, K, -1 )
*
*           Only need a workspace for cung2l in case of bailout
*           and for the panel factorization
*
            LWKOPT = N
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
         IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      NX = MAX( 0, ILAENV( 3, 'CUNGQL', ' ', M, N, K, -1 ))
      IWS = N
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        We use blocked code for the entire construction
*
         KK = K
      ELSE
         KK = 0
      END IF
*
*     Possibly bail to the unblocked code.
*
      IF( KK.EQ.0 ) THEN
         CALL CUNG2L( M, N, K, A, LDA, TAU, WORK, IINFO )
      END IF
*
      IF( KK.GT.0 ) THEN
*
*        Factor the first block assuming that our first application
*        will be on the Identity matrix
*
         I = 1
         IB = NB
*
*        Form the triangular factor of the block reflector
*        H = H(i+ib-1) . . . H(i+1) H(i)
*
         CALL CLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                  A( 1, N-K+I ), LDA, TAU( I ),
     $                  A( M-K+I, N-K+I ), LDA)
*
*        Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
*        Exploit the fact that we are applying to an identity 
*
         CALL CLARFB0C2(.TRUE., 'Left', 'No Transpose', 'Backward', 
     $         'Columnwise', M-K+I+IB-1, N-K+I-1, IB, A(1, N-K+I), 
     $         LDA, A( M-K+I, N-K+I ), LDA, A, LDA)
*
*        Apply H to rows 1:m-k+i+ib-1 of current block
*
         CALL CUNGKL( M-K+I+IB-1, IB, A( 1, N-K+I ), LDA)

*        Use blocked code on the remaining blocks if there are any.
*
         DO I = NB+1, K, NB
*
*           The last block may be less than size NB
*
            IB = MIN(NB, K-I+1)
*
*           Form the triangular factor of the block reflector
*           H = H(i+ib-1) . . . H(i+1) H(i)
*
            CALL CLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                  A( 1, N-K+I ), LDA, TAU( I ), 
     $                  A( M-K+I, N-K+I ), LDA )
*
*           Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
*
            CALL CLARFB0C2(.FALSE., 'Left', 'No Transpose',
     $            'Backward', 'Columnwise', M-K+I+IB-1, N-K+I-1, IB, 
     $            A(1, N-K+I), LDA, A( M-K+I, N-K+I ), LDA, A, LDA)
*
*           Apply H to rows 1:m-k+i+ib-1 of current block
*
            CALL CUNGKL( M-K+I+IB-1, IB, A( 1, N-K+I ), LDA)
         END DO
      END IF
*
      WORK( 1 ) = CMPLX( IWS )
      RETURN
*
*     End of CUNGQL
*
      END
