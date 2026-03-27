*> \brief \b DORGR2 generates all or part of the orthogonal matrix Q from an RQ factorization determined by sgerqf (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DORGR2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgr2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgr2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgr2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORGR2( M, N, K, A, LDA, TAU, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, M, N
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
*> DORGR2 generates an m by n real matrix Q with orthonormal rows,
*> which is defined as the last m rows of a product of k elementary
*> reflectors of order n
*>
*>       Q  =  H(1) H(2) . . . H(k)
*>
*> as returned by DGERQF.
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the (m-k+i)-th row must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by DGERQF in the last k rows of its array argument
*>          A.
*>          On exit, the m by n matrix Q.
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
*>          reflector H(i), as returned by DGERQF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array. No longer referenced
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument has an illegal value
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
*> \ingroup ungr2
*
*  =====================================================================
      SUBROUTINE DORGR2( M, N, K, A, LDA, TAU, WORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF0C2, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGR2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
*     Note that if N=0, then M must also be 0, so it's sufficient to only test
*     M=0. If we have 0 reflectors, then we define the matrix Q to be the
*     m\times n `identity'
*
      IF( M.LE.0 ) THEN
         RETURN
      ELSE IF( K.LE.0 ) THEN
         CALL DLASET('All', M, N, ZERO, ZERO, A, LDA)
         DO I = N-M+1, N-K
            A(M-N+I,I) = ONE
         END DO
         RETURN
      END IF
*
*     Apply H(1) to the assumed identity matrix from the right
*
      CALL DLARF0C2('Identity', 'Right', 'Backward', 'Rowwise',
     $      M-K, N-K+1, TAU(1), A(M-K+1, 1), LDA, A, LDA)
*
*     Apply H(1) to v_1
*
      CALL DSCAL(N-K, -TAU(1), A(M-K+1, 1), LDA)
      A( M-K+1, N-K+1) = ONE - TAU(1)
      IF( K.GT.1 ) THEN
         DO I = 2, K
*
*           Apply H(i) to A from the right
*
            CALL DLARF0C2('General', 'Right', 'Backward', 'Rowwise',
     $            M-K+I-1, N-K+I, TAU(I), A(M-K+I, 1), LDA, A, LDA)
*
*           Apply H(i) to v_i
*
            CALL DSCAL(N-K+I-1, -TAU(I), A(M-K+I, 1), LDA)
            A( M-K+I, N-K+I) = ONE - TAU(I)
         END DO
      END IF
      RETURN
*
*     End of DORGR2
*
      END
