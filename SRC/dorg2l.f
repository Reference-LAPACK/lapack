*> \brief \b DORG2L generates all or part of the orthogonal matrix Q from a QL factorization determined by sgeqlf (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DORG2L + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2l.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2l.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2l.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
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
*> DORG2L generates an m by n real matrix Q with orthonormal columns,
*> which is defined as the last n columns of a product of k elementary
*> reflectors of order m
*>
*>       Q  =  H(k) . . . H(2) H(1)
*>
*> as returned by DGEQLF.
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
*>          On entry, the (n-k+i)-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by DGEQLF in the last k columns of its array
*>          argument A.
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
*>          reflector H(i), as returned by DGEQLF.
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
*> \ingroup ung2l
*
*  =====================================================================
      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
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
      INTEGER            J, L
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
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2L', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
*     Note that if M=0, then N must also be 0, so it's sufficient to only test
*     N=0. If we have 0 reflectors, then we define the matrix Q to be the
*     m\times n `identity'
*
      IF( N.LE.0 ) THEN
         RETURN
      ELSE IF( K.LE.0 ) THEN
         CALL DLASET('ALL', M, N, ZERO, ZERO, A, LDA)
         DO J = 1, N
            A(M-N+J,J) = ONE
         END DO
         RETURN
      END IF
*
*     Apply H(1) to the assumed identity matrix from the left
*
      CALL DLARF0C2('Identity', 'Left', 'Backward', 'Columnwise',
     $      M-K+1, N-K, TAU(1), A(1, N-K+1), 1, A, LDA)
*
*     Apply H(1) to v_1
*
      CALL DSCAL(M-K, -TAU(1), A(1, N-K+1), 1)
      A(M-K+1,N-K+1) = ONE - TAU(1)
      DO L = M - K + 1 + 1, M
         A( L, N-K+1 ) = ZERO
      END DO
      IF( K.GT.1 ) THEN
         DO J = 2, K
*
*           Apply H(j) to the leading n-k+j-1 columns of Q from the left
*
            CALL DLARF0C2('General', 'Left', 'Backward',
     $         'Columnwise', M-K+J, N-K+J-1, TAU(J), A(1, N-K+J), 1,
     $         A, LDA)
*
*           Apply H(j) to v_j
*
            CALL DSCAL(M-K+J-1, -TAU(J), A(1,N-K+J), 1)
            A(M-K+J,N-K+J) = ONE - TAU(J)
            DO L = M - K + J + 1, M
               A( L, N-K+J ) = ZERO
            END DO
         END DO
      END IF
      RETURN
*
*     End of DORG2L
*
      END
