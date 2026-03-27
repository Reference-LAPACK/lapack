*> \brief \b SORGL2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, M, N
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
*> SORGL2 generates an m by n real matrix Q with orthonormal rows,
*> which is defined as the first m rows of a product of k elementary
*> reflectors of order n
*>
*>       Q  =  H(k) . . . H(2) H(1)
*>
*> as returned by SGELQF.
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
*>          On entry, the i-th row must contain the vector which defines
*>          the elementary reflector H(i), for i = 1,2,...,k, as returned
*>          by SGELQF in the first k rows of its array argument A.
*>          On exit, the m-by-n matrix Q.
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
*>          reflector H(i), as returned by SGELQF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array. No longer referenced
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
*> \ingroup ungl2
*
*  =====================================================================
      SUBROUTINE SORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
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
      REAL               A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARF0C2, SSCAL, XERBLA
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
         CALL XERBLA( 'SORGL2', -INFO )
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
         CALL SLASET('All', M, N, ZERO, ONE, A, LDA)
         RETURN
      END IF
*
*     Apply the first (kth) reflector to the assumed identity matrix from
*     the right. Note that if m=k, we do nothing
*
      CALL SLARF0C2('Identity', 'Right', 'Forward', 'Rowwise',
     $   M-K, N-K+1, TAU(K), A(K,K+1), LDA, A(K+1,K), LDA)
*
*     Now we compute the 1st non-zero row of H, which is given by
*     A(k,k:n) = (e_k - tau*v_k)'
*     Analagous to orglk for n=1 (but T is not used as it is a scalar)
*
      A(K,K) = ONE - TAU(K)
      CALL SSCAL(N-K, -TAU(K), A(K,K+1), LDA)
      IF( K.GT.1 ) THEN
         DO I = K-1, 1, -1
            CALL SLARF0C2('General', 'Right', 'Forward', 'Rowwise',
     $         M-I, N-I+1, TAU(I), A(I,I+1), LDA, A(I+1,I), LDA)
*
*           A(i,i:n) = (e_i - tau*v_i)'
*           Analagous to orglk for n=1 (but T is not used as it is a scalar)
*
            A(I,I) = ONE - TAU(I)
            CALL SSCAL(N-I, -TAU(I), A(I,I+1), LDA)
         END DO
      END IF
      RETURN
*
*     End of SORGL2
*
      END
