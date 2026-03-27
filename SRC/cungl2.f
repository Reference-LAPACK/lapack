*> \brief \b CUNGL2 generates all or part of the unitary matrix Q from an LQ factorization determined by cgelqf (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CUNGL2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cungl2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cungl2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cungl2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CUNGL2( M, N, K, A, LDA, TAU, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, M, N
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
*> CUNGL2 generates an m-by-n complex matrix Q with orthonormal rows,
*> which is defined as the first m rows of a product of k elementary
*> reflectors of order n
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
*>          TAU is COMPLEX array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by CGELQF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array. No longer referenced
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
      SUBROUTINE CUNGL2( M, N, K, A, LDA, TAU, WORK, INFO )
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
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = (1.0E+0, 0.0E+0),
     $                   ZERO = (0.0E+0, 0.0E+0) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARF0C2, CSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
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
         CALL XERBLA( 'CUNGL2', -INFO )
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
         CALL CLASET('All', M, N, ZERO, ONE, A, LDA)
         RETURN
      END IF
*
*     Apply the first (kth) reflector to the assumed identity matrix from
*     the right. Note that if m=k, we do nothing
*
      CALL CLARF0C2('Identity', 'Right', 'Forward', 'Rowwise',
     $   M-K, N-K+1, CONJG(TAU(K)), A(K,K+1), LDA, A(K+1,K), LDA)
*
*     Now we compute the 1st non-zero row of H, which is given by
*     A(k,k:n) = (e_k - tau*v_k)'
*     Analagous to orglk for n=1 (but T is not used as it is a scalar)
*
      A(K,K) = ONE - CONJG(TAU(K))
      CALL CSCAL(N-K, -CONJG(TAU(K)), A(K,K+1), LDA)
      IF( K.GT.1 ) THEN
         DO I = K-1, 1, -1
            CALL CLARF0C2('General', 'Right', 'Forward', 'Rowwise',
     $         M-I, N-I+1, CONJG(TAU(I)), A(I,I+1), LDA, A(I+1,I), LDA)
*
*           A(i,i:n) = (e_i - tau*v_i)'
*           Analagous to orglk for n=1 (but T is not used as it is a scalar)
*
            A(I,I) = ONE - CONJG(TAU(I))
            CALL CSCAL(N-I, -CONJG(TAU(I)), A(I,I+1), LDA)
         END DO
      END IF
      RETURN
*
*     End of CUNGL2
*
      END
