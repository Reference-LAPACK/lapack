*> \brief \b DORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by dgeqrf (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DORG2R + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2r.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2r.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2r.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
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
*> DORG2R generates an m by n real matrix Q with orthonormal columns,
*> which is defined as the first n columns of a product of k elementary
*> reflectors of order m
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
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DGEQRF.
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
*> \ingroup ung2r
*
*  =====================================================================
      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
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
      INTEGER            J
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
         CALL XERBLA( 'DORG2R', -INFO )
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
         CALL DLASET('All', M, N, ZERO, ONE, A, LDA)
         RETURN
      END IF
*
*     Apply the first (kth) reflector to the assumed identity matrix from
*     the left. Note that if n=k, we do nothing
*
      CALL DLARF0C2('Identity', 'Left', 'Forward', 'Columnwise',
     $   M-K+1, N-K, TAU(K), A(K+1,K), 1, A(K,K+1), LDA)
*
*     Now we compute the 1st non-zero column of H, which is given by
*     A(k:m,k) = e_k - tau*v_k
*     Analagous to orgkr for n=1 (but T is not used as it is a scalar)
*
      A(K,K) = ONE - TAU(K)
      CALL DSCAL(M-K, -TAU(K), A(K+1,K), 1)
*
*     Now we apply columns 1:k-1 of V to A
*
      IF( K.GT.1 ) THEN
         DO J = K-1, 1, -1
           CALL DLARF0C2('General', 'Left', 'Forward',
     $       'Columnwise',  M-J+1, N-J, TAU(J), A(J+1,J), 1,
     $       A(J,J+1), LDA)
*
*          A(i:m,i) = e_i - tau*v_i
*          Analagous to orgkr for n=1 (but T is not used as it is a scalar)
*
           A(J,J) = ONE - TAU(J)
           CALL DSCAL(M-J, -TAU(J), A(J+1,J), 1)
         END DO
      END IF
      RETURN
*
*     End of DORG2R
*
      END
