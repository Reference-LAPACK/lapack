*> \brief \b CGEQL2 computes the QL factorization of a general rectangular matrix using an unblocked algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CGEQL2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeql2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeql2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeql2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CGEQL2( M, N, A, LDA, TAU, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
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
*> CGEQL2 computes a QL factorization of a complex m by n matrix A:
*> A = Q * L.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the m by n matrix A.
*>          On exit, if m >= n, the lower triangle of the subarray
*>          A(m-n+1:m,1:n) contains the n by n lower triangular matrix L;
*>          if m <= n, the elements on and below the (n-m)-th
*>          superdiagonal contain the m by n lower trapezoidal matrix L;
*>          the remaining elements, with the array TAU, represent the
*>          unitary matrix Q as a product of elementary reflectors
*>          (see Further Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
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
*> \ingroup geql2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of elementary reflectors
*>
*>     Q = H(k) . . . H(2) H(1), where k = min(m,n).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**H
*>
*>  where tau is a complex scalar, and v is a complex vector with
*>  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in
*>  A(1:m-k+i-1,n-k+i), and tau in TAU(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CGEQL2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, K
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARF1L, CLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEQL2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = K, 1, -1
*
*        Generate elementary reflector H(i) to annihilate
*        A(1:m-k+i-1,n-k+i)
*
         CALL CLARFG( M-K+I, A( M-K+I, N-K+I ), A( 1, N-K+I ), 1,
     $                TAU( I ) )
*
*        Apply H(i)**H to A(1:m-k+i,1:n-k+i-1) from the left
*
         CALL CLARF1L( 'Left', M-K+I, N-K+I-1, A( 1, N-K+I ), 1,
     $                 CONJG( TAU( I ) ), A, LDA, WORK )
   10 CONTINUE
      RETURN
*
*     End of CGEQL2
*
      END
