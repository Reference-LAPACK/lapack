      SUBROUTINE CGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK routine (version 3.?) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- July 2011 --
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDT, M, N, NB
*     ..
*     .. Array Arguments ..
      COMPLEX A( LDA, * ), T( LDT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CGEQRT computes a blocked QR factorization of a complex M-by-N matrix A
*  using the compact WY representation of Q.  
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NB      (input) INTEGER
*          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*          upper triangular if M >= N); the elements below the diagonal
*          are the columns of V.  See below for further details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  T       (output) COMPLEX array, dimension (LDT,MIN(M,N))
*          The upper triangular block reflectors stored in compact form
*          as a sequence of upper triangular blocks.  See below
*          for further details.
*          
*  LDT     (input) INTEGER
*          The leading dimension of the array T.  LDT >= NB.
*
*  WORK    (workspace) COMPLEX array, dimension (NB*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*  
*  The matrix V stores the elementary reflectors H(i) in the i-th column
*  below the diagonal. For example, if M=5 and N=3, the matrix V is
*
*               V = (  1       )
*                   ( v1  1    )
*                   ( v1 v2  1 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  where the vi's represent the vectors which define H(i), which are returned
*  in the matrix A.  The 1's along the diagonal of V are not stored in A.
*
*  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each
*  block is of order NB except for the last block, which is of order 
*  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block
*  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB 
*  for the last block) T's are stored in the NB-by-N matrix T as
*
*               T = (T1 T2 ... TB).
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, IINFO, K
      LOGICAL    USE_RECURSIVE_QR
      PARAMETER( USE_RECURSIVE_QR=.TRUE. )
*     ..
*     .. External Subroutines ..
      EXTERNAL   CGEQRT2, CGEQRT3, CLARFB, XERBLA
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
      ELSE IF( NB.LT.1 .OR. NB.GT.MIN(M,N) ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEQRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) RETURN
*
*     Blocked loop of length K
*
      DO I = 1, K,  NB
         IB = MIN( K-I+1, NB )
*     
*     Compute the QR factorization of the current block A(I:M,I:I+IB-1)
*
         IF( USE_RECURSIVE_QR ) THEN
            CALL CGEQRT3( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         ELSE
            CALL CGEQRT2( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         END IF
         IF( I+IB.LE.N ) THEN
*
*     Update by applying H**H to A(I:M,I+IB:N) from the left
*
            CALL CLARFB( 'L', 'C', 'F', 'C', M-I+1, N-I-IB+1, IB,
     $                   A( I, I ), LDA, T( 1, I ), LDT, 
     $                   A( I, I+IB ), LDA, WORK , N-I-IB+1 )
         END IF
      END DO
      RETURN
*     
*     End of CGEQRT
*
      END
