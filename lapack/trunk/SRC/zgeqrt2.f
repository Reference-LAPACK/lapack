      SUBROUTINE ZGEQRT2( M, N, A, LDA, T, LDT, INFO )
*
*  -- LAPACK routine (version 3.?) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- July 2011                                                      --
*
*     .. Scalar Arguments ..
      INTEGER   INFO, LDA, LDT, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16   A( LDA, * ), T( LDT, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEQRT2 computes a QR factorization of a complex M-by-N matrix A, 
*  using the compact WY representation of Q. 
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= N.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the complex M-by-N matrix A.  On exit, the elements on and
*          above the diagonal contain the N-by-N upper triangular matrix R; the
*          elements below the diagonal are the columns of V.  See below for
*          further details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  T       (output) COMPLEX*16 array, dimension (LDT,N)
*          The N-by-N upper triangular factor of the block reflector.
*          The elements on and above the diagonal contain the block
*          reflector T; the elements below the diagonal are not used.
*          See below for further details.
*
*  LDT     (intput) INTEGER
*          The leading dimension of the array T.  LDT >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
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
*  in the matrix A.  The 1's along the diagonal of V are not stored in A.  The
*  block reflector H is then given by
*
*               H = I - V * T * V**H
*
*  where V**H is the conjugate transpose of V.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16  ONE, ZERO
      PARAMETER( ONE = (1.0D+00,0.0D+00), ZERO = (0.0D+00,0.0D+00) )
*     ..
*     .. Local Scalars ..
      INTEGER   I, K
      COMPLEX*16   AII, ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL  ZLARFG, ZGEMV, ZGERC, ZTRMV, XERBLA
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
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQRT2', -INFO )
         RETURN
      END IF
*      
      K = MIN( M, N )
*
      DO I = 1, K
*
*        Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1)
*
         CALL ZLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                T( I, 1 ) )
         IF( I.LT.N ) THEN
*
*           Apply H(i) to A(I:M,I+1:N) from the left
*
            AII = A( I, I )
            A( I, I ) = ONE
*
*           W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)]
*
            CALL ZGEMV( 'C',M-I+1, N-I, ONE, A( I, I+1 ), LDA, 
     $                  A( I, I ), 1, ZERO, T( 1, N ), 1 )
*
*           A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H
*
            ALPHA = -CONJG(T( I, 1 ))
            CALL ZGERC( M-I+1, N-I, ALPHA, A( I, I ), 1, 
     $           T( 1, N ), 1, A( I, I+1 ), LDA )
            A( I, I ) = AII
         END IF
      END DO
*
      DO I = 2, N
         AII = A( I, I )
         A( I, I ) = ONE
*
*        T(1:I-1,I) := alpha * A(I:M,1:I-1)**H * A(I:M,I)
*
         ALPHA = -T( I, 1 )
         CALL ZGEMV( 'C', M-I+1, I-1, ALPHA, A( I, 1 ), LDA, 
     $               A( I, I ), 1, ZERO, T( 1, I ), 1 )
         A( I, I ) = AII
*
*        T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)
*
         CALL ZTRMV( 'U', 'N', 'N', I-1, T, LDT, T( 1, I ), 1 )
*
*           T(I,I) = tau(I)
*
            T( I, I ) = T( I, 1 )
            T( I, 1) = ZERO
      END DO
   
*
*     End of ZGEQRT2
*
      END
