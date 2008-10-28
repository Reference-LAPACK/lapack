      SUBROUTINE CLATRZ( M, N, L, A, LDA, TAU, WORK )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            L, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CLATRZ factors the M-by-(M+L) complex upper trapezoidal matrix
*  [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z by means
*  of unitary transformations, where  Z is an (M+L)-by-(M+L) unitary
*  matrix and, R and A1 are M-by-M upper triangular matrices.
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
*  L       (input) INTEGER
*          The number of columns of the matrix A containing the
*          meaningful part of the Householder vectors. N-M >= L >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the leading M-by-N upper trapezoidal part of the
*          array A must contain the matrix to be factorized.
*          On exit, the leading M-by-M upper triangular part of A
*          contains the upper triangular matrix R, and elements N-L+1 to
*          N of the first M rows of A, with the array TAU, represent the
*          unitary matrix Z as a product of M elementary reflectors.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX array, dimension (M)
*          The scalar factors of the elementary reflectors.
*
*  WORK    (workspace) COMPLEX array, dimension (M)
*
*  Further Details
*  ===============
*
*  Based on contributions by
*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
*  The factorization is obtained by Householder's method.  The kth
*  transformation matrix, Z( k ), which is used to introduce zeros into
*  the ( m - k + 1 )th row of A, is given in the form
*
*     Z( k ) = ( I     0   ),
*              ( 0  T( k ) )
*
*  where
*
*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),
*                                                 (   0    )
*                                                 ( z( k ) )
*
*  tau is a scalar and z( k ) is an l element vector. tau and z( k )
*  are chosen to annihilate the elements of the kth row of A2.
*
*  The scalar tau is returned in the kth element of TAU and the vector
*  u( k ) in the kth row of A2, such that the elements of z( k ) are
*  in  a( k, l + 1 ), ..., a( k, n ). The elements of R are returned in
*  the upper triangular part of A1.
*
*  Z is given by
*
*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX            ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACGV, CLARFG, CLARZ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.EQ.0 ) THEN
         RETURN
      ELSE IF( M.EQ.N ) THEN
         DO 10 I = 1, N
            TAU( I ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
      DO 20 I = M, 1, -1
*
*        Generate elementary reflector H(i) to annihilate
*        [ A(i,i) A(i,n-l+1:n) ]
*
         CALL CLACGV( L, A( I, N-L+1 ), LDA )
         ALPHA = CONJG( A( I, I ) )
         CALL CLARFG( L+1, ALPHA, A( I, N-L+1 ), LDA, TAU( I ) )
         TAU( I ) = CONJG( TAU( I ) )
*
*        Apply H(i) to A(1:i-1,i:n) from the right
*
         CALL CLARZ( 'Right', I-1, N-I+1, L, A( I, N-L+1 ), LDA,
     $               CONJG( TAU( I ) ), A( 1, I ), LDA, WORK )
         A( I, I ) = CONJG( ALPHA )
*
   20 CONTINUE
*
      RETURN
*
*     End of CLATRZ
*
      END
