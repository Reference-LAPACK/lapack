      SUBROUTINE ZLAKF2( M, N, A, LDA, B, D, E, Z, LDZ )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDZ, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDA, * ), D( LDA, * ),
     $                   E( LDA, * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  Form the 2*M*N by 2*M*N matrix
*
*         Z = [ kron(In, A)  -kron(B', Im) ]
*             [ kron(In, D)  -kron(E', Im) ],
*
*  where In is the identity matrix of size n and X' is the transpose
*  of X. kron(X, Y) is the Kronecker product between the matrices X
*  and Y.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          Size of matrix, must be >= 1.
*
*  N       (input) INTEGER
*          Size of matrix, must be >= 1.
*
*  A       (input) COMPLEX*16, dimension ( LDA, M )
*          The matrix A in the output matrix Z.
*
*  LDA     (input) INTEGER
*          The leading dimension of A, B, D, and E. ( LDA >= M+N )
*
*  B       (input) COMPLEX*16, dimension ( LDA, N )
*  D       (input) COMPLEX*16, dimension ( LDA, M )
*  E       (input) COMPLEX*16, dimension ( LDA, N )
*          The matrices used in forming the output matrix Z.
*
*  Z       (output) COMPLEX*16, dimension ( LDZ, 2*M*N )
*          The resultant Kronecker M*N*2 by M*N*2 matrix (see above.)
*
*  LDZ     (input) INTEGER
*          The leading dimension of Z. ( LDZ >= 2*M*N )
*
*  ====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IK, J, JK, L, MN, MN2
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASET
*     ..
*     .. Executable Statements ..
*
*     Initialize Z
*
      MN = M*N
      MN2 = 2*MN
      CALL ZLASET( 'Full', MN2, MN2, ZERO, ZERO, Z, LDZ )
*
      IK = 1
      DO 50 L = 1, N
*
*        form kron(In, A)
*
         DO 20 I = 1, M
            DO 10 J = 1, M
               Z( IK+I-1, IK+J-1 ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
*
*        form kron(In, D)
*
         DO 40 I = 1, M
            DO 30 J = 1, M
               Z( IK+MN+I-1, IK+J-1 ) = D( I, J )
   30       CONTINUE
   40    CONTINUE
*
         IK = IK + M
   50 CONTINUE
*
      IK = 1
      DO 90 L = 1, N
         JK = MN + 1
*
         DO 80 J = 1, N
*
*           form -kron(B', Im)
*
            DO 60 I = 1, M
               Z( IK+I-1, JK+I-1 ) = -B( J, L )
   60       CONTINUE
*
*           form -kron(E', Im)
*
            DO 70 I = 1, M
               Z( IK+MN+I-1, JK+I-1 ) = -E( J, L )
   70       CONTINUE
*
            JK = JK + M
   80    CONTINUE
*
         IK = IK + M
   90 CONTINUE
*
      RETURN
*
*     End of ZLAKF2
*
      END
