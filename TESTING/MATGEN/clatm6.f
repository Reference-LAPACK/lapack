      SUBROUTINE CLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA,
     $                   BETA, WX, WY, S, DIF )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     June 2010
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDX, LDY, N, TYPE
      COMPLEX            ALPHA, BETA, WX, WY
*     ..
*     .. Array Arguments ..
      REAL               DIF( * ), S( * )
      COMPLEX            A( LDA, * ), B( LDA, * ), X( LDX, * ),
     $                   Y( LDY, * )
*     ..
*
*  Purpose
*  =======
*
*  CLATM6 generates test matrices for the generalized eigenvalue
*  problem, their corresponding right and left eigenvector matrices,
*  and also reciprocal condition numbers for all eigenvalues and
*  the reciprocal condition numbers of eigenvectors corresponding to
*  the 1th and 5th eigenvalues.
*
*  Test Matrices
*  =============
*
*  Two kinds of test matrix pairs
*           (A, B) = inverse(YH) * (Da, Db) * inverse(X)
*  are used in the tests:
*
*  Type 1:
*     Da = 1+a   0    0    0    0    Db = 1   0   0   0   0
*           0   2+a   0    0    0         0   1   0   0   0
*           0    0   3+a   0    0         0   0   1   0   0
*           0    0    0   4+a   0         0   0   0   1   0
*           0    0    0    0   5+a ,      0   0   0   0   1
*  and Type 2:
*     Da = 1+i   0    0       0       0    Db = 1   0   0   0   0
*           0   1-i   0       0       0         0   1   0   0   0
*           0    0    1       0       0         0   0   1   0   0
*           0    0    0 (1+a)+(1+b)i  0         0   0   0   1   0
*           0    0    0       0 (1+a)-(1+b)i,   0   0   0   0   1 .
*
*  In both cases the same inverse(YH) and inverse(X) are used to compute
*  (A, B), giving the exact eigenvectors to (A,B) as (YH, X):
*
*  YH:  =  1    0   -y    y   -y    X =  1   0  -x  -x   x
*          0    1   -y    y   -y         0   1   x  -x  -x
*          0    0    1    0    0         0   0   1   0   0
*          0    0    0    1    0         0   0   0   1   0
*          0    0    0    0    1,        0   0   0   0   1 , where
*
*  a, b, x and y will have all values independently of each other.
*
*  Arguments
*  =========
*
*  TYPE    (input) INTEGER
*          Specifies the problem type (see futher details).
*
*  N       (input) INTEGER
*          Size of the matrices A and B.
*
*  A       (output) COMPLEX array, dimension (LDA, N).
*          On exit A N-by-N is initialized according to TYPE.
*
*  LDA     (input) INTEGER
*          The leading dimension of A and of B.
*
*  B       (output) COMPLEX array, dimension (LDA, N).
*          On exit B N-by-N is initialized according to TYPE.
*
*  X       (output) COMPLEX array, dimension (LDX, N).
*          On exit X is the N-by-N matrix of right eigenvectors.
*
*  LDX     (input) INTEGER
*          The leading dimension of X.
*
*  Y       (output) COMPLEX array, dimension (LDY, N).
*          On exit Y is the N-by-N matrix of left eigenvectors.
*
*  LDY     (input) INTEGER
*          The leading dimension of Y.
*
*  ALPHA   (input) COMPLEX
*
*  BETA    (input) COMPLEX
*          Weighting constants for matrix A.
*
*  WX      (input) COMPLEX
*          Constant for right eigenvector matrix.
*
*  WY      (input) COMPLEX
*          Constant for left eigenvector matrix.
*
*  S       (output) REAL array, dimension (N)
*          S(i) is the reciprocal condition number for eigenvalue i.
*
*  DIF     (output) REAL array, dimension (N)
*          DIF(i) is the reciprocal condition number for eigenvector i.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               RONE, TWO, THREE
      PARAMETER          ( RONE = 1.0E+0, TWO = 2.0E+0, THREE = 3.0E+0 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),
     $                   ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
*     ..
*     .. Local Arrays ..
      REAL               RWORK( 50 )
      COMPLEX            WORK( 26 ), Z( 8, 8 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CABS, CMPLX, CONJG, REAL, SQRT
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGESVD, CLACPY, CLAKF2
*     ..
*     .. Executable Statements ..
*
*     Generate test problem ...
*     (Da, Db) ...
*
      DO 20 I = 1, N
         DO 10 J = 1, N
*
            IF( I.EQ.J ) THEN
               A( I, I ) = CMPLX( I ) + ALPHA
               B( I, I ) = ONE
            ELSE
               A( I, J ) = ZERO
               B( I, J ) = ZERO
            END IF
*
   10    CONTINUE
   20 CONTINUE
      IF( TYPE.EQ.2 ) THEN
         A( 1, 1 ) = CMPLX( RONE, RONE )
         A( 2, 2 ) = CONJG( A( 1, 1 ) )
         A( 3, 3 ) = ONE
         A( 4, 4 ) = CMPLX( REAL( ONE+ALPHA ), REAL( ONE+BETA ) )
         A( 5, 5 ) = CONJG( A( 4, 4 ) )
      END IF
*
*     Form X and Y
*
      CALL CLACPY( 'F', N, N, B, LDA, Y, LDY )
      Y( 3, 1 ) = -CONJG( WY )
      Y( 4, 1 ) = CONJG( WY )
      Y( 5, 1 ) = -CONJG( WY )
      Y( 3, 2 ) = -CONJG( WY )
      Y( 4, 2 ) = CONJG( WY )
      Y( 5, 2 ) = -CONJG( WY )
*
      CALL CLACPY( 'F', N, N, B, LDA, X, LDX )
      X( 1, 3 ) = -WX
      X( 1, 4 ) = -WX
      X( 1, 5 ) = WX
      X( 2, 3 ) = WX
      X( 2, 4 ) = -WX
      X( 2, 5 ) = -WX
*
*     Form (A, B)
*
      B( 1, 3 ) = WX + WY
      B( 2, 3 ) = -WX + WY
      B( 1, 4 ) = WX - WY
      B( 2, 4 ) = WX - WY
      B( 1, 5 ) = -WX + WY
      B( 2, 5 ) = WX + WY
      A( 1, 3 ) = WX*A( 1, 1 ) + WY*A( 3, 3 )
      A( 2, 3 ) = -WX*A( 2, 2 ) + WY*A( 3, 3 )
      A( 1, 4 ) = WX*A( 1, 1 ) - WY*A( 4, 4 )
      A( 2, 4 ) = WX*A( 2, 2 ) - WY*A( 4, 4 )
      A( 1, 5 ) = -WX*A( 1, 1 ) + WY*A( 5, 5 )
      A( 2, 5 ) = WX*A( 2, 2 ) + WY*A( 5, 5 )
*
*     Compute condition numbers
*
      S( 1 ) = RONE / SQRT( ( RONE+THREE*CABS( WY )*CABS( WY ) ) /
     $         ( RONE+CABS( A( 1, 1 ) )*CABS( A( 1, 1 ) ) ) )
      S( 2 ) = RONE / SQRT( ( RONE+THREE*CABS( WY )*CABS( WY ) ) /
     $         ( RONE+CABS( A( 2, 2 ) )*CABS( A( 2, 2 ) ) ) )
      S( 3 ) = RONE / SQRT( ( RONE+TWO*CABS( WX )*CABS( WX ) ) /
     $         ( RONE+CABS( A( 3, 3 ) )*CABS( A( 3, 3 ) ) ) )
      S( 4 ) = RONE / SQRT( ( RONE+TWO*CABS( WX )*CABS( WX ) ) /
     $         ( RONE+CABS( A( 4, 4 ) )*CABS( A( 4, 4 ) ) ) )
      S( 5 ) = RONE / SQRT( ( RONE+TWO*CABS( WX )*CABS( WX ) ) /
     $         ( RONE+CABS( A( 5, 5 ) )*CABS( A( 5, 5 ) ) ) )
*
      CALL CLAKF2( 1, 4, A, LDA, A( 2, 2 ), B, B( 2, 2 ), Z, 8 )
      CALL CGESVD( 'N', 'N', 8, 8, Z, 8, RWORK, WORK, 1, WORK( 2 ), 1,
     $             WORK( 3 ), 24, RWORK( 9 ), INFO )
      DIF( 1 ) = RWORK( 8 )
*
      CALL CLAKF2( 4, 1, A, LDA, A( 5, 5 ), B, B( 5, 5 ), Z, 8 )
      CALL CGESVD( 'N', 'N', 8, 8, Z, 8, RWORK, WORK, 1, WORK( 2 ), 1,
     $             WORK( 3 ), 24, RWORK( 9 ), INFO )
      DIF( 5 ) = RWORK( 8 )
*
      RETURN
*
*     End of CLATM6
*
      END
