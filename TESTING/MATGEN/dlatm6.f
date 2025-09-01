*> \brief \b DLATM6
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA,
*                          BETA, WX, WY, S, DIF )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDX, LDY, N, TYPE
*       DOUBLE PRECISION   ALPHA, BETA, WX, WY
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDA, * ), DIF( * ), S( * ),
*      $                   X( LDX, * ), Y( LDY, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLATM6 generates test matrices for the generalized eigenvalue
*> problem, their corresponding right and left eigenvector matrices,
*> and also reciprocal condition numbers for all eigenvalues and
*> the reciprocal condition numbers of eigenvectors corresponding to
*> the 1th and 5th eigenvalues.
*>
*> Test Matrices
*> =============
*>
*> Two kinds of test matrix pairs
*>
*>       (A, B) = inverse(YH) * (Da, Db) * inverse(X)
*>
*> are used in the tests:
*>
*> Type 1:
*>    Da = 1+a   0    0    0    0    Db = 1   0   0   0   0
*>          0   2+a   0    0    0         0   1   0   0   0
*>          0    0   3+a   0    0         0   0   1   0   0
*>          0    0    0   4+a   0         0   0   0   1   0
*>          0    0    0    0   5+a ,      0   0   0   0   1 , and
*>
*> Type 2:
*>    Da =  1   -1    0    0    0    Db = 1   0   0   0   0
*>          1    1    0    0    0         0   1   0   0   0
*>          0    0    1    0    0         0   0   1   0   0
*>          0    0    0   1+a  1+b        0   0   0   1   0
*>          0    0    0  -1-b  1+a ,      0   0   0   0   1 .
*>
*> In both cases the same inverse(YH) and inverse(X) are used to compute
*> (A, B), giving the exact eigenvectors to (A,B) as (YH, X):
*>
*> YH:  =  1    0   -y    y   -y    X =  1   0  -x  -x   x
*>         0    1   -y    y   -y         0   1   x  -x  -x
*>         0    0    1    0    0         0   0   1   0   0
*>         0    0    0    1    0         0   0   0   1   0
*>         0    0    0    0    1,        0   0   0   0   1 ,
*>
*> where a, b, x and y will have all values independently of each other.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TYPE
*> \verbatim
*>          TYPE is INTEGER
*>          Specifies the problem type (see further details).
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          Size of the matrices A and B.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N).
*>          On exit A N-by-N is initialized according to TYPE.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A and of B.
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDA, N).
*>          On exit B N-by-N is initialized according to TYPE.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension (LDX, N).
*>          On exit X is the N-by-N matrix of right eigenvectors.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of X.
*> \endverbatim
*>
*> \param[out] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array, dimension (LDY, N).
*>          On exit Y is the N-by-N matrix of left eigenvectors.
*> \endverbatim
*>
*> \param[in] LDY
*> \verbatim
*>          LDY is INTEGER
*>          The leading dimension of Y.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION
*>
*>          Weighting constants for matrix A.
*> \endverbatim
*>
*> \param[in] WX
*> \verbatim
*>          WX is DOUBLE PRECISION
*>          Constant for right eigenvector matrix.
*> \endverbatim
*>
*> \param[in] WY
*> \verbatim
*>          WY is DOUBLE PRECISION
*>          Constant for left eigenvector matrix.
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is DOUBLE PRECISION array, dimension (N)
*>          S(i) is the reciprocal condition number for eigenvalue i.
*> \endverbatim
*>
*> \param[out] DIF
*> \verbatim
*>          DIF is DOUBLE PRECISION array, dimension (N)
*>          DIF(i) is the reciprocal condition number for eigenvector i.
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
*> \ingroup double_matgen
*
*  =====================================================================
      SUBROUTINE DLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA,
     $                   BETA, WX, WY, S, DIF )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDX, LDY, N, TYPE
      DOUBLE PRECISION   ALPHA, BETA, WX, WY
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDA, * ), DIF( * ), S( * ),
     $                   X( LDX, * ), Y( LDY, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                   THREE = 3.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   WORK( 100 ), Z( 12, 12 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, SQRT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGESVD, DLACPY, DLAKF2
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
               A( I, I ) = DBLE( I ) + ALPHA
               B( I, I ) = ONE
            ELSE
               A( I, J ) = ZERO
               B( I, J ) = ZERO
            END IF
*
   10    CONTINUE
   20 CONTINUE
*
*     Form X and Y
*
      CALL DLACPY( 'F', N, N, B, LDA, Y, LDY )
      Y( 3, 1 ) = -WY
      Y( 4, 1 ) = WY
      Y( 5, 1 ) = -WY
      Y( 3, 2 ) = -WY
      Y( 4, 2 ) = WY
      Y( 5, 2 ) = -WY
*
      CALL DLACPY( 'F', N, N, B, LDA, X, LDX )
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
      IF( TYPE.EQ.1 ) THEN
         A( 1, 3 ) = WX*A( 1, 1 ) + WY*A( 3, 3 )
         A( 2, 3 ) = -WX*A( 2, 2 ) + WY*A( 3, 3 )
         A( 1, 4 ) = WX*A( 1, 1 ) - WY*A( 4, 4 )
         A( 2, 4 ) = WX*A( 2, 2 ) - WY*A( 4, 4 )
         A( 1, 5 ) = -WX*A( 1, 1 ) + WY*A( 5, 5 )
         A( 2, 5 ) = WX*A( 2, 2 ) + WY*A( 5, 5 )
      ELSE IF( TYPE.EQ.2 ) THEN
         A( 1, 3 ) = TWO*WX + WY
         A( 2, 3 ) = WY
         A( 1, 4 ) = -WY*( TWO+ALPHA+BETA )
         A( 2, 4 ) = TWO*WX - WY*( TWO+ALPHA+BETA )
         A( 1, 5 ) = -TWO*WX + WY*( ALPHA-BETA )
         A( 2, 5 ) = WY*( ALPHA-BETA )
         A( 1, 1 ) = ONE
         A( 1, 2 ) = -ONE
         A( 2, 1 ) = ONE
         A( 2, 2 ) = A( 1, 1 )
         A( 3, 3 ) = ONE
         A( 4, 4 ) = ONE + ALPHA
         A( 4, 5 ) = ONE + BETA
         A( 5, 4 ) = -A( 4, 5 )
         A( 5, 5 ) = A( 4, 4 )
      END IF
*
*     Compute condition numbers
*
      IF( TYPE.EQ.1 ) THEN
*
         S( 1 ) = ONE / SQRT( ( ONE+THREE*WY*WY ) /
     $            ( ONE+A( 1, 1 )*A( 1, 1 ) ) )
         S( 2 ) = ONE / SQRT( ( ONE+THREE*WY*WY ) /
     $            ( ONE+A( 2, 2 )*A( 2, 2 ) ) )
         S( 3 ) = ONE / SQRT( ( ONE+TWO*WX*WX ) /
     $            ( ONE+A( 3, 3 )*A( 3, 3 ) ) )
         S( 4 ) = ONE / SQRT( ( ONE+TWO*WX*WX ) /
     $            ( ONE+A( 4, 4 )*A( 4, 4 ) ) )
         S( 5 ) = ONE / SQRT( ( ONE+TWO*WX*WX ) /
     $            ( ONE+A( 5, 5 )*A( 5, 5 ) ) )
*
         CALL DLAKF2( 1, 4, A, LDA, A( 2, 2 ), B, B( 2, 2 ), Z, 12 )
         CALL DGESVD( 'N', 'N', 8, 8, Z, 12, WORK, WORK( 9 ), 1,
     $                WORK( 10 ), 1, WORK( 11 ), 40, INFO )
         DIF( 1 ) = WORK( 8 )
*
         CALL DLAKF2( 4, 1, A, LDA, A( 5, 5 ), B, B( 5, 5 ), Z, 12 )
         CALL DGESVD( 'N', 'N', 8, 8, Z, 12, WORK, WORK( 9 ), 1,
     $                WORK( 10 ), 1, WORK( 11 ), 40, INFO )
         DIF( 5 ) = WORK( 8 )
*
      ELSE IF( TYPE.EQ.2 ) THEN
*
         S( 1 ) = ONE / SQRT( ONE / THREE+WY*WY )
         S( 2 ) = S( 1 )
         S( 3 ) = ONE / SQRT( ONE / TWO+WX*WX )
         S( 4 ) = ONE / SQRT( ( ONE+TWO*WX*WX ) /
     $            ( ONE+( ONE+ALPHA )*( ONE+ALPHA )+( ONE+BETA )*( ONE+
     $            BETA ) ) )
         S( 5 ) = S( 4 )
*
         CALL DLAKF2( 2, 3, A, LDA, A( 3, 3 ), B, B( 3, 3 ), Z, 12 )
         CALL DGESVD( 'N', 'N', 12, 12, Z, 12, WORK, WORK( 13 ), 1,
     $                WORK( 14 ), 1, WORK( 15 ), 60, INFO )
         DIF( 1 ) = WORK( 12 )
*
         CALL DLAKF2( 3, 2, A, LDA, A( 4, 4 ), B, B( 4, 4 ), Z, 12 )
         CALL DGESVD( 'N', 'N', 12, 12, Z, 12, WORK, WORK( 13 ), 1,
     $                WORK( 14 ), 1, WORK( 15 ), 60, INFO )
         DIF( 5 ) = WORK( 12 )
*
      END IF
*
      RETURN
*
*     End of DLATM6
*
      END
