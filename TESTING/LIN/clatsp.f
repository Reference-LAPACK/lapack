*> \brief \b CLATSP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLATSP( UPLO, N, X, ISEED )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            N
*       ..
*       .. Array Arguments ..
*       INTEGER            ISEED( * )
*       COMPLEX            X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLATSP generates a special test matrix for the complex symmetric
*> (indefinite) factorization for packed matrices.  The pivot blocks of
*> the generated matrix will be in the following order:
*>    2x2 pivot block, non diagonalizable
*>    1x1 pivot block
*>    2x2 pivot block, diagonalizable
*>    (cycle repeats)
*> A row interchange is required for each non-diagonalizable 2x2 block.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER
*>          Specifies whether the generated matrix is to be upper or
*>          lower triangular.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The dimension of the matrix to be generated.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX array, dimension (N*(N+1)/2)
*>          The generated matrix in packed storage format.  The matrix
*>          consists of 3x3 and 2x2 diagonal blocks which result in the
*>          pivot sequence given above.  The matrix outside these
*>          diagonal blocks is zero.
*> \endverbatim
*>
*> \param[in,out] ISEED
*> \verbatim
*>          ISEED is INTEGER array, dimension (4)
*>          On entry, the seed for the random number generator.  The last
*>          of the four integers must be odd.  (modified on exit)
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
*> \ingroup complex_lin
*
*  =====================================================================
      SUBROUTINE CLATSP( UPLO, N, X, ISEED )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( * )
      COMPLEX            X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            EYE
      PARAMETER          ( EYE = ( 0.0, 1.0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JJ, N5
      REAL               ALPHA, ALPHA3, BETA
      COMPLEX            A, B, C, R
*     ..
*     .. External Functions ..
      COMPLEX            CLARND
      EXTERNAL           CLARND
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Initialize constants
*
      ALPHA = ( 1.+SQRT( 17. ) ) / 8.
      BETA = ALPHA - 1. / 1000.
      ALPHA3 = ALPHA*ALPHA*ALPHA
*
*     Fill the matrix with zeros.
*
      DO 10 J = 1, N*( N+1 ) / 2
         X( J ) = 0.0
   10 CONTINUE
*
*     UPLO = 'U':  Upper triangular storage
*
      IF( UPLO.EQ.'U' ) THEN
         N5 = N / 5
         N5 = N - 5*N5 + 1
*
         JJ = N*( N+1 ) / 2
         DO 20 J = N, N5, -5
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ-2 ) = B
            JJ = JJ - J
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ-1 ) = R
            JJ = JJ - ( J-1 )
            X( JJ ) = C
            JJ = JJ - ( J-2 )
            X( JJ ) = CLARND( 2, ISEED )
            JJ = JJ - ( J-3 )
            X( JJ ) = CLARND( 2, ISEED )
            IF( ABS( X( JJ+( J-3 ) ) ).GT.ABS( X( JJ ) ) ) THEN
               X( JJ+( J-4 ) ) = 2.0*X( JJ+( J-3 ) )
            ELSE
               X( JJ+( J-4 ) ) = 2.0*X( JJ )
            END IF
            JJ = JJ - ( J-4 )
   20    CONTINUE
*
*        Clean-up for N not a multiple of 5.
*
         J = N5 - 1
         IF( J.GT.2 ) THEN
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ-2 ) = B
            JJ = JJ - J
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ-1 ) = R
            JJ = JJ - ( J-1 )
            X( JJ ) = C
            JJ = JJ - ( J-2 )
            J = J - 3
         END IF
         IF( J.GT.1 ) THEN
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ-J ) = CLARND( 2, ISEED )
            IF( ABS( X( JJ ) ).GT.ABS( X( JJ-J ) ) ) THEN
               X( JJ-1 ) = 2.0*X( JJ )
            ELSE
               X( JJ-1 ) = 2.0*X( JJ-J )
            END IF
            JJ = JJ - J - ( J-1 )
            J = J - 2
         ELSE IF( J.EQ.1 ) THEN
            X( JJ ) = CLARND( 2, ISEED )
            J = J - 1
         END IF
*
*     UPLO = 'L':  Lower triangular storage
*
      ELSE
         N5 = N / 5
         N5 = N5*5
*
         JJ = 1
         DO 30 J = 1, N5, 5
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ+2 ) = B
            JJ = JJ + ( N-J+1 )
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ+1 ) = R
            JJ = JJ + ( N-J )
            X( JJ ) = C
            JJ = JJ + ( N-J-1 )
            X( JJ ) = CLARND( 2, ISEED )
            JJ = JJ + ( N-J-2 )
            X( JJ ) = CLARND( 2, ISEED )
            IF( ABS( X( JJ-( N-J-2 ) ) ).GT.ABS( X( JJ ) ) ) THEN
               X( JJ-( N-J-2 )+1 ) = 2.0*X( JJ-( N-J-2 ) )
            ELSE
               X( JJ-( N-J-2 )+1 ) = 2.0*X( JJ )
            END IF
            JJ = JJ + ( N-J-3 )
   30    CONTINUE
*
*        Clean-up for N not a multiple of 5.
*
         J = N5 + 1
         IF( J.LT.N-1 ) THEN
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ+2 ) = B
            JJ = JJ + ( N-J+1 )
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ+1 ) = R
            JJ = JJ + ( N-J )
            X( JJ ) = C
            JJ = JJ + ( N-J-1 )
            J = J + 3
         END IF
         IF( J.LT.N ) THEN
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ+( N-J+1 ) ) = CLARND( 2, ISEED )
            IF( ABS( X( JJ ) ).GT.ABS( X( JJ+( N-J+1 ) ) ) ) THEN
               X( JJ+1 ) = 2.0*X( JJ )
            ELSE
               X( JJ+1 ) = 2.0*X( JJ+( N-J+1 ) )
            END IF
            JJ = JJ + ( N-J+1 ) + ( N-J )
            J = J + 2
         ELSE IF( J.EQ.N ) THEN
            X( JJ ) = CLARND( 2, ISEED )
            JJ = JJ + ( N-J+1 )
            J = J + 1
         END IF
      END IF
*
      RETURN
*
*     End of CLATSP
*
      END
