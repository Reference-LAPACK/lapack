*> \brief \b ZHET01_ROOK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZHET01_ROOK( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC,
*                               RWORK, RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDA, LDAFAC, LDC, N
*       DOUBLE PRECISION   RESID
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZHET01_ROOK reconstructs a complex Hermitian indefinite matrix A from its
*> block L*D*L' or U*D*U' factorization and computes the residual
*>    norm( C - A ) / ( N * norm(A) * EPS ),
*> where C is the reconstructed matrix, EPS is the machine epsilon,
*> L' is the transpose of L, and U' is the transpose of U.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          complex Hermitian matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows and columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The original complex Hermitian matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N)
*> \endverbatim
*>
*> \param[in] AFAC
*> \verbatim
*>          AFAC is COMPLEX*16 array, dimension (LDAFAC,N)
*>          The factored form of the matrix A.  AFAC contains the block
*>          diagonal matrix D and the multipliers used to obtain the
*>          factor L or U from the block L*D*L' or U*D*U' factorization
*>          as computed by CSYTRF_ROOK.
*> \endverbatim
*>
*> \param[in] LDAFAC
*> \verbatim
*>          LDAFAC is INTEGER
*>          The leading dimension of the array AFAC.  LDAFAC >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices from CSYTRF_ROOK.
*> \endverbatim
*>
*> \param[out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C.  LDC >= max(1,N).
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is DOUBLE PRECISION
*>          If UPLO = 'L', norm(L*D*L' - A) / ( N * norm(A) * EPS )
*>          If UPLO = 'U', norm(U*D*U' - A) / ( N * norm(A) * EPS )
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
*> \ingroup complex16_lin
*
*  =====================================================================
      SUBROUTINE ZHET01_ROOK( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C,
     $                        LDC, RWORK, RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDAFAC, LDC, N
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   ANORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   ZLANHE, DLAMCH
      EXTERNAL           LSAME, ZLANHE, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASET, ZLAVHE_ROOK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DIMAG, DBLE
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Determine EPS and the norm of A.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
*
*     Check the imaginary parts of the diagonal elements and return with
*     an error code if any are nonzero.
*
      DO 10 J = 1, N
         IF( DIMAG( AFAC( J, J ) ).NE.ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
   10 CONTINUE
*
*     Initialize C to the identity matrix.
*
      CALL ZLASET( 'Full', N, N, CZERO, CONE, C, LDC )
*
*     Call ZLAVHE_ROOK to form the product D * U' (or D * L' ).
*
      CALL ZLAVHE_ROOK( UPLO, 'Conjugate', 'Non-unit', N, N, AFAC,
     $                  LDAFAC, IPIV, C, LDC, INFO )
*
*     Call ZLAVHE_ROOK again to multiply by U (or L ).
*
      CALL ZLAVHE_ROOK( UPLO, 'No transpose', 'Unit', N, N, AFAC,
     $                  LDAFAC, IPIV, C, LDC, INFO )
*
*     Compute the difference  C - A .
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 30 J = 1, N
            DO 20 I = 1, J - 1
               C( I, J ) = C( I, J ) - A( I, J )
   20       CONTINUE
            C( J, J ) = C( J, J ) - DBLE( A( J, J ) )
   30    CONTINUE
      ELSE
         DO 50 J = 1, N
            C( J, J ) = C( J, J ) - DBLE( A( J, J ) )
            DO 40 I = J + 1, N
               C( I, J ) = C( I, J ) - A( I, J )
   40       CONTINUE
   50    CONTINUE
      END IF
*
*     Compute norm( C - A ) / ( N * norm(A) * EPS )
*
      RESID = ZLANHE( '1', UPLO, N, C, LDC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID/DBLE( N ) )/ANORM ) / EPS
      END IF
*
      RETURN
*
*     End of ZHET01_ROOK
*
      END
