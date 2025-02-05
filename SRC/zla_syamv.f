*> \brief \b ZLA_SYAMV computes a matrix-vector product using a symmetric indefinite matrix to calculate error bounds.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZLA_SYAMV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_syamv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_syamv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_syamv.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y,
*                             INCY )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   ALPHA, BETA
*       INTEGER            INCX, INCY, LDA, N
*       INTEGER            UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), X( * )
*       DOUBLE PRECISION   Y( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLA_SYAMV  performs the matrix-vector operation
*>
*>         y := alpha*abs(A)*abs(x) + beta*abs(y),
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> n by n symmetric matrix.
*>
*> This function is primarily used in calculating error bounds.
*> To protect against underflow during evaluation, components in
*> the resulting vector are perturbed away from zero by (N+1)
*> times the underflow threshold.  To prevent unnecessarily large
*> errors for block-structure embedded in general matrices,
*> "symbolically" zero components are not perturbed.  A zero
*> entry is considered "symbolic" if all multiplications involved
*> in computing that entry have at least one zero multiplicand.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is INTEGER
*>           On entry, UPLO specifies whether the upper or lower
*>           triangular part of the array A is to be referenced as
*>           follows:
*>
*>              UPLO = BLAS_UPPER   Only the upper triangular part of A
*>                                  is to be referenced.
*>
*>              UPLO = BLAS_LOWER   Only the lower triangular part of A
*>                                  is to be referenced.
*>
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION .
*>           On entry, ALPHA specifies the scalar alpha.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, n ).
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, n ).
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) )
*>           Before entry, the incremented array X must contain the
*>           vector x.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION .
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array, dimension
*>           ( 1 + ( n - 1 )*abs( INCY ) )
*>           Before entry with BETA non-zero, the incremented array Y
*>           must contain the vector y. On exit, Y is overwritten by the
*>           updated vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*>           Unchanged on exit.
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
*> \ingroup la_heamv
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*>  -- Modified for the absolute-value product, April 2006
*>     Jason Riedy, UC Berkeley
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y,
     $                      INCY )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      INTEGER            UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * )
      DOUBLE PRECISION   Y( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            SYMB_ZERO
      DOUBLE PRECISION   TEMP, SAFE1
      INTEGER            I, INFO, IY, J, JX, KX, KY
      COMPLEX*16         ZDUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DLAMCH
      DOUBLE PRECISION   DLAMCH
*     ..
*     .. External Functions ..
      EXTERNAL           ILAUPLO
      INTEGER            ILAUPLO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, ABS, SIGN, REAL, DIMAG
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE ( ZDUM ) ) + ABS( DIMAG ( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( UPLO.NE.ILAUPLO( 'U' ) .AND.
     $         UPLO.NE.ILAUPLO( 'L' ) )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZLA_SYAMV', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
*
*     Set SAFE1 essentially to be the underflow threshold times the
*     number of additions in each row.
*
      SAFE1 = DLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1
*
*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y).
*
*     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to
*     the inexact flag.  Still doesn't help change the iteration order
*     to per-column.
*
      IY = KY
      IF ( INCX.EQ.1 ) THEN
         IF ( UPLO .EQ. ILAUPLO( 'U' ) ) THEN
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
                  DO J = I+1, N
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO )
     $              Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         ELSE
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
                  DO J = I+1, N
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO )
     $              Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         END IF
      ELSE
         IF ( UPLO .EQ. ILAUPLO( 'U' ) ) THEN
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               JX = KX
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
                  DO J = I+1, N
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO )
     $              Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         ELSE
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               JX = KX
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
                  DO J = I+1, N
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO )
     $              Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         END IF

      END IF
*
      RETURN
*
*     End of ZLA_SYAMV
*
      END
