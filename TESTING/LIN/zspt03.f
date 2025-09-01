*> \brief \b ZSPT03
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSPT03( UPLO, N, A, AINV, WORK, LDW, RWORK, RCOND,
*                          RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDW, N
*       DOUBLE PRECISION   RCOND, RESID
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( * ), AINV( * ), WORK( LDW, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZSPT03 computes the residual for a complex symmetric packed matrix
*> times its inverse:
*>    norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS ),
*> where EPS is the machine epsilon.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          complex symmetric matrix A is stored:
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
*>          A is COMPLEX*16 array, dimension (N*(N+1)/2)
*>          The original complex symmetric matrix A, stored as a packed
*>          triangular matrix.
*> \endverbatim
*>
*> \param[in] AINV
*> \verbatim
*>          AINV is COMPLEX*16 array, dimension (N*(N+1)/2)
*>          The (symmetric) inverse of the matrix A, stored as a packed
*>          triangular matrix.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (LDW,N)
*> \endverbatim
*>
*> \param[in] LDW
*> \verbatim
*>          LDW is INTEGER
*>          The leading dimension of the array WORK.  LDW >= max(1,N).
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] RCOND
*> \verbatim
*>          RCOND is DOUBLE PRECISION
*>          The reciprocal of the condition number of A, computed as
*>          ( 1/norm(A) ) / norm(AINV).
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is DOUBLE PRECISION
*>          norm(I - A*AINV) / ( N * norm(A) * norm(AINV) * EPS )
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
      SUBROUTINE ZSPT03( UPLO, N, A, AINV, WORK, LDW, RWORK, RCOND,
     $                   RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDW, N
      DOUBLE PRECISION   RCOND, RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( * ), AINV( * ), WORK( LDW, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOL, J, JCOL, K, KCOL, NALL
      DOUBLE PRECISION   AINVNM, ANORM, EPS
      COMPLEX*16         T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANSP
      COMPLEX*16         ZDOTU
      EXTERNAL           LSAME, DLAMCH, ZLANGE, ZLANSP, ZDOTU
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RCOND = ONE
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANSP( '1', UPLO, N, A, RWORK )
      AINVNM = ZLANSP( '1', UPLO, N, AINV, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
*     Case where both A and AINV are upper triangular:
*     Each element of - A * AINV is computed by taking the dot product
*     of a row of A with a column of AINV.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 70 I = 1, N
            ICOL = ( ( I-1 )*I ) / 2 + 1
*
*           Code when J <= I
*
            DO 30 J = 1, I
               JCOL = ( ( J-1 )*J ) / 2 + 1
               T = ZDOTU( J, A( ICOL ), 1, AINV( JCOL ), 1 )
               JCOL = JCOL + 2*J - 1
               KCOL = ICOL - 1
               DO 10 K = J + 1, I
                  T = T + A( KCOL+K )*AINV( JCOL )
                  JCOL = JCOL + K
   10          CONTINUE
               KCOL = KCOL + 2*I
               DO 20 K = I + 1, N
                  T = T + A( KCOL )*AINV( JCOL )
                  KCOL = KCOL + K
                  JCOL = JCOL + K
   20          CONTINUE
               WORK( I, J ) = -T
   30       CONTINUE
*
*           Code when J > I
*
            DO 60 J = I + 1, N
               JCOL = ( ( J-1 )*J ) / 2 + 1
               T = ZDOTU( I, A( ICOL ), 1, AINV( JCOL ), 1 )
               JCOL = JCOL - 1
               KCOL = ICOL + 2*I - 1
               DO 40 K = I + 1, J
                  T = T + A( KCOL )*AINV( JCOL+K )
                  KCOL = KCOL + K
   40          CONTINUE
               JCOL = JCOL + 2*J
               DO 50 K = J + 1, N
                  T = T + A( KCOL )*AINV( JCOL )
                  KCOL = KCOL + K
                  JCOL = JCOL + K
   50          CONTINUE
               WORK( I, J ) = -T
   60       CONTINUE
   70    CONTINUE
      ELSE
*
*        Case where both A and AINV are lower triangular
*
         NALL = ( N*( N+1 ) ) / 2
         DO 140 I = 1, N
*
*           Code when J <= I
*
            ICOL = NALL - ( ( N-I+1 )*( N-I+2 ) ) / 2 + 1
            DO 100 J = 1, I
               JCOL = NALL - ( ( N-J )*( N-J+1 ) ) / 2 - ( N-I )
               T = ZDOTU( N-I+1, A( ICOL ), 1, AINV( JCOL ), 1 )
               KCOL = I
               JCOL = J
               DO 80 K = 1, J - 1
                  T = T + A( KCOL )*AINV( JCOL )
                  JCOL = JCOL + N - K
                  KCOL = KCOL + N - K
   80          CONTINUE
               JCOL = JCOL - J
               DO 90 K = J, I - 1
                  T = T + A( KCOL )*AINV( JCOL+K )
                  KCOL = KCOL + N - K
   90          CONTINUE
               WORK( I, J ) = -T
  100       CONTINUE
*
*           Code when J > I
*
            ICOL = NALL - ( ( N-I )*( N-I+1 ) ) / 2
            DO 130 J = I + 1, N
               JCOL = NALL - ( ( N-J+1 )*( N-J+2 ) ) / 2 + 1
               T = ZDOTU( N-J+1, A( ICOL-N+J ), 1, AINV( JCOL ), 1 )
               KCOL = I
               JCOL = J
               DO 110 K = 1, I - 1
                  T = T + A( KCOL )*AINV( JCOL )
                  JCOL = JCOL + N - K
                  KCOL = KCOL + N - K
  110          CONTINUE
               KCOL = KCOL - I
               DO 120 K = I, J - 1
                  T = T + A( KCOL+K )*AINV( JCOL )
                  JCOL = JCOL + N - K
  120          CONTINUE
               WORK( I, J ) = -T
  130       CONTINUE
  140    CONTINUE
      END IF
*
*     Add the identity matrix to WORK .
*
      DO 150 I = 1, N
         WORK( I, I ) = WORK( I, I ) + ONE
  150 CONTINUE
*
*     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = ZLANGE( '1', N, N, WORK, LDW, RWORK )
*
      RESID = ( ( RESID*RCOND ) / EPS ) / DBLE( N )
*
      RETURN
*
*     End of ZSPT03
*
      END
