*> \brief \b ZBDT01
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK,
*                          RWORK, RESID )
*
*       .. Scalar Arguments ..
*       INTEGER            KD, LDA, LDPT, LDQ, M, N
*       DOUBLE PRECISION   RESID
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   D( * ), E( * ), RWORK( * )
*       COMPLEX*16         A( LDA, * ), PT( LDPT, * ), Q( LDQ, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZBDT01 reconstructs a general matrix A from its bidiagonal form
*>    A = Q * B * P**H
*> where Q (m by min(m,n)) and P**H (min(m,n) by n) are unitary
*> matrices and B is bidiagonal.
*>
*> The test ratio to test the reduction is
*>    RESID = norm(A - Q * B * P**H) / ( n * norm(A) * EPS )
*> where EPS is the machine precision.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrices A and Q.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and P**H.
*> \endverbatim
*>
*> \param[in] KD
*> \verbatim
*>          KD is INTEGER
*>          If KD = 0, B is diagonal and the array E is not referenced.
*>          If KD = 1, the reduction was performed by xGEBRD; B is upper
*>          bidiagonal if M >= N, and lower bidiagonal if M < N.
*>          If KD = -1, the reduction was performed by xGBBRD; B is
*>          always upper bidiagonal.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDQ,N)
*>          The m by min(m,n) unitary matrix Q in the reduction
*>          A = Q * B * P**H.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.  LDQ >= max(1,M).
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (min(M,N))
*>          The diagonal elements of the bidiagonal matrix B.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (min(M,N)-1)
*>          The superdiagonal elements of the bidiagonal matrix B if
*>          m >= n, or the subdiagonal elements of B if m < n.
*> \endverbatim
*>
*> \param[in] PT
*> \verbatim
*>          PT is COMPLEX*16 array, dimension (LDPT,N)
*>          The min(m,n) by n unitary matrix P**H in the reduction
*>          A = Q * B * P**H.
*> \endverbatim
*>
*> \param[in] LDPT
*> \verbatim
*>          LDPT is INTEGER
*>          The leading dimension of the array PT.
*>          LDPT >= max(1,min(M,N)).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (M+N)
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (M)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is DOUBLE PRECISION
*>          The test ratio:
*>          norm(A - Q * B * P**H) / ( n * norm(A) * EPS )
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
*> \ingroup complex16_eig
*
*  =====================================================================
      SUBROUTINE ZBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK,
     $                   RWORK, RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            KD, LDA, LDPT, LDQ, M, N
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), RWORK( * )
      COMPLEX*16         A( LDA, * ), PT( LDPT, * ), Q( LDQ, * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ANORM, EPS
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANGE
      EXTERNAL           DLAMCH, DZASUM, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Compute A - Q * B * P**H one column at a time.
*
      RESID = ZERO
      IF( KD.NE.0 ) THEN
*
*        B is bidiagonal.
*
         IF( KD.NE.0 .AND. M.GE.N ) THEN
*
*           B is upper bidiagonal and M >= N.
*
            DO 20 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 10 I = 1, N - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   10          CONTINUE
               WORK( M+N ) = D( N )*PT( N, J )
               CALL ZGEMV( 'No transpose', M, N, -DCMPLX( ONE ), Q, LDQ,
     $                     WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   20       CONTINUE
         ELSE IF( KD.LT.0 ) THEN
*
*           B is upper bidiagonal and M < N.
*
            DO 40 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 30 I = 1, M - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   30          CONTINUE
               WORK( M+M ) = D( M )*PT( M, J )
               CALL ZGEMV( 'No transpose', M, M, -DCMPLX( ONE ), Q, LDQ,
     $                     WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   40       CONTINUE
         ELSE
*
*           B is lower bidiagonal.
*
            DO 60 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               WORK( M+1 ) = D( 1 )*PT( 1, J )
               DO 50 I = 2, M
                  WORK( M+I ) = E( I-1 )*PT( I-1, J ) +
     $                          D( I )*PT( I, J )
   50          CONTINUE
               CALL ZGEMV( 'No transpose', M, M, -DCMPLX( ONE ), Q, LDQ,
     $                     WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   60       CONTINUE
         END IF
      ELSE
*
*        B is diagonal.
*
         IF( M.GE.N ) THEN
            DO 80 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 70 I = 1, N
                  WORK( M+I ) = D( I )*PT( I, J )
   70          CONTINUE
               CALL ZGEMV( 'No transpose', M, N, -DCMPLX( ONE ), Q, LDQ,
     $                     WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   80       CONTINUE
         ELSE
            DO 100 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 90 I = 1, M
                  WORK( M+I ) = D( I )*PT( I, J )
   90          CONTINUE
               CALL ZGEMV( 'No transpose', M, M, -DCMPLX( ONE ), Q, LDQ,
     $                     WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
  100       CONTINUE
         END IF
      END IF
*
*     Compute norm(A - Q * B * P**H) / ( n * norm(A) * EPS )
*
      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      EPS = DLAMCH( 'Precision' )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         IF( ANORM.GE.RESID ) THEN
            RESID = ( RESID / ANORM ) / ( DBLE( N )*EPS )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, DBLE( N )*ANORM ) / ANORM ) /
     $                 ( DBLE( N )*EPS )
            ELSE
               RESID = MIN( RESID / ANORM, DBLE( N ) ) /
     $                 ( DBLE( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of ZBDT01
*
      END
