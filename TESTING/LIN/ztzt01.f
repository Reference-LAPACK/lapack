      DOUBLE PRECISION FUNCTION ZTZT01( M, N, A, AF, LDA, TAU, WORK,
     $                 LWORK )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  ZTZT01 returns
*       || A - R*Q || / ( M * eps * ||A|| )
*  for an upper trapezoidal A that was factored with ZTZRQF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrices A and AF.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and AF.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The original upper trapezoidal M by N matrix A.
*
*  AF      (input) COMPLEX*16 array, dimension (LDA,N)
*          The output of ZTZRQF for input matrix A.
*          The lower triangle is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A and AF.
*
*  TAU     (input) COMPLEX*16 array, dimension (M)
*          Details of the  Householder transformations as returned by
*          ZTZRQF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= m*n + m.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   NORMA
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZLASET, ZLATZM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX
*     ..
*     .. Executable Statements ..
*
      ZTZT01 = ZERO
*
      IF( LWORK.LT.M*N+M ) THEN
         CALL XERBLA( 'ZTZT01', 8 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, RWORK )
*
*     Copy upper triangle R
*
      CALL ZLASET( 'Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), WORK,
     $             M )
      DO 20 J = 1, M
         DO 10 I = 1, J
            WORK( ( J-1 )*M+I ) = AF( I, J )
   10    CONTINUE
   20 CONTINUE
*
*     R = R * P(1) * ... *P(m)
*
      DO 30 I = 1, M
         CALL ZLATZM( 'Right', I, N-M+1, AF( I, M+1 ), LDA, TAU( I ),
     $                WORK( ( I-1 )*M+1 ), WORK( M*M+1 ), M,
     $                WORK( M*N+1 ) )
   30 CONTINUE
*
*     R = R - A
*
      DO 40 I = 1, N
         CALL ZAXPY( M, DCMPLX( -ONE ), A( 1, I ), 1,
     $               WORK( ( I-1 )*M+1 ), 1 )
   40 CONTINUE
*
      ZTZT01 = ZLANGE( 'One-norm', M, N, WORK, M, RWORK )
*
      ZTZT01 = ZTZT01 / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N ) ) )
      IF( NORMA.NE.ZERO )
     $   ZTZT01 = ZTZT01 / NORMA
*
      RETURN
*
*     End of ZTZT01
*
      END
