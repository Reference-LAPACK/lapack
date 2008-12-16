      DOUBLE PRECISION FUNCTION DRZT01( M, N, A, AF, LDA, TAU, WORK,
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
      DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DRZT01 returns
*       || A - R*Q || / ( M * eps * ||A|| )
*  for an upper trapezoidal A that was factored with DTZRZF.
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
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The original upper trapezoidal M by N matrix A.
*
*  AF      (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The output of DTZRZF for input matrix A.
*          The lower triangle is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A and AF.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (M)
*          Details of the Householder transformations as returned by
*          DTZRZF.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= m*n + m*nb.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   NORMA
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DLASET, DORMRZ, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Executable Statements ..
*
      DRZT01 = ZERO
*
      IF( LWORK.LT.M*N+M ) THEN
         CALL XERBLA( 'DRZT01', 8 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      NORMA = DLANGE( 'One-norm', M, N, A, LDA, RWORK )
*
*     Copy upper triangle R
*
      CALL DLASET( 'Full', M, N, ZERO, ZERO, WORK, M )
      DO 20 J = 1, M
         DO 10 I = 1, J
            WORK( ( J-1 )*M+I ) = AF( I, J )
   10    CONTINUE
   20 CONTINUE
*
*     R = R * P(1) * ... *P(m)
*
      CALL DORMRZ( 'Right', 'No tranpose', M, N, M, N-M, AF, LDA, TAU,
     $             WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO )
*
*     R = R - A
*
      DO 30 I = 1, N
         CALL DAXPY( M, -ONE, A( 1, I ), 1, WORK( ( I-1 )*M+1 ), 1 )
   30 CONTINUE
*
      DRZT01 = DLANGE( 'One-norm', M, N, WORK, M, RWORK )
*
      DRZT01 = DRZT01 / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N ) ) )
      IF( NORMA.NE.ZERO )
     $   DRZT01 = DRZT01 / NORMA
*
      RETURN
*
*     End of DRZT01
*
      END
