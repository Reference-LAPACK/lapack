      SUBROUTINE CLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDC, M, N
*     ..
*     .. Array Arguments ..
      REAL               B( LDB, * ), RWORK( * )
      COMPLEX            A( LDA, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  CLACRM performs a very simple matrix-matrix multiplication:
*           C := A * B,
*  where A is M by N and complex; B is N by N and real;
*  C is M by N and complex.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A and of the matrix C.
*          M >= 0.
*
*  N       (input) INTEGER
*          The number of columns and rows of the matrix B and
*          the number of columns of the matrix C.
*          N >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA, N)
*          A contains the M by N matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >=max(1,M).
*
*  B       (input) REAL array, dimension (LDB, N)
*          B contains the N by N matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >=max(1,N).
*
*  C       (input) COMPLEX array, dimension (LDC, N)
*          C contains the M by N matrix C.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >=max(1,N).
*
*  RWORK   (workspace) REAL array, dimension (2*M*N)
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, CMPLX, REAL
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
*
      DO 20 J = 1, N
         DO 10 I = 1, M
            RWORK( ( J-1 )*M+I ) = REAL( A( I, J ) )
   10    CONTINUE
   20 CONTINUE
*
      L = M*N + 1
      CALL SGEMM( 'N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO,
     $            RWORK( L ), M )
      DO 40 J = 1, N
         DO 30 I = 1, M
            C( I, J ) = RWORK( L+( J-1 )*M+I-1 )
   30    CONTINUE
   40 CONTINUE
*
      DO 60 J = 1, N
         DO 50 I = 1, M
            RWORK( ( J-1 )*M+I ) = AIMAG( A( I, J ) )
   50    CONTINUE
   60 CONTINUE
      CALL SGEMM( 'N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO,
     $            RWORK( L ), M )
      DO 80 J = 1, N
         DO 70 I = 1, M
            C( I, J ) = CMPLX( REAL( C( I, J ) ),
     $                  RWORK( L+( J-1 )*M+I-1 ) )
   70    CONTINUE
   80 CONTINUE
*
      RETURN
*
*     End of CLACRM
*
      END
