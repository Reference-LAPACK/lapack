      SUBROUTINE DTIMMG( IFLAG, M, N, A, LDA, KL, KU )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      INTEGER            IFLAG, KL, KU, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTIMMG generates a real test matrix whose type is given by IFLAG.
*  All the matrices are Toeplitz (constant along a diagonal), with
*  random elements on each diagonal.
*
*  Arguments
*  =========
*
*  IFLAG   (input) INTEGER
*          The type of matrix to be generated.
*          = 0 or 1:   General matrix
*          = 2 or -2:  General banded matrix
*          = 3 or -3:  Symmetric positive definite matrix
*          = 4 or -4:  Symmetric positive definite packed
*          = 5 or -5:  Symmetric positive definite banded
*          = 6 or -6:  Symmetric indefinite matrix
*          = 7 or -7:  Symmetric indefinite packed
*          = 8 or -8:  Symmetric indefinite banded
*          = 9 or -9:  Triangular
*          = 10 or -10:  Triangular packed
*          = 11 or -11:  Triangular banded
*          = 12:         General tridiagonal
*          = 13 or -13:  Positive definite tridiagonal
*          For symmetric or triangular matrices, IFLAG > 0 indicates
*          upper triangular storage and IFLAG < 0 indicates lower
*          triangular storage.
*
*  M       (input) INTEGER
*          The number of rows of the matrix to be generated.
*
*  N       (input) INTEGER
*          The number of columns of the matrix to be generated.
*
*  A       (output) DOUBLE PRECISION array, dimension (LDA,N)
*          The generated matrix.
*
*          If the absolute value of IFLAG is 1, 3, or 6, the leading
*          M x N (or N x N) subblock is used to store the matrix.
*          If the matrix is symmetric, only the upper or lower triangle
*          of this block is referenced.
*
*          If the absolute value of IFLAG is 4 or 7, the matrix is
*          symmetric and packed storage is used for the upper or lower
*          triangle.  The triangular matrix is stored columnwise as a
*          inear array, and the array A is treated as a vector of
*          length LDA.  LDA must be set to at least N*(N+1)/2.
*
*          If the absolute value of IFLAG is 2 or 5, the matrix is
*          returned in band format.  The columns of the matrix are
*          specified in the columns of A and the diagonals of the
*          matrix are specified in the rows of A, with the leading
*          diagonal in row
*              KL + KU + 1,  if IFLAG = 2
*              KU + 1,       if IFLAG = 5 or -2
*              1,            if IFLAG = -5
*          If IFLAG = 2, the first KL rows are not used to leave room
*          for pivoting in DGBTRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  If the generated matrix is
*          packed, LDA >= N*(N+1)/2, otherwise LDA >= max(1,M).
*
*  KL      (input) INTEGER
*          The number of subdiagonals if IFLAG = 2, 5, or -5.
*
*  KU      (input) INTEGER
*          The number of superdiagonals if IFLAG = 2, 5, or -5.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, JJ, JN, K, MJ, MU
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLARNV
*     ..
*     .. Data statements ..
      DATA               ISEED / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RETURN
*
      ELSE IF( IFLAG.EQ.0 .OR. IFLAG.EQ.1 ) THEN
*
*        General matrix
*
*        Set first column and row to random values.
*
         CALL DLARNV( 2, ISEED, M, A( 1, 1 ) )
         DO 10 J = 2, N, M
            MJ = MIN( M, N-J+1 )
            CALL DLARNV( 2, ISEED, MJ, A( 1, J ) )
            IF( MJ.GT.1 )
     $         CALL DCOPY( MJ-1, A( 2, J ), 1, A( 1, J+1 ), LDA )
   10    CONTINUE
*
*        Fill in the rest of the matrix.
*
         DO 30 J = 2, N
            DO 20 I = 2, M
               A( I, J ) = A( I-1, J-1 )
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( IFLAG.EQ.2 .OR. IFLAG.EQ.-2 ) THEN
*
*        General band matrix
*
         IF( IFLAG.EQ.2 ) THEN
            K = KL + KU + 1
         ELSE
            K = KU + 1
         END IF
         CALL DLARNV( 2, ISEED, MIN( M, KL+1 ), A( K, 1 ) )
         MU = MIN( N-1, KU )
         CALL DLARNV( 2, ISEED, MU+1, A( K-MU, N ) )
         DO 40 J = 2, N - 1
            MU = MIN( J-1, KU )
            CALL DCOPY( MU, A( K-MU, N ), 1, A( K-MU, J ), 1 )
            CALL DCOPY( MIN( M-J+1, KL+1 ), A( K, 1 ), 1, A( K, J ), 1 )
   40    CONTINUE
*
      ELSE IF( IFLAG.EQ.3 ) THEN
*
*        Symmetric positive definite, upper triangle
*
         CALL DLARNV( 2, ISEED, N-1, A( 1, N ) )
         A( N, N ) = DBLE( N )
         DO 50 J = N - 1, 1, -1
            CALL DCOPY( J, A( N-J+1, N ), 1, A( 1, J ), 1 )
   50    CONTINUE
*
      ELSE IF( IFLAG.EQ.-3 ) THEN
*
*        Symmetric positive definite, lower triangle
*
         A( 1, 1 ) = DBLE( N )
         IF( N.GT.1 )
     $      CALL DLARNV( 2, ISEED, N-1, A( 2, 1 ) )
         DO 60 J = 2, N
            CALL DCOPY( N-J+1, A( 1, 1 ), 1, A( J, J ), 1 )
   60    CONTINUE
*
      ELSE IF( IFLAG.EQ.4 ) THEN
*
*        Symmetric positive definite packed, upper triangle
*
         JN = ( N-1 )*N / 2 + 1
         CALL DLARNV( 2, ISEED, N-1, A( JN, 1 ) )
         A( JN+N-1, 1 ) = DBLE( N )
         JJ = JN
         DO 70 J = N - 1, 1, -1
            JJ = JJ - J
            JN = JN + 1
            CALL DCOPY( J, A( JN, 1 ), 1, A( JJ, 1 ), 1 )
   70    CONTINUE
*
      ELSE IF( IFLAG.EQ.-4 ) THEN
*
*        Symmetric positive definite packed, lower triangle
*
         A( 1, 1 ) = DBLE( N )
         IF( N.GT.1 )
     $      CALL DLARNV( 2, ISEED, N-1, A( 2, 1 ) )
         JJ = N + 1
         DO 80 J = 2, N
            CALL DCOPY( N-J+1, A( 1, 1 ), 1, A( JJ, 1 ), 1 )
            JJ = JJ + N - J + 1
   80    CONTINUE
*
      ELSE IF( IFLAG.EQ.5 ) THEN
*
*        Symmetric positive definite banded, upper triangle
*
         K = KL
         MU = MIN( N-1, K )
         CALL DLARNV( 2, ISEED, MU, A( K+1-MU, N ) )
         A( K+1, N ) = DBLE( N )
         DO 90 J = N - 1, 1, -1
            MU = MIN( J, K+1 )
            CALL DCOPY( MU, A( K+2-MU, N ), 1, A( K+2-MU, J ), 1 )
   90    CONTINUE
*
      ELSE IF( IFLAG.EQ.-5 ) THEN
*
*        Symmetric positive definite banded, lower triangle
*
         K = KL
         A( 1, 1 ) = DBLE( N )
         CALL DLARNV( 2, ISEED, MIN( N-1, K ), A( 2, 1 ) )
         DO 100 J = 2, N
            CALL DCOPY( MIN( N-J+1, K+1 ), A( 1, 1 ), 1, A( 1, J ), 1 )
  100    CONTINUE
*
      ELSE IF( IFLAG.EQ.6 ) THEN
*
*        Symmetric indefinite, upper triangle
*
         CALL DLARNV( 2, ISEED, N, A( 1, N ) )
         DO 110 J = N - 1, 1, -1
            CALL DCOPY( J, A( N-J+1, N ), 1, A( 1, J ), 1 )
  110    CONTINUE
*
      ELSE IF( IFLAG.EQ.-6 ) THEN
*
*        Symmetric indefinite, lower triangle
*
         CALL DLARNV( 2, ISEED, N, A( 1, 1 ) )
         DO 120 J = 2, N
            CALL DCOPY( N-J+1, A( 1, 1 ), 1, A( J, J ), 1 )
  120    CONTINUE
*
      ELSE IF( IFLAG.EQ.7 ) THEN
*
*        Symmetric indefinite packed, upper triangle
*
         JN = ( N-1 )*N / 2 + 1
         CALL DLARNV( 2, ISEED, N, A( JN, 1 ) )
         JJ = JN
         DO 130 J = N - 1, 1, -1
            JJ = JJ - J
            JN = JN + 1
            CALL DCOPY( J, A( JN, 1 ), 1, A( JJ, 1 ), 1 )
  130    CONTINUE
*
      ELSE IF( IFLAG.EQ.-7 ) THEN
*
*        Symmetric indefinite packed, lower triangle
*
         CALL DLARNV( 2, ISEED, N, A( 1, 1 ) )
         JJ = N + 1
         DO 140 J = 2, N
            CALL DCOPY( N-J+1, A( 1, 1 ), 1, A( JJ, 1 ), 1 )
            JJ = JJ + N - J + 1
  140    CONTINUE
*
      ELSE IF( IFLAG.EQ.8 ) THEN
*
*        Symmetric indefinite banded, upper triangle
*
         K = KL
         MU = MIN( N, K+1 )
         CALL DLARNV( 2, ISEED, MU, A( K+2-MU, N ) )
         DO 150 J = N - 1, 1, -1
            MU = MIN( J, K+1 )
            CALL DCOPY( MU, A( K+2-MU, N ), 1, A( K+2-MU, J ), 1 )
  150    CONTINUE
*
      ELSE IF( IFLAG.EQ.-8 ) THEN
*
*        Symmetric indefinite banded, lower triangle
*
         K = KL
         CALL DLARNV( 2, ISEED, MIN( N, K+1 ), A( 1, 1 ) )
         DO 160 J = 2, N
            CALL DCOPY( MIN( N-J+1, K+1 ), A( 1, 1 ), 1, A( 1, J ), 1 )
  160    CONTINUE
*
      ELSE IF( IFLAG.EQ.9 ) THEN
*
*        Upper triangular
*
         CALL DLARNV( 2, ISEED, N, A( 1, N ) )
         A( N, N ) = SIGN( DBLE( N ), A( N, N ) )
         DO 170 J = N - 1, 1, -1
            CALL DCOPY( J, A( N-J+1, N ), 1, A( 1, J ), 1 )
  170    CONTINUE
*
      ELSE IF( IFLAG.EQ.-9 ) THEN
*
*        Lower triangular
*
         CALL DLARNV( 2, ISEED, N, A( 1, 1 ) )
         A( 1, 1 ) = SIGN( DBLE( N ), A( 1, 1 ) )
         DO 180 J = 2, N
            CALL DCOPY( N-J+1, A( 1, 1 ), 1, A( J, J ), 1 )
  180    CONTINUE
*
      ELSE IF( IFLAG.EQ.10 ) THEN
*
*        Upper triangular packed
*
         JN = ( N-1 )*N / 2 + 1
         CALL DLARNV( 2, ISEED, N, A( JN, 1 ) )
         A( JN+N-1, 1 ) = SIGN( DBLE( N ), A( JN+N-1, 1 ) )
         JJ = JN
         DO 190 J = N - 1, 1, -1
            JJ = JJ - J
            JN = JN + 1
            CALL DCOPY( J, A( JN, 1 ), 1, A( JJ, 1 ), 1 )
  190    CONTINUE
*
      ELSE IF( IFLAG.EQ.-10 ) THEN
*
*        Lower triangular packed
*
         CALL DLARNV( 2, ISEED, N, A( 1, 1 ) )
         A( 1, 1 ) = SIGN( DBLE( N ), A( 1, 1 ) )
         JJ = N + 1
         DO 200 J = 2, N
            CALL DCOPY( N-J+1, A( 1, 1 ), 1, A( JJ, 1 ), 1 )
            JJ = JJ + N - J + 1
  200    CONTINUE
*
      ELSE IF( IFLAG.EQ.11 ) THEN
*
*        Upper triangular banded
*
         K = KL
         MU = MIN( N, K+1 )
         CALL DLARNV( 2, ISEED, MU, A( K+2-MU, N ) )
         A( K+1, N ) = SIGN( DBLE( K+1 ), A( K+1, N ) )
         DO 210 J = N - 1, 1, -1
            MU = MIN( J, K+1 )
            CALL DCOPY( MU, A( K+2-MU, N ), 1, A( K+2-MU, J ), 1 )
  210    CONTINUE
*
      ELSE IF( IFLAG.EQ.-11 ) THEN
*
*        Lower triangular banded
*
         K = KL
         CALL DLARNV( 2, ISEED, MIN( N, K+1 ), A( 1, 1 ) )
         A( 1, 1 ) = SIGN( DBLE( K+1 ), A( 1, 1 ) )
         DO 220 J = 2, N
            CALL DCOPY( MIN( N-J+1, K+1 ), A( 1, 1 ), 1, A( 1, J ), 1 )
  220    CONTINUE
*
      ELSE IF( IFLAG.EQ.12 ) THEN
*
*        General tridiagonal
*
         CALL DLARNV( 2, ISEED, 3*N-2, A )
*
      ELSE IF( IFLAG.EQ.13 .OR. IFLAG.EQ.-13 ) THEN
*
*        Positive definite tridiagonal
*
         DO 230 J = 1, N
            A( J, 1 ) = 2.0D0
  230    CONTINUE
         CALL DLARNV( 2, ISEED, N-1, A( N+1, 1 ) )
      END IF
*
      RETURN
*
*     End of DTIMMG
*
      END
