      SUBROUTINE CGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO )
*
*  -- LAPACK deprecated driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  This routine is deprecated and has been replaced by routine CGEQP3.
*
*  CGEQPF computes a QR factorization with column pivoting of a
*  complex M-by-N matrix A: A*P = Q*R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A. N >= 0
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the upper triangle of the array contains the
*          min(M,N)-by-N upper triangular matrix R; the elements
*          below the diagonal, together with the array TAU,
*          represent the unitary matrix Q as a product of
*          min(m,n) elementary reflectors.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
*          to the front of A*P (a leading column); if JPVT(i) = 0,
*          the i-th column of A is a free column.
*          On exit, if JPVT(i) = k, then the i-th column of A*P
*          was the k-th column of A.
*
*  TAU     (output) COMPLEX array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors.
*
*  WORK    (workspace) COMPLEX array, dimension (N)
*
*  RWORK   (workspace) REAL array, dimension (2*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(n)
*
*  Each H(i) has the form
*
*     H = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).
*
*  The matrix P is represented in jpvt as follows: If
*     jpvt(j) = i
*  then the jth column of P is the ith canonical unit vector.
*
*  Partial column norm updating strategy modified by
*    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,
*    University of Zagreb, Croatia.
*    June 2006.
*  For more details see LAPACK Working Note 176.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITEMP, J, MA, MN, PVT
      REAL               TEMP, TEMP2, TOL3Z
      COMPLEX            AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEQR2, CLARF, CLARFG, CSWAP, CUNM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, CONJG, MAX, MIN, SQRT
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SCNRM2, SLAMCH
      EXTERNAL           ISAMAX, SCNRM2, SLAMCH
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEQPF', -INFO )
         RETURN
      END IF
*
      MN = MIN( M, N )
      TOL3Z = SQRT(SLAMCH('Epsilon'))
*
*     Move initial columns up front
*
      ITEMP = 1
      DO 10 I = 1, N
         IF( JPVT( I ).NE.0 ) THEN
            IF( I.NE.ITEMP ) THEN
               CALL CSWAP( M, A( 1, I ), 1, A( 1, ITEMP ), 1 )
               JPVT( I ) = JPVT( ITEMP )
               JPVT( ITEMP ) = I
            ELSE
               JPVT( I ) = I
            END IF
            ITEMP = ITEMP + 1
         ELSE
            JPVT( I ) = I
         END IF
   10 CONTINUE
      ITEMP = ITEMP - 1
*
*     Compute the QR factorization and update remaining columns
*
      IF( ITEMP.GT.0 ) THEN
         MA = MIN( ITEMP, M )
         CALL CGEQR2( M, MA, A, LDA, TAU, WORK, INFO )
         IF( MA.LT.N ) THEN
            CALL CUNM2R( 'Left', 'Conjugate transpose', M, N-MA, MA, A,
     $                   LDA, TAU, A( 1, MA+1 ), LDA, WORK, INFO )
         END IF
      END IF
*
      IF( ITEMP.LT.MN ) THEN
*
*        Initialize partial column norms. The first n elements of
*        work store the exact column norms.
*
         DO 20 I = ITEMP + 1, N
            RWORK( I ) = SCNRM2( M-ITEMP, A( ITEMP+1, I ), 1 )
            RWORK( N+I ) = RWORK( I )
   20    CONTINUE
*
*        Compute factorization
*
         DO 40 I = ITEMP + 1, MN
*
*           Determine ith pivot column and swap if necessary
*
            PVT = ( I-1 ) + ISAMAX( N-I+1, RWORK( I ), 1 )
*
            IF( PVT.NE.I ) THEN
               CALL CSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
               ITEMP = JPVT( PVT )
               JPVT( PVT ) = JPVT( I )
               JPVT( I ) = ITEMP
               RWORK( PVT ) = RWORK( I )
               RWORK( N+PVT ) = RWORK( N+I )
            END IF
*
*           Generate elementary reflector H(i)
*
            AII = A( I, I )
            CALL CLARFG( M-I+1, AII, A( MIN( I+1, M ), I ), 1,
     $                   TAU( I ) )
            A( I, I ) = AII
*
            IF( I.LT.N ) THEN
*
*              Apply H(i) to A(i:m,i+1:n) from the left
*
               AII = A( I, I )
               A( I, I ) = CMPLX( ONE )
               CALL CLARF( 'Left', M-I+1, N-I, A( I, I ), 1,
     $                     CONJG( TAU( I ) ), A( I, I+1 ), LDA, WORK )
               A( I, I ) = AII
            END IF
*
*           Update partial column norms
*
            DO 30 J = I + 1, N
               IF( RWORK( J ).NE.ZERO ) THEN
*
*                 NOTE: The following 4 lines follow from the analysis in
*                 Lapack Working Note 176.
*                 
                  TEMP = ABS( A( I, J ) ) / RWORK( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = TEMP*( RWORK( J ) / RWORK( N+J ) )**2
                  IF( TEMP2 .LE. TOL3Z ) THEN 
                     IF( M-I.GT.0 ) THEN
                        RWORK( J ) = SCNRM2( M-I, A( I+1, J ), 1 )
                        RWORK( N+J ) = RWORK( J )
                     ELSE
                        RWORK( J ) = ZERO
                        RWORK( N+J ) = ZERO
                     END IF
                  ELSE
                     RWORK( J ) = RWORK( J )*SQRT( TEMP )
                  END IF
               END IF
   30       CONTINUE
*
   40    CONTINUE
      END IF
      RETURN
*
*     End of CGEQPF
*
      END
