      SUBROUTINE SLASD3( NL, NR, SQRE, K, D, Q, LDQ, DSIGMA, U, LDU, U2,
     $                   LDU2, VT, LDVT, VT2, LDVT2, IDXC, CTOT, Z,
     $                   INFO )
*
*  -- LAPACK auxiliary routine (instrumented to count ops, version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDQ, LDU, LDU2, LDVT, LDVT2, NL, NR,
     $                   SQRE
*     ..
*     .. Array Arguments ..
      INTEGER            CTOT( * ), IDXC( * )
      REAL               D( * ), DSIGMA( * ), Q( LDQ, * ), U( LDU, * ),
     $                   U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ),
     $                   Z( * )
*     ..
*     .. Common block to return operation count ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      REAL               ITCNT, OPS
*     ..
*
*  Purpose
*  =======
*
*  SLASD3 finds all the square roots of the roots of the secular
*  equation, as defined by the values in D and Z.  It makes the
*  appropriate calls to SLASD4 and then updates the singular
*  vectors by matrix multiplication.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  SLASD3 is called from SLASD1.
*
*  Arguments
*  =========
*
*  NL     (input) INTEGER
*         The row dimension of the upper block.  NL >= 1.
*
*  NR     (input) INTEGER
*         The row dimension of the lower block.  NR >= 1.
*
*  SQRE   (input) INTEGER
*         = 0: the lower block is an NR-by-NR square matrix.
*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
*
*         The bidiagonal matrix has N = NL + NR + 1 rows and
*         M = N + SQRE >= N columns.
*
*  K      (input) INTEGER
*         The size of the secular equation, 1 =< K = < N.
*
*  D      (output) REAL array, dimension(K)
*         On exit the square roots of the roots of the secular equation,
*         in ascending order.
*
*  Q      (workspace) REAL array,
*                     dimension at least (LDQ,K).
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= K.
*
*  DSIGMA (input) REAL array, dimension(K)
*         The first K elements of this array contain the old roots
*         of the deflated updating problem.  These are the poles
*         of the secular equation.
*
*  U      (input) REAL array, dimension (LDU, N)
*         The last N - K columns of this matrix contain the deflated
*         left singular vectors.
*
*  LDU    (input) INTEGER
*         The leading dimension of the array U.  LDU >= N.
*
*  U2     (input) REAL array, dimension (LDU2, N)
*         The first K columns of this matrix contain the non-deflated
*         left singular vectors for the split problem.
*
*  LDU2   (input) INTEGER
*         The leading dimension of the array U2.  LDU2 >= N.
*
*  VT     (input) REAL array, dimension (LDVT, M)
*         The last M - K columns of VT' contain the deflated
*         right singular vectors.
*
*  LDVT   (input) INTEGER
*         The leading dimension of the array VT.  LDVT >= N.
*
*  VT2    (input) REAL array, dimension (LDVT2, N)
*         The first K columns of VT2' contain the non-deflated
*         right singular vectors for the split problem.
*
*  LDVT2  (input) INTEGER
*         The leading dimension of the array VT2.  LDVT2 >= N.
*
*  IDXC   (input) INTEGER array, dimension ( N )
*         The permutation used to arrange the columns of U (and rows of
*         VT) into three groups:  the first group contains non-zero
*         entries only at and above (or before) NL +1; the second
*         contains non-zero entries only at and below (or after) NL+2;
*         and the third is dense. The first column of U and the row of
*         VT are treated separately, however.
*
*         The rows of the singular vectors found by SLASD4
*         must be likewise permuted before the matrix multiplies can
*         take place.
*
*  CTOT   (input) INTEGER array, dimension ( 4 )
*         A count of the total number of the various types of columns
*         in U (or rows in VT), as described in IDXC. The fourth column
*         type is any column which has been deflated.
*
*  Z      (input) REAL array, dimension (K)
*         The first K elements of this array contain the components
*         of the deflation-adjusted updating row vector.
*
*  INFO   (output) INTEGER
*         = 0:  successful exit.
*         < 0:  if INFO = -i, the i-th argument had an illegal value.
*         > 0:  if INFO = 1, an singular value did not converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ming Gu and Huan Ren, Computer Science Division, University of
*     California at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO, NEGONE
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0, NEGONE = -1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CTEMP, I, J, JC, KTEMP, M, N, NLP1, NLP2, NRP1
      REAL               RHO, TEMP
*     ..
*     .. External Functions ..
      REAL               SLAMC3, SNRM2, SOPBL3
      EXTERNAL           SLAMC3, SNRM2, SOPBL3
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMM, SLACPY, SLASCL, SLASD4, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( NL.LT.1 ) THEN
         INFO = -1
      ELSE IF( NR.LT.1 ) THEN
         INFO = -2
      ELSE IF( ( SQRE.NE.1 ) .AND. ( SQRE.NE.0 ) ) THEN
         INFO = -3
      END IF
*
      N = NL + NR + 1
      M = N + SQRE
      NLP1 = NL + 1
      NLP2 = NL + 2
*
      IF( ( K.LT.1 ) .OR. ( K.GT.N ) ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.K ) THEN
         INFO = -7
      ELSE IF( LDU.LT.N ) THEN
         INFO = -10
      ELSE IF( LDU2.LT.N ) THEN
         INFO = -12
      ELSE IF( LDVT.LT.M ) THEN
         INFO = -14
      ELSE IF( LDVT2.LT.M ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASD3', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( K.EQ.1 ) THEN
         D( 1 ) = ABS( Z( 1 ) )
         CALL SCOPY( M, VT2( 1, 1 ), LDVT2, VT( 1, 1 ), LDVT )
         IF( Z( 1 ).GT.ZERO ) THEN
            CALL SCOPY( N, U2( 1, 1 ), 1, U( 1, 1 ), 1 )
         ELSE
            DO 10 I = 1, N
               U( I, 1 ) = -U2( I, 1 )
   10       CONTINUE
         END IF
         RETURN
      END IF
*
*     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
*     be computed with high relative accuracy (barring over/underflow).
*     This is a problem on machines without a guard digit in
*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
*     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
*     which on any of these machines zeros out the bottommost
*     bit of DSIGMA(I) if it is 1; this makes the subsequent
*     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
*     occurs. On binary machines with a guard digit (almost all
*     machines) it does not change DSIGMA(I) at all. On hexadecimal
*     and decimal machines with a guard digit, it slightly
*     changes the bottommost bits of DSIGMA(I). It does not account
*     for hexadecimal or decimal machines without guard digits
*     (we know of none). We use a subroutine call to compute
*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
*     this code.
*
      DO 20 I = 1, K
         DSIGMA( I ) = SLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
   20 CONTINUE
*
*     Keep a copy of Z.
*
      CALL SCOPY( K, Z, 1, Q, 1 )
*
*     Normalize Z.
*
      OPS = OPS + REAL( K*3 + 1)
      RHO = SNRM2( K, Z, 1 )
      CALL SLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
*
*     Find the new singular values.
*
      DO 30 J = 1, K
         CALL SLASD4( K, J, DSIGMA, Z, U( 1, J ), RHO, D( J ),
     $                VT( 1, J ), INFO )
*
*        If the zero finder fails, the computation is terminated.
*
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
   30 CONTINUE
*
*     Compute updated Z.
*
      OPS = OPS + REAL( K*2 )
      DO 60 I = 1, K
         Z( I ) = U( I, K )*VT( I, K )
         OPS = OPS + REAL( (I-1)*6 )
         DO 40 J = 1, I - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) /
     $               ( DSIGMA( I )-DSIGMA( J ) ) /
     $               ( DSIGMA( I )+DSIGMA( J ) ) )
   40    CONTINUE
         OPS = OPS + REAL( (K-I)*6 )
         DO 50 J = I, K - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) /
     $               ( DSIGMA( I )-DSIGMA( J+1 ) ) /
     $               ( DSIGMA( I )+DSIGMA( J+1 ) ) )
   50    CONTINUE
         Z( I ) = SIGN( SQRT( ABS( Z( I ) ) ), Q( I, 1 ) )
   60 CONTINUE
*
*     Compute left singular vectors of the modified diagonal matrix,
*     and store related information for the right singular vectors.
*
      OPS = OPS + REAL( K*(3+K*2) + MAX(0,(K-1)*4) )
      DO 90 I = 1, K
         VT( 1, I ) = Z( 1 ) / U( 1, I ) / VT( 1, I )
         U( 1, I ) = NEGONE
         DO 70 J = 2, K
            VT( J, I ) = Z( J ) / U( J, I ) / VT( J, I )
            U( J, I ) = DSIGMA( J )*VT( J, I )
   70    CONTINUE
         TEMP = SNRM2( K, U( 1, I ), 1 )
         Q( 1, I ) = U( 1, I ) / TEMP
         DO 80 J = 2, K
            JC = IDXC( J )
            Q( J, I ) = U( JC, I ) / TEMP
   80    CONTINUE
   90 CONTINUE
*
*     Update the left singular vector matrix.
*
      IF( K.EQ.2 ) THEN
         OPS = OPS + SOPBL3( 'SGEMM ', N, K, K ) 
         CALL SGEMM( 'N', 'N', N, K, K, ONE, U2, LDU2, Q, LDQ, ZERO, U,
     $               LDU )
         GO TO 100
      END IF
      IF( CTOT( 1 ).GT.0 ) THEN
         OPS = OPS + SOPBL3( 'SGEMM ', NL, K, CTOT( 1 ) )
         CALL SGEMM( 'N', 'N', NL, K, CTOT( 1 ), ONE, U2( 1, 2 ), LDU2,
     $               Q( 2, 1 ), LDQ, ZERO, U( 1, 1 ), LDU )
         IF( CTOT( 3 ).GT.0 ) THEN
            KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
            OPS = OPS + SOPBL3( 'SGEMM ', NL, K, CTOT( 3 ) )
            CALL SGEMM( 'N', 'N', NL, K, CTOT( 3 ), ONE, U2( 1, KTEMP ),
     $                  LDU2, Q( KTEMP, 1 ), LDQ, ONE, U( 1, 1 ), LDU )
         END IF
      ELSE IF( CTOT( 3 ).GT.0 ) THEN
         KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
         OPS = OPS + SOPBL3( 'SGEMM ', NL, K, CTOT( 3 ) )
         CALL SGEMM( 'N', 'N', NL, K, CTOT( 3 ), ONE, U2( 1, KTEMP ),
     $               LDU2, Q( KTEMP, 1 ), LDQ, ZERO, U( 1, 1 ), LDU )
      ELSE
         CALL SLACPY( 'F', NL, K, U2, LDU2, U, LDU )
      END IF
      CALL SCOPY( K, Q( 1, 1 ), LDQ, U( NLP1, 1 ), LDU )
      KTEMP = 2 + CTOT( 1 )
      CTEMP = CTOT( 2 ) + CTOT( 3 )
      OPS = OPS + SOPBL3( 'SGEMM ', NR, K, CTEMP )
      CALL SGEMM( 'N', 'N', NR, K, CTEMP, ONE, U2( NLP2, KTEMP ), LDU2,
     $            Q( KTEMP, 1 ), LDQ, ZERO, U( NLP2, 1 ), LDU )
*
*     Generate the right singular vectors.
*
  100 CONTINUE
      OPS = OPS + REAL( K*(K*2+1) + MAX(0,K-1) )
      DO 120 I = 1, K
         TEMP = SNRM2( K, VT( 1, I ), 1 )
         Q( I, 1 ) = VT( 1, I ) / TEMP
         DO 110 J = 2, K
            JC = IDXC( J )
            Q( I, J ) = VT( JC, I ) / TEMP
  110    CONTINUE
  120 CONTINUE
*
*     Update the right singular vector matrix.
*
      IF( K.EQ.2 ) THEN
         OPS = OPS + SOPBL3( 'SGEMM ', K, M, K ) 
         CALL SGEMM( 'N', 'N', K, M, K, ONE, Q, LDQ, VT2, LDVT2, ZERO,
     $               VT, LDVT )
         RETURN
      END IF
      KTEMP = 1 + CTOT( 1 )
      OPS = OPS + SOPBL3( 'SGEMM ', K, NLP1, KTEMP )
      CALL SGEMM( 'N', 'N', K, NLP1, KTEMP, ONE, Q( 1, 1 ), LDQ,
     $            VT2( 1, 1 ), LDVT2, ZERO, VT( 1, 1 ), LDVT )
      KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
      OPS = OPS + SOPBL3( 'SGEMM ', K, NLP1, CTOT( 3 ) )
      IF( KTEMP.LE.LDVT2 )
     $   CALL SGEMM( 'N', 'N', K, NLP1, CTOT( 3 ), ONE, Q( 1, KTEMP ),
     $               LDQ, VT2( KTEMP, 1 ), LDVT2, ONE, VT( 1, 1 ),
     $               LDVT )
*
      KTEMP = CTOT( 1 ) + 1
      NRP1 = NR + SQRE
      IF( KTEMP.GT.1 ) THEN
         DO 130 I = 1, K
            Q( I, KTEMP ) = Q( I, 1 )
  130    CONTINUE
         DO 140 I = NLP2, M
            VT2( KTEMP, I ) = VT2( 1, I )
  140    CONTINUE
      END IF
      CTEMP = 1 + CTOT( 2 ) + CTOT( 3 )
      OPS = OPS + SOPBL3( 'SGEMM ', K, NRP1, CTEMP ) 
      CALL SGEMM( 'N', 'N', K, NRP1, CTEMP, ONE, Q( 1, KTEMP ), LDQ,
     $            VT2( KTEMP, NLP2 ), LDVT2, ZERO, VT( 1, NLP2 ), LDVT )
*
      RETURN
*
*     End of SLASD3
*
      END
