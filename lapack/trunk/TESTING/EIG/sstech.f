      SUBROUTINE SSTECH( N, A, B, EIG, TOL, WORK, INFO )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
      REAL               TOL
*     ..
*     .. Array Arguments ..
      REAL               A( * ), B( * ), EIG( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*     Let T be the tridiagonal matrix with diagonal entries A(1) ,...,
*     A(N) and offdiagonal entries B(1) ,..., B(N-1)).  SSTECH checks to
*     see if EIG(1) ,..., EIG(N) are indeed accurate eigenvalues of T.
*     It does this by expanding each EIG(I) into an interval
*     [SVD(I) - EPS, SVD(I) + EPS], merging overlapping intervals if
*     any, and using Sturm sequences to count and verify whether each
*     resulting interval has the correct number of eigenvalues (using
*     SSTECT).  Here EPS = TOL*MACHEPS*MAXEIG, where MACHEPS is the
*     machine precision and MAXEIG is the absolute value of the largest
*     eigenvalue. If each interval contains the correct number of
*     eigenvalues, INFO = 0 is returned, otherwise INFO is the index of
*     the first eigenvalue in the first bad interval.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The dimension of the tridiagonal matrix T.
*
*  A       (input) REAL array, dimension (N)
*          The diagonal entries of the tridiagonal matrix T.
*
*  B       (input) REAL array, dimension (N-1)
*          The offdiagonal entries of the tridiagonal matrix T.
*
*  EIG     (input) REAL array, dimension (N)
*          The purported eigenvalues to be checked.
*
*  TOL     (input) REAL
*          Error tolerance for checking, a multiple of the
*          machine precision.
*
*  WORK    (workspace) REAL array, dimension (N)
*
*  INFO    (output) INTEGER
*          0  if the eigenvalues are all correct (to within
*             1 +- TOL*MACHEPS*MAXEIG)
*          >0 if the interval containing the INFO-th eigenvalue
*             contains the incorrect number of eigenvalues.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            BPNT, COUNT, I, ISUB, J, NUML, NUMU, TPNT
      REAL               EMIN, EPS, LOWER, MX, TUPPR, UNFLEP, UPPER
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSTECT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Check input parameters
*
      INFO = 0
      IF( N.EQ.0 )
     $   RETURN
      IF( N.LT.0 ) THEN
         INFO = -1
         RETURN
      END IF
      IF( TOL.LT.ZERO ) THEN
         INFO = -5
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      UNFLEP = SLAMCH( 'Safe minimum' ) / EPS
      EPS = TOL*EPS
*
*     Compute maximum absolute eigenvalue, error tolerance
*
      MX = ABS( EIG( 1 ) )
      DO 10 I = 2, N
         MX = MAX( MX, ABS( EIG( I ) ) )
   10 CONTINUE
      EPS = MAX( EPS*MX, UNFLEP )
*
*     Sort eigenvalues from EIG into WORK
*
      DO 20 I = 1, N
         WORK( I ) = EIG( I )
   20 CONTINUE
      DO 40 I = 1, N - 1
         ISUB = 1
         EMIN = WORK( 1 )
         DO 30 J = 2, N + 1 - I
            IF( WORK( J ).LT.EMIN ) THEN
               ISUB = J
               EMIN = WORK( J )
            END IF
   30    CONTINUE
         IF( ISUB.NE.N+1-I ) THEN
            WORK( ISUB ) = WORK( N+1-I )
            WORK( N+1-I ) = EMIN
         END IF
   40 CONTINUE
*
*     TPNT points to singular value at right endpoint of interval
*     BPNT points to singular value at left  endpoint of interval
*
      TPNT = 1
      BPNT = 1
*
*     Begin loop over all intervals
*
   50 CONTINUE
      UPPER = WORK( TPNT ) + EPS
      LOWER = WORK( BPNT ) - EPS
*
*     Begin loop merging overlapping intervals
*
   60 CONTINUE
      IF( BPNT.EQ.N )
     $   GO TO 70
      TUPPR = WORK( BPNT+1 ) + EPS
      IF( TUPPR.LT.LOWER )
     $   GO TO 70
*
*     Merge
*
      BPNT = BPNT + 1
      LOWER = WORK( BPNT ) - EPS
      GO TO 60
   70 CONTINUE
*
*     Count singular values in interval [ LOWER, UPPER ]
*
      CALL SSTECT( N, A, B, LOWER, NUML )
      CALL SSTECT( N, A, B, UPPER, NUMU )
      COUNT = NUMU - NUML
      IF( COUNT.NE.BPNT-TPNT+1 ) THEN
*
*        Wrong number of singular values in interval
*
         INFO = TPNT
         GO TO 80
      END IF
      TPNT = BPNT + 1
      BPNT = TPNT
      IF( TPNT.LE.N )
     $   GO TO 50
   80 CONTINUE
      RETURN
*
*     End of SSTECH
*
      END
