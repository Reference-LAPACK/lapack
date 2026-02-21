*> \brief \b DKTEBZ
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DKTEBZ + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dktebz.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dktebz.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dktebz.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DKTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, E,
*                          M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          ORDER, RANGE
*       INTEGER            IL, INFO, IU, M, N, NSPLIT
*       DOUBLE PRECISION   ABSTOL, VL, VU
*       ..
*       .. Array Arguments ..
*       INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
*       DOUBLE PRECISION   E( * ), W( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DKTEBZ computes the eigenvalues of a skew-symmetric tridiagonal
*> matrix T.  The user may ask for all eigenvalues, all eigenvalues
*> in the half-open interval (VL,VU], or the IL-th through IU-th
*> eigenvalues.
*>
*> Since the eigenvalues of skew-symmetric matrix are conjugate purely
*> imaginary pairs, to be consistent with output layout of subroutine
*> dkyev, the output layout according to VL/VU and IL/IU abides the
*> following rules (suppose ORDER = 'E' here, and the same for ORDER = 'B'):
*>
*> If RANGE = 'A', the same as layout of dkyev.
*>
*> If RANGE = 'V', based on the output of RANGE = 'A', but excluding the
*> values not in the interval (if the value is positive, also excluding the
*> 0 followed by it).
*>
*> If RANGE = 'I' and N is odd, sort all the eigenvalues in ascending order
*> of their imaginary parts. The index starts from the central element,
*> which occupy the index 1 exclusively. Then the index increases in the
*> positive direction, and each index corresponds to two eigenvalues (one
*> is at the position of index, and the other is at the symmetric position
*> of former). These two eigenvalues are either two zeros, or a pair of conjugate
*> purely imaginary values, and are organized in the same layout as dkyev.
*>
*> If RANGE = 'I' and N is even, the same as N being odd, except that the index
*> starts from the first element beside the center in the positive direction,
*> and all the indexes correspond to two eigenvalues.
*>
*> This subroutine reuses the bisection mathod of symmetric matrix, since the
*> symmetrix and skew-symmetrix tridiagonal with the same E have consistent
*> eigenvalues except the scale of i. Range and indexes are adjusted to produce
*> only the position values. The organization of zeros is derived from the
*> parity of block sizes.
*>
*> To avoid overflow, the matrix must be scaled so that its
*> largest element is no greater than overflow**(1/2) * underflow**(1/4)
*> in absolute value, and for greatest accuracy, it should not 
*> be much smaller than that.
*>
*> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
*> Matrix", Report CS41, Computer Science Dept., Stanford
*> University, July 21, 1966.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] RANGE
*> \verbatim
*>          RANGE is CHARACTER*1
*>          = 'A': ("All")   all eigenvalues will be found.
*>          = 'V': ("Value") all eigenvalues in the half-open interval
*>                           (VL,VU] will be found.
*>          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
*>                           entire matrix) will be found.
*> \endverbatim
*>
*> \param[in] ORDER
*> \verbatim
*>          ORDER is CHARACTER*1
*>          = 'B': ("By Block") the eigenvalues will be grouped by
*>                              split-off block (see IBLOCK, ISPLIT) and
*>                              ordered from smallest to largest within
*>                              the block.
*>          = 'E': ("Entire matrix")
*>                              the eigenvalues for the entire matrix
*>                              will be ordered from smallest to
*>                              largest.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the tridiagonal matrix T.  N >= 0.
*> \endverbatim
*>
*> \param[in] VL
*> \verbatim
*>          VL is DOUBLE PRECISION
*>
*>          If RANGE='V', the lower bound of the interval to
*>          be searched for eigenvalues.  Eigenvalues with positive
*>          imaginary part less than or equal to VL, or greater
*>          than VU, will not be returned.  VL < VU. If VL < 0,
*>          zero eigenvalues will be included.
*>          Not referenced if RANGE = 'A' or 'I'.
*> \endverbatim
*>
*> \param[in] VU
*> \verbatim
*>          VU is DOUBLE PRECISION
*>
*>          If RANGE='V', the upper bound of the interval to
*>          be searched for eigenvalues.  Eigenvalues with positive
*>          imaginary part less than or equal to VL, or greater
*>          than VU, will not be returned.  VL < VU, VU >= 0.
*>          Not referenced if RANGE = 'A' or 'I'.
*> \endverbatim
*>
*> \param[in] IL
*> \verbatim
*>          IL is INTEGER
*>
*>          If RANGE='I', the index of eigenvalue with the
*>          smallest positive imaginary part to be returned.
*>          1 <= IL <= IU <= ceil(N/2), if N > 0; IL = 1 and
*>          IU = 0 if N = 0.
*>          Not referenced if RANGE = 'A' or 'V'.
*> \endverbatim
*>
*> \param[in] IU
*> \verbatim
*>          IU is INTEGER
*>
*>          If RANGE='I', the index of eigenvalue with the
*>          largest positive imaginary part to be returned.
*>          1 <= IL <= IU <= ceil(N/2), if N > 0; IL = 1 and
*》         IU = 0 if N = 0.
*>          Not referenced if RANGE = 'A' or 'V'.
*> \endverbatim
*>
*> \param[in] ABSTOL
*> \verbatim
*>          ABSTOL is DOUBLE PRECISION
*>          The absolute tolerance for the eigenvalues.  An eigenvalue
*>          (or cluster) is considered to be located if it has been
*>          determined to lie in an interval whose width is ABSTOL or
*>          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
*>          will be used, where |T| means the 1-norm of T.
*>
*>          Eigenvalues will be computed most accurately when ABSTOL is
*>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          The (n-1) off-diagonal elements of the tridiagonal matrix T.
*> \endverbatim
*>
*> \param[out] M
*> \verbatim
*>          M is INTEGER
*>          The actual number of eigenvalues found. 0 <= M <= N.
*>          (See also the description of INFO=2,3.)
*> \endverbatim
*>
*> \param[out] NSPLIT
*> \verbatim
*>          NSPLIT is INTEGER
*>          The number of diagonal blocks in the matrix T.
*>          1 <= NSPLIT <= N.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is DOUBLE PRECISION array, dimension (N)
*>          On exit, the first M elements of W will contain the
*>          eigenvalues.  (DKTEBZ may use the remaining N-M elements as
*>          workspace.)
*> \endverbatim
*>
*> \param[out] IBLOCK
*> \verbatim
*>          IBLOCK is INTEGER array, dimension (N)
*>          At each row/column j where E(j) is zero or small, the
*>          matrix T is considered to split into a block diagonal
*>          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
*>          block (from 1 to the number of blocks) the eigenvalue W(i)
*>          belongs.  (DKTEBZ may use the remaining N-M elements as
*>          workspace.)
*> \endverbatim
*>
*> \param[out] ISPLIT
*> \verbatim
*>          ISPLIT is INTEGER array, dimension (N)
*>          The splitting points, at which T breaks up into submatrices.
*>          The first submatrix consists of rows/columns 1 to ISPLIT(1),
*>          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
*>          etc., and the NSPLIT-th consists of rows/columns
*>          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
*>          (Only the first NSPLIT elements will actually be used, but
*>          since the user cannot know a priori what value NSPLIT will
*>          have, N words must be reserved for ISPLIT.)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (5*N)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (3*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  some or all of the eigenvalues failed to converge or
*>                were not computed:
*>                =1 or 3: Bisection failed to converge for some
*>                        eigenvalues; these eigenvalues are flagged by a
*>                        negative block number.  The effect is that the
*>                        eigenvalues may not be as accurate as the
*>                        absolute and relative tolerances.  This is
*>                        generally caused by unexpectedly inaccurate
*>                        arithmetic.
*>                =2 or 3: RANGE='I' only: Not all of the eigenvalues
*>                        IL:IU were found.
*>                        Effect: M < IU+1-IL
*>                        Cause:  non-monotonic arithmetic, causing the
*>                                Sturm sequence to be non-monotonic.
*>                        Cure:   recalculate, using RANGE='A', and pick
*>                                out eigenvalues IL:IU.  In some cases,
*>                                increasing the PARAMETER "FUDGE" may
*>                                make things work.
*>                = 4:    RANGE='I', and the Gershgorin interval
*>                        initially used was too small.  No eigenvalues
*>                        were computed.
*>                        Probable cause: your machine has sloppy
*>                                        floating-point arithmetic.
*>                        Cure: Increase the PARAMETER "FUDGE",
*>                              recompile, and try again.
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  RELFAC  DOUBLE PRECISION, default = 2.0e0
*>          The relative tolerance.  An interval (a,b] lies within
*>          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
*>          where "ulp" is the machine precision (distance from 1 to
*>          the next larger floating point number.)
*>
*>  FUDGE   DOUBLE PRECISION, default = 2
*>          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
*>          a value of 1 should work, but on machines with sloppy
*>          arithmetic, this needs to be larger.  The default for
*>          publicly released versions should be large enough to handle
*>          the worst machine around.  Note that this has no effect
*>          on accuracy of the solution.
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
*> \ingroup stebz
*
*  =====================================================================
      SUBROUTINE DKTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, E,
     $                   M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,
     $                   INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          ORDER, RANGE
      INTEGER            IL, INFO, IU, M, N, NSPLIT
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
      DOUBLE PRECISION   E( * ), W( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   HALF = 1.0D0 / TWO )
      DOUBLE PRECISION   FUDGE, RELFAC
      PARAMETER          ( FUDGE = 2.1D0, RELFAC = 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NCNVRG, TOOFEW
      INTEGER            IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO,
     $                   IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX,
     $                   ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWL,
     $                   NWU, ILCVT, IUCVT, ODDCT, ITMP2
      DOUBLE PRECISION   ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN,
     $                   TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL
*     ..
*     .. Local Arrays ..
      INTEGER            IDUMMA( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, ILAENV, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAEBZ, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Decode RANGE
*
      IF( LSAME( RANGE, 'A' ) ) THEN
         IRANGE = 1
      ELSE IF( LSAME( RANGE, 'V' ) ) THEN
         IRANGE = 2
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         IRANGE = 3
      ELSE
         IRANGE = 0
      END IF
*
*     Decode ORDER
*
      IF( LSAME( ORDER, 'B' ) ) THEN
         IORDER = 2
      ELSE IF( LSAME( ORDER, 'E' ) ) THEN
         IORDER = 1
      ELSE
         IORDER = 0
      END IF
*
*     Check for Errors
*
      IF( IRANGE.LE.0 ) THEN
         INFO = -1
      ELSE IF( IORDER.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( IRANGE.EQ.2 ) THEN
         IF( VL.GE.VU .OR. VU.LT.ZERO ) INFO = -5
      ELSE IF( IRANGE.EQ.3 .AND. ( IL.LT.1 .OR.
     $         IL.GT.MAX( 1, (N+1)/2 ) ) )  THEN
         INFO = -6
      ELSE IF( IRANGE.EQ.3 .AND. ( IU.LT.MIN( N, IL ).OR.
     $         IU.GT.(N+1)/2 ) )  THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DKTEBZ', -INFO )
         RETURN
      END IF
*
*     Initialize error flags
*
      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 )
     $   RETURN
*
*     Simplifications:
*
      IF( IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.(N+1)/2 )
     $   IRANGE = 1
*
*     Get machine constants
*     NB is the minimum vector length for vector bisection, or 0
*     if only scalar is to be done.
*
      SAFEMN = DLAMCH( 'S' )
      ULP = DLAMCH( 'P' )
      RTOLI = ULP*RELFAC
      NB = ILAENV( 1, 'DKTEBZ', ' ', N, -1, -1, -1 )
      IF( NB.LE.1 )
     $   NB = 0
*
*     Special Case when N=1
*
      IF( N.EQ.1 ) THEN
         NSPLIT = 1
         ISPLIT( 1 ) = 1
         IF( IRANGE.EQ.2 .AND. VL.GE.ZERO ) THEN
            M = 0
         ELSE
            W( 1 ) = ZERO
            IBLOCK( 1 ) = 1
            M = 1
         END IF
         RETURN
      ELSEIF( N.EQ.2 ) THEN
         NSPLIT = 1
         ISPLIT( 1 ) = 2
         IF( IRANGE.EQ.2 .AND. ( VL.GE.ABS( E( 1 ) ) .OR.
     $   VU.LT.ABS( E( 1 ) ) ) ) THEN
            M = 0
         ELSE
            W( 1 ) = ABS( E( 1 ) )
            W( 2 ) = ZERO
            IBLOCK( 1 ) = 1
            IBLOCK( 2 ) = 1
            M = 2
         END IF
         RETURN
      END IF
*
      DO J = 1, N
         WORK( 4*N + J ) = ZERO
      END DO
*
*     Compute Splitting Points
*
      NSPLIT = 1
      WORK( N ) = ZERO
      PIVMIN = ONE
      ODDCT = 0
*
      DO 10 J = 2, N
         TMP1 = E( J-1 )**2
         IF ( J.EQ.2 ) THEN
            TMP2 = E( J )**2
         ELSEIF ( J.EQ.N ) THEN
            TMP2 = E( J-2 )**2
         ELSE
            TMP2 = E( J-2 )*E( J )
         END IF
         IF( ABS( TMP2 )*ULP**2+SAFEMN.GT.TMP1 ) THEN
            ISPLIT( NSPLIT ) = J - 1
            IF( NSPLIT.EQ.1 ) THEN
               ODDCT = ODDCT + MOD( ISPLIT( NSPLIT ), 2 )
            ELSE
               ODDCT = ODDCT + MOD( ISPLIT( NSPLIT ) -
     $                 ISPLIT( NSPLIT-1 ), 2 )
            END IF
            NSPLIT = NSPLIT + 1
            WORK( J-1 ) = ZERO
         ELSE
            WORK( J-1 ) = TMP1
            PIVMIN = MAX( PIVMIN, TMP1 )
         END IF
   10 CONTINUE
      ISPLIT( NSPLIT ) = N
      IF( NSPLIT.EQ.1 ) THEN
         ODDCT = ODDCT + MOD( ISPLIT( NSPLIT ), 2 )
      ELSE
         ODDCT = ODDCT + MOD( ISPLIT( NSPLIT ) -
     $           ISPLIT( NSPLIT-1 ), 2 )
      END IF
      PIVMIN = PIVMIN*SAFEMN
*
*     Convert IL and IU
*
      IF( IRANGE.EQ.3 ) THEN
         ILCVT = MAX( IL, (ODDCT+1)/2+1) + N/2
         IUCVT = IU + N/2
         IF( IUCVT.LT.ILCVT ) THEN
            GO TO 160
         END IF
      END IF
*
*     Compute Gershgorin interval for entire (split) matrix
*     and use it as the initial interval
*
      GU = ZERO
      GL = ZERO
      TMP1 = ZERO
*
      DO 20 J = 1, N - 1
         TMP2 = SQRT( WORK( J ) )
         GU = MAX( GU, TMP1+TMP2 )
         GL = MIN( GL, -TMP1-TMP2 )
         TMP1 = TMP2
   20 CONTINUE
*
      GU = MAX( GU, TMP1 )
      GL = MIN( GL, -TMP1 )
      TNORM = MAX( ABS( GL ), ABS( GU ) )
      GL = GL - FUDGE*TNORM*ULP*DBLE( N ) - FUDGE*TWO*PIVMIN
      GU = GU + FUDGE*TNORM*ULP*DBLE( N ) + FUDGE*PIVMIN
*
*     Compute Interval and ATOLI
*
      IF( IRANGE.EQ.3 ) THEN
*
*        RANGE='I': Compute the interval containing eigenvalues
*                   IL through IU.
*
*        Compute Iteration parameters
*
         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) /
     $           LOG( TWO ) ) + 2
         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
*
         WORK( N+1 ) = GL
         WORK( N+2 ) = GL
         WORK( N+3 ) = GU
         WORK( N+4 ) = GU
         WORK( N+5 ) = GL
         WORK( N+6 ) = GU
         IWORK( 1 ) = -1
         IWORK( 2 ) = -1
         IWORK( 3 ) = N + 1
         IWORK( 4 ) = N + 1
         IWORK( 5 ) = ILCVT - 1
         IWORK( 6 ) = IUCVT
*
         CALL DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN,
     $                WORK( 4*N+1 ), E, WORK, IWORK( 5 ),
     $                WORK( N+1 ), WORK( N+5 ), IOUT, IWORK, W,
     $                IBLOCK, IINFO )
*
         IF( IWORK( 6 ).EQ.IUCVT ) THEN
            WL = WORK( N+1 )
            WLU = WORK( N+3 )
            NWL = IWORK( 1 )
            WU = WORK( N+4 )
            WUL = WORK( N+2 )
            NWU = IWORK( 4 )
         ELSE
            WL = WORK( N+2 )
            WLU = WORK( N+4 )
            NWL = IWORK( 2 )
            WU = WORK( N+3 )
            WUL = WORK( N+1 )
            NWU = IWORK( 3 )
         END IF
*
         IF( NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N ) THEN
            INFO = 4
            RETURN
         END IF
      ELSE
*
*        RANGE='A' or 'V' -- Set ATOLI
*
         TNORM = MAX( ABS( E( 1 ) ), ABS( E( N-1 ) ) )
*
         DO 30 J = 2, N - 1
            TNORM = MAX( TNORM, ABS( E( J-1 ) )+ABS( E( J ) ) )
   30    CONTINUE
*
         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
*
         IF( IRANGE.EQ.1 ) THEN
            IRANGE = 2
            WL = ZERO
            WU = GU
         ELSEIF( IRANGE.EQ.2 .AND. VL.LT.ZERO ) THEN
            WL = ZERO
            WU = VU
         ELSEIF( IRANGE.EQ.2 .AND. VL.GE.ZERO ) THEN
            WL = VL
            WU = VU
         ELSE
            WL = ZERO
            WU = ZERO
         END IF
      END IF
*
*     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
*     NWL accumulates the number of eigenvalues .le. WL,
*     NWU accumulates the number of eigenvalues .le. WU
*
      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0
*
      DO 70 JB = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT( JB )
         IN = IEND - IOFF
*
         IF( IN.EQ.1 ) THEN
*
*           Special Case -- IN=1
*
            IF( WL.GE.-PIVMIN )
     $         NWL = NWL + 1
            IF( WU.GE.-PIVMIN )
     $         NWU = NWU + 1
            IF( WL.LT.-PIVMIN .AND. WU.GE.-PIVMIN ) THEN
               M = M + 1
               W( M ) = ZERO
               IBLOCK( M ) = JB
            END IF
         ELSE
*
*           General Case -- IN > 1
*
*           Compute Gershgorin Interval
*           and use it as the initial interval
*
            GU = ZERO
            GL = ZERO
            TMP1 = ZERO
*
            DO 40 J = IBEGIN, IEND - 1
               TMP2 = ABS( E( J ) )
               GU = MAX( GU, TMP1+TMP2 )
               GL = MIN( GL, -TMP1-TMP2 )
               TMP1 = TMP2
   40       CONTINUE
*
            GU = MAX( GU, TMP1 )
            GL = MIN( GL, -TMP1 )
            BNORM = MAX( ABS( GL ), ABS( GU ) )
            GL = GL - FUDGE*BNORM*ULP*DBLE( IN ) - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*DBLE( IN ) + FUDGE*PIVMIN
*
*           Compute ATOLI for the current submatrix
*
            IF( ABSTOL.LE.ZERO ) THEN
               ATOLI = ULP*MAX( ABS( GL ), ABS( GU ) )
            ELSE
               ATOLI = ABSTOL
            END IF
*
            IF( GU.LT.WL ) THEN
               NWL = NWL + IN
               NWU = NWU + IN
               GO TO 70
            END IF
            GL = MAX( GL, WL )
            GU = MIN( GU, WU )
            IF( GL.GE.GU )
     $         GO TO 70
*
*           Set Up Initial Interval
*
            WORK( N+1 ) = GL
            WORK( N+IN+1 ) = GU
            CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
     $                   WORK( 4*N+IBEGIN ), E( IBEGIN ),
     $                   WORK( IBEGIN ), IDUMMA, WORK( N+1 ),
     $                   WORK( N+2*IN+1 ), IM, IWORK, W( M+1 ),
     $                   IBLOCK( M+1 ), IINFO )
*
            NWL = NWL + IWORK( 1 )
            NWU = NWU + IWORK( IN+1 )
            IWOFF = M - IWORK( 1 )
*
*           Compute Eigenvalues
*
            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) /
     $              LOG( TWO ) ) + 2
            CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI,
     $                   PIVMIN, WORK( 4*N+IBEGIN ), E( IBEGIN ),
     $                   WORK( IBEGIN ), IDUMMA, WORK( N+1 ),
     $                   WORK( N+2*IN+1 ), IOUT, IWORK, W( M+1 ),
     $                   IBLOCK( M+1 ), IINFO )
*
*           Copy Eigenvalues Into W and IBLOCK
*           Use -JB for block number for unconverged eigenvalues.
*
            DO 60 J = 1, IOUT
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) )
*
*              Flag non-convergence.
*
               IF( J.GT.IOUT-IINFO ) THEN
                  NCNVRG = .TRUE.
                  IB = -JB
               ELSE
                  IB = JB
               END IF
               DO 50 JE = IWORK( J ) + 1 + IWOFF,
     $                 IWORK( J+IN ) + IWOFF
                  IF (TMP1.GE.ZERO) THEN
                     W( JE ) = TMP1
                  ELSE
                     W( JE ) = -TMP1
                  END IF
                  IBLOCK( JE ) = IB
   50          CONTINUE
   60       CONTINUE
*
            M = M + IM
         END IF
   70 CONTINUE
*
*     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
*     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
*
      IF( IRANGE.EQ.3 ) THEN
         IM = 0
         IDISCL = ILCVT - 1 - NWL
         IDISCU = NWU - IUCVT
*
         IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
            DO 80 JE = 1, M
               IF( W( JE ).LE.WLU .AND. IDISCL.GT.0 ) THEN
                  IDISCL = IDISCL - 1
               ELSE IF( W( JE ).GE.WUL .AND. IDISCU.GT.0 ) THEN
                  IDISCU = IDISCU - 1
               ELSE
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
   80       CONTINUE
            M = IM
         END IF
         IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
*
*           Code to deal with effects of bad arithmetic:
*           Some low eigenvalues to be discarded are not in (WL,WLU],
*           or high eigenvalues to be discarded are not in (WUL,WU]
*           so just kill off the smallest IDISCL/largest IDISCU
*           eigenvalues, by simply finding the smallest/largest
*           eigenvalue(s).
*
*           (If N(w) is monotone non-decreasing, this should never
*               happen.)
*
            IF( IDISCL.GT.0 ) THEN
               WKILL = WU
               DO 100 JDISC = 1, IDISCL
                  IW = 0
                  DO 90 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.
     $                   ( W( JE ).LT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
   90             CONTINUE
                  IBLOCK( IW ) = 0
  100          CONTINUE
            END IF
            IF( IDISCU.GT.0 ) THEN
*
               WKILL = WL
               DO 120 JDISC = 1, IDISCU
                  IW = 0
                  DO 110 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.
     $                   ( W( JE ).GT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
  110             CONTINUE
                  IBLOCK( IW ) = 0
  120          CONTINUE
            END IF
            IM = 0
            DO 130 JE = 1, M
               IF( IBLOCK( JE ).NE.0 ) THEN
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
  130       CONTINUE
            M = IM
         END IF
         IF( IDISCL.LT.0 .OR. IDISCU.LT.0 ) THEN
            TOOFEW = .TRUE.
         END IF
      END IF
*
*     If ORDER='B', sort the eigenvalues from largest to smallest for each block
*     If ORDER='E', sort the eigenvalues from largest to smallest entirely
*
      DO 150 JE = 1, M - 1
         IE = 0
         TMP1 = W( JE )
         DO 140 J = JE + 1, M
            IF( W( J ).GE.TMP1 ) THEN
               IE = J
               TMP1 = W( J )
            END IF
  140       CONTINUE
*
         IF( IE.NE.0 ) THEN
            ITMP1 = IBLOCK( IE )
            W( IE ) = W( JE )
            IBLOCK( IE ) = IBLOCK( JE )
            W( JE ) = TMP1
            IBLOCK( JE ) = ITMP1
         END IF
  150 CONTINUE
  160 CONTINUE
*
*     Expand M, W, and IBLOCK to standard form
*
      IF( IRANGE.EQ.2 .AND. WL.EQ.ZERO ) THEN
         DO 170 J = 1, M
            WORK( 4*N+2*J-1 ) = W(J)
            WORK( 4*N+2*J ) = ZERO
            IWORK( 2*J-1 ) = IBLOCK(J)
            IWORK( 2*J ) = IBLOCK(J)
  170    CONTINUE
         ITMP1 = 0
         IEND = 0
         DO 180 JB = 1, NSPLIT
            IOFF = IEND
            IBEGIN = IOFF + 1
            IEND = ISPLIT( JB )
            IN = IEND - IOFF
            IF( MOD( IN, 2 ).NE.0 ) THEN
               ITMP1 = ITMP1 + 1
               WORK( 4*N+2*M+ITMP1 ) = ZERO
               IWORK( 2*M+ITMP1 ) = JB
            END IF
  180    CONTINUE
         M = 2*M+ITMP1
      ELSEIF( IRANGE.EQ.3 .AND. IUCVT.LT.ILCVT ) THEN
         IF( IL.EQ.1 .AND. MOD( N, 2 ).NE.0 ) THEN
            M = 2*( IU-IL )+1
         ELSE
            M = 2*( IU-IL+1 )
         END IF
         ITMP1 = M
         IEND = 0
         DO 190 JB = 1, NSPLIT
            IOFF = IEND
            IBEGIN = IOFF + 1
            IEND = ISPLIT( JB )
            IN = IEND - IOFF
            IF( MOD( IN, 2 ).NE.0 .AND. ITMP1.GE.1 ) THEN
               WORK( 4*N+M-ITMP1+1 ) = ZERO
               IWORK( M-ITMP1+1 ) = JB
               ITMP1 = ITMP1 - 1
            END IF
  190    CONTINUE
      ELSEIF( IRANGE.EQ.3 .AND. (ODDCT+1)/2-IL.GE.0 ) THEN
         DO 200 J = 1, M
            WORK( 4*N+2*J-1 ) = W(J)
            WORK( 4*N+2*J ) = ZERO
            IWORK( 2*J-1 ) = IBLOCK(J)
            IWORK( 2*J ) = IBLOCK(J)
  200    CONTINUE
         IF( IL.EQ.1 .AND. MOD( N, 2 ).NE.0 ) THEN
            ITMP2 = 2*( (ODDCT+1)/2-IL )+1
         ELSE
            ITMP2 = 2*( (ODDCT+1)/2-IL+1 )
         END IF
         ITMP1 = 0
         IEND = 0
         DO 210 JB = 1, NSPLIT
            IOFF = IEND
            IBEGIN = IOFF + 1
            IEND = ISPLIT( JB )
            IN = IEND - IOFF
            IF( MOD( IN, 2 ).NE.0 .AND. ITMP1.LT.ITMP2 ) THEN
               ITMP1 = ITMP1 + 1
               WORK( 4*N+2*M+ITMP1 ) = ZERO
               IWORK( 2*M+ITMP1 ) = JB
            END IF
  210    CONTINUE
         M = 2*M+ITMP1
      ELSE
         DO 220 J = 1, M
            WORK( 4*N+2*J-1 ) = W(J)
            WORK( 4*N+2*J ) = ZERO
            IWORK( 2*J-1 ) = IBLOCK(J)
            IWORK( 2*J ) = IBLOCK(J)
  220    CONTINUE
         M = 2*M
      END IF
*
      IF ( IORDER.EQ.2 .AND. NSPLIT.GT.1 ) THEN
         ITMP1 = 1
         DO JB = 1, NSPLIT
            DO J = 1, M
               IF ( IWORK( J ).EQ.JB .AND. ITMP1.LE.M ) THEN
                  W( ITMP1 ) = WORK( 4*N+J )
                  IBLOCK( ITMP1 ) = IWORK( J )
                  ITMP1 = ITMP1+1
               END IF
            END DO
         END DO
      ELSE
         DO J = 1, M
            W( J ) = WORK( 4*N+J )
            IBLOCK( J ) = IWORK( J )
         END DO
      END IF
*
      INFO = 0
      IF( NCNVRG )
     $   INFO = INFO + 1
      IF( TOOFEW )
     $   INFO = INFO + 2
      RETURN
*
*     End of DKTEBZ
*
      END
