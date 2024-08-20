*> \brief \b DLATMS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLATMS( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX,
*                          KL, KU, PACK, A, LDA, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIST, PACK, SYM
*       INTEGER            INFO, KL, KU, LDA, M, MODE, N
*       DOUBLE PRECISION   COND, DMAX
*       ..
*       .. Array Arguments ..
*       INTEGER            ISEED( 4 )
*       DOUBLE PRECISION   A( LDA, * ), D( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DLATMS generates random matrices with specified singular values
*>    (or symmetric/hermitian/skew-symmetric with specified eigenvalues)
*>    for testing LAPACK programs.
*>
*>    DLATMS operates by applying the following sequence of
*>    operations:
*>
*>      Set the diagonal to D, where D may be input or
*>         computed according to MODE, COND, DMAX, and SYM
*>         as described below.
*>
*>      Generate a matrix with the appropriate band structure, by one
*>         of two methods:
*>
*>      Method A:
*>          Generate a dense M x N matrix by multiplying D on the left
*>              and the right by random unitary matrices, then:
*>
*>          Reduce the bandwidth according to KL and KU, using
*>          Householder transformations.
*>
*>      Method B:
*>          Convert the bandwidth-0 (i.e., diagonal) matrix to a
*>              bandwidth-1 matrix using Givens rotations, "chasing"
*>              out-of-band elements back, much as in QR; then
*>              convert the bandwidth-1 to a bandwidth-2 matrix, etc.
*>              Note that for reasonably small bandwidths (relative to
*>              M and N) this requires less storage, as a dense matrix
*>              is not generated.  Also, for symmetric matrices, only
*>              one triangle is generated.
*>
*>      Method A is chosen if the bandwidth is a large fraction of the
*>          order of the matrix, and LDA is at least M (so a dense
*>          matrix can be stored.)  Method B is chosen if the bandwidth
*>          is small (< 1/2 N for symmetric, < .3 N+M for
*>          non-symmetric), or LDA is less than M and not less than the
*>          bandwidth.
*>
*>      Pack the matrix if desired. Options specified by PACK are:
*>         no packing
*>         zero out upper half (if symmetric/skew-symmetric)
*>         zero out lower half (if symmetric/skew-symmetric)
*>         store the upper half columnwise (if symmetric or upper
*>               triangular)
*>         store the lower half columnwise (if symmetric or lower
*>               triangular)
*>         store the lower triangle in banded format (if symmetric
*>               or lower triangular)
*>         store the upper triangle in banded format (if symmetric
*>               or upper triangular)
*>         store the entire matrix in banded format
*>      If Method B is chosen, and band format is specified, then the
*>         matrix will be generated in the band format, so no repacking
*>         will be necessary.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           The number of rows of A. Not modified.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           The number of columns of A. Not modified.
*> \endverbatim
*>
*> \param[in] DIST
*> \verbatim
*>          DIST is CHARACTER*1
*>           On entry, DIST specifies the type of distribution to be used
*>           to generate the random eigen-/singular values.
*>           'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )
*>           'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric/skew-symmetric )
*>           'N' => NORMAL( 0, 1 )   ( 'N' for normal )
*>           Not modified.
*> \endverbatim
*>
*> \param[in,out] ISEED
*> \verbatim
*>          ISEED is INTEGER array, dimension ( 4 )
*>           On entry ISEED specifies the seed of the random number
*>           generator. They should lie between 0 and 4095 inclusive,
*>           and ISEED(4) should be odd. The random number generator
*>           uses a linear congruential sequence limited to small
*>           integers, and so should produce machine independent
*>           random numbers. The values of ISEED are changed on
*>           exit, and can be used in the next call to DLATMS
*>           to continue the same random number sequence.
*>           Changed on exit.
*> \endverbatim
*>
*> \param[in] SYM
*> \verbatim
*>          SYM is CHARACTER*1
*>           If SYM='S' or 'H', the generated matrix is symmetric, with
*>             eigenvalues specified by D, COND, MODE, and DMAX; they
*>             may be positive, negative, or zero.
*>           If SYM='K', the generated matrix is skew-symmetric, with
*>             eigenvalues specified by D, COND, MODE, and DMAX; they
*>             may be positive, negative, or zero.
*>           If SYM='P', the generated matrix is symmetric, with
*>             eigenvalues (= singular values) specified by D, COND,
*>             MODE, and DMAX; they will not be negative.
*>           If SYM='N', the generated matrix is nonsymmetric, with
*>             singular values specified by D, COND, MODE, and DMAX;
*>             they will not be negative.
*>           Not modified.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension ( MIN( M , N ) )
*>           This array is used to specify the singular values or
*>           eigenvalues of A (see SYM, above.)  If MODE=0, then D is
*>           assumed to contain the singular/eigenvalues, otherwise
*>           they will be computed according to MODE, COND, and DMAX,
*>           and placed in D.
*>           Modified if MODE is nonzero.
*> \endverbatim
*>
*> \param[in] MODE
*> \verbatim
*>          MODE is INTEGER
*>           On entry this describes how the singular/eigenvalues are to
*>           be specified:
*>           MODE = 0 means use D as input
*>           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
*>           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
*>           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
*>           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
*>           MODE = 5 sets D to random numbers in the range
*>                    ( 1/COND , 1 ) such that their logarithms
*>                    are uniformly distributed.
*>           MODE = 6 set D to random numbers from same distribution
*>                    as the rest of the matrix.
*>           MODE < 0 has the same meaning as ABS(MODE), except that
*>              the order of the elements of D is reversed.
*>           Thus if MODE is positive, D has entries ranging from
*>              1 to 1/COND, if negative, from 1/COND to 1,
*>           If SYM='S' or 'H', and MODE is neither 0, 6, nor -6, then
*>              the elements of D will also be multiplied by a random
*>              sign (i.e., +1 or -1.)
*>           Not modified.
*> \endverbatim
*>
*> \param[in] COND
*> \verbatim
*>          COND is DOUBLE PRECISION
*>           On entry, this is used as described under MODE above.
*>           If used, it must be >= 1. Not modified.
*> \endverbatim
*>
*> \param[in] DMAX
*> \verbatim
*>          DMAX is DOUBLE PRECISION
*>           If MODE is neither -6, 0 nor 6, the contents of D, as
*>           computed according to MODE and COND, will be scaled by
*>           DMAX / max(abs(D(i))); thus, the maximum absolute eigen- or
*>           singular value (which is to say the norm) will be abs(DMAX).
*>           Note that DMAX need not be positive: if DMAX is negative
*>           (or zero), D will be scaled by a negative number (or zero).
*>           Not modified.
*> \endverbatim
*>
*> \param[in] KL
*> \verbatim
*>          KL is INTEGER
*>           This specifies the lower bandwidth of the  matrix. For
*>           example, KL=0 implies upper triangular, KL=1 implies upper
*>           Hessenberg, and KL being at least M-1 means that the matrix
*>           has full lower bandwidth.  KL must equal KU if the matrix
*>           is symmetric/skew-symmetric.
*>           Not modified.
*> \endverbatim
*>
*> \param[in] KU
*> \verbatim
*>          KU is INTEGER
*>           This specifies the upper bandwidth of the  matrix. For
*>           example, KU=0 implies lower triangular, KU=1 implies lower
*>           Hessenberg, and KU being at least N-1 means that the matrix
*>           has full upper bandwidth.  KL must equal KU if the matrix
*>           is symmetric/skew-symmetric.
*>           Not modified.
*> \endverbatim
*>
*> \param[in] PACK
*> \verbatim
*>          PACK is CHARACTER*1
*>           This specifies packing of matrix as follows:
*>           'N' => no packing
*>           'U' => zero out all subdiagonal entries (if symmetric/skew-symmetric)
*>           'L' => zero out all superdiagonal entries (if symmetric/skew-symmetric)
*>           'C' => store the upper triangle columnwise
*>                  (only if the matrix is symmetric or upper triangular)
*>           'R' => store the lower triangle columnwise
*>                  (only if the matrix is symmetric or lower triangular)
*>           'B' => store the lower triangle in band storage scheme
*>                  (only if matrix symmetric or lower triangular)
*>           'Q' => store the upper triangle in band storage scheme
*>                  (only if matrix symmetric or upper triangular)
*>           'Z' => store the entire matrix in band storage scheme
*>                      (pivoting can be provided for by using this
*>                      option to store A in the trailing rows of
*>                      the allocated storage)
*>
*>           Using these options, the various LAPACK packed and banded
*>           storage schemes can be obtained:
*>           GB               - use 'Z'
*>           PB, SB or TB     - use 'B' or 'Q'
*>           PP, SP or TP     - use 'C' or 'R'
*>
*>           If two calls to DLATMS differ only in the PACK parameter,
*>           they will generate mathematically equivalent matrices.
*>           Not modified.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, N )
*>           On exit A is the desired test matrix.  A is first generated
*>           in full (unpacked) form, and then packed, if so specified
*>           by PACK.  Thus, the first M elements of the first N
*>           columns will always be modified.  If PACK specifies a
*>           packed or banded storage scheme, all LDA elements of the
*>           first N columns will be modified; the elements of the
*>           array which do not correspond to elements of the generated
*>           matrix are set to zero.
*>           Modified.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           LDA specifies the first dimension of A as declared in the
*>           calling program.  If PACK='N', 'U', 'L', 'C', or 'R', then
*>           LDA must be at least M.  If PACK='B' or 'Q', then LDA must
*>           be at least MIN( KL, M-1) (which is equal to MIN(KU,N-1)).
*>           If PACK='Z', LDA must be large enough to hold the packed
*>           array: MIN( KU, N-1) + MIN( KL, M-1) + 1.
*>           Not modified.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension ( 3*MAX( N , M ) )
*>           Workspace.
*>           Modified.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           Error code.  On exit, INFO will be set to one of the
*>           following values:
*>             0 => normal return
*>            -1 => M negative or unequal to N and SYM='S', 'H', 'K', or 'P'
*>            -2 => N negative
*>            -3 => DIST illegal string
*>            -5 => SYM illegal string
*>            -7 => MODE not in range -6 to 6
*>            -8 => COND less than 1.0, and MODE neither -6, 0 nor 6
*>           -10 => KL negative
*>           -11 => KU negative, or SYM='S' or 'H' and KU not equal to KL
*>           -12 => PACK illegal string, or PACK='U' or 'L', and SYM='N';
*>                  or PACK='C' or 'Q' and SYM='N' and KL is not zero;
*>                  or PACK='R' or 'B' and SYM='N' and KU is not zero;
*>                  or PACK='U', 'L', 'C', 'R', 'B', or 'Q', and M is not
*>                  N.
*>           -14 => LDA is less than M, or PACK='Z' and LDA is less than
*>                  MIN(KU,N-1) + MIN(KL,M-1) + 1.
*>            1  => Error return from DLATM1
*>            2  => Cannot scale to DMAX (max. sing. value is 0)
*>            3  => Error return from DLAGGE or SLAGSY
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
*> \ingroup double_matgen
*
*  =====================================================================
      SUBROUTINE DLATMS( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX,
     $                   KL, KU, PACK, A, LDA, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIST, PACK, SYM
      INTEGER            INFO, KL, KU, LDA, M, MODE, N
      DOUBLE PRECISION   COND, DMAX
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), D( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWOPI
      PARAMETER  ( TWOPI = 6.28318530717958647692528676655900576839D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            GIVENS, ILEXTR, ILTEMP, TOPDWN
      INTEGER            I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA,
     $                   IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2,
     $                   IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH,
     $                   JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC,
     $                   UUB, MNMINNEW
      DOUBLE PRECISION   ALPHA, ANGLE, C, DUMMY, EXTRA, S, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLARND
      EXTERNAL           LSAME, DLARND
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAGGE, DLAGSY, DLAGKY,
     $                   DLAROT, DLARTG, DLASET,
     $                   DLATM1, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, COS, DBLE, MAX, MIN, MOD, SIN
*     ..
*     .. Executable Statements ..
*
*     1)      Decode and Test the input parameters.
*             Initialize flags & seed.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Decode DIST
*
      IF( LSAME( DIST, 'U' ) ) THEN
         IDIST = 1
      ELSE IF( LSAME( DIST, 'S' ) ) THEN
         IDIST = 2
      ELSE IF( LSAME( DIST, 'N' ) ) THEN
         IDIST = 3
      ELSE
         IDIST = -1
      END IF
*
*     Decode SYM
*
      IF( LSAME( SYM, 'N' ) ) THEN
         ISYM = 1
         IRSIGN = 0
      ELSE IF( LSAME( SYM, 'P' ) ) THEN
         ISYM = 2
         IRSIGN = 0
      ELSE IF( LSAME( SYM, 'S' ) ) THEN
         ISYM = 2
         IRSIGN = 1
      ELSE IF( LSAME( SYM, 'H' ) ) THEN
         ISYM = 2
         IRSIGN = 1
      ELSE IF( LSAME( SYM, 'K' ) ) THEN
         ISYM = 3
         IRSIGN = 1
      ELSE
         ISYM = -1
      END IF
*
*     Decode PACK
*
      ISYMPK = 0
      IF( LSAME( PACK, 'N' ) ) THEN
         IPACK = 0
      ELSE IF( LSAME( PACK, 'U' ) ) THEN
         IPACK = 1
         ISYMPK = 1
      ELSE IF( LSAME( PACK, 'L' ) ) THEN
         IPACK = 2
         ISYMPK = 1
      ELSE IF( LSAME( PACK, 'C' ) ) THEN
         IPACK = 3
         ISYMPK = 2
      ELSE IF( LSAME( PACK, 'R' ) ) THEN
         IPACK = 4
         ISYMPK = 3
      ELSE IF( LSAME( PACK, 'B' ) ) THEN
         IPACK = 5
         ISYMPK = 3
      ELSE IF( LSAME( PACK, 'Q' ) ) THEN
         IPACK = 6
         ISYMPK = 2
      ELSE IF( LSAME( PACK, 'Z' ) ) THEN
         IPACK = 7
      ELSE
         IPACK = -1
      END IF
*
*     Set certain internal parameters
*
      MNMIN = MIN( M, N )
      LLB = MIN( KL, M-1 )
      UUB = MIN( KU, N-1 )
      MR = MIN( M, N+LLB )
      NC = MIN( N, M+UUB )
*
      IF( IPACK.EQ.5 .OR. IPACK.EQ.6 ) THEN
         MINLDA = UUB + 1
      ELSE IF( IPACK.EQ.7 ) THEN
         MINLDA = LLB + UUB + 1
      ELSE
         MINLDA = M
      END IF
*
*     Use Givens rotation method if bandwidth small enough,
*     or if LDA is too small to store the matrix unpacked.
*
      GIVENS = .FALSE.
      IF( ISYM.EQ.1 ) THEN
         IF( DBLE( LLB+UUB ).LT.0.3D0*DBLE( MAX( 1, MR+NC ) ) )
     $      GIVENS = .TRUE.
      ELSE
         IF( 2*LLB.LT.M )
     $      GIVENS = .TRUE.
      END IF
      IF( LDA.LT.M .AND. LDA.GE.MINLDA )
     $   GIVENS = .TRUE.
      IF( ISYM.EQ.3 ) THEN
         GIVENS = .FALSE.
      END IF 
*
*     Set INFO if an error
*
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.NE.N .AND. ISYM.NE.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( IDIST.EQ.-1 ) THEN
         INFO = -3
      ELSE IF( ISYM.EQ.-1 ) THEN
         INFO = -5
      ELSE IF( ABS( MODE ).GT.6 ) THEN
         INFO = -7
      ELSE IF( ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) .AND. COND.LT.ONE )
     $          THEN
         INFO = -8
      ELSE IF( KL.LT.0 ) THEN
         INFO = -10
      ELSE IF( KU.LT.0 .OR. ( ISYM.NE.1 .AND. KL.NE.KU ) ) THEN
         INFO = -11
      ELSE IF( IPACK.EQ.-1 .OR. ( ISYMPK.EQ.1 .AND. ISYM.EQ.1 ) .OR.
     $         ( ISYMPK.EQ.2 .AND. ISYM.EQ.1 .AND. KL.GT.0 ) .OR.
     $         ( ISYMPK.EQ.3 .AND. ISYM.EQ.1 .AND. KU.GT.0 ) .OR.
     $         ( ISYMPK.NE.0 .AND. M.NE.N ) ) THEN
         INFO = -12
      ELSE IF( LDA.LT.MAX( 1, MINLDA ) ) THEN
         INFO = -14
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATMS', -INFO )
         RETURN
      END IF
*
*     Initialize random number generator
*
      DO 10 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   10 CONTINUE
*
      IF( MOD( ISEED( 4 ), 2 ).NE.1 )
     $   ISEED( 4 ) = ISEED( 4 ) + 1
*
*     2)      Set up D  if indicated.
*
*             Compute D according to COND and MODE
*
      IF( ISYM.EQ.3 ) THEN
         MNMINNEW = MNMIN / 2
      ELSE
         MNMINNEW = MNMIN
      END IF
      CALL DLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D,
     $             MNMINNEW, IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      IF( ISYM.EQ.3 .AND. MNMIN.EQ.1 ) THEN
         D(1) = ONE
      END IF
*
*     Choose Top-Down if D is (apparently) increasing,
*     Bottom-Up if D is (apparently) decreasing.
*
      IF( ABS( D( 1 ) ).LE.ABS( D( MNMINNEW ) ) ) THEN
         TOPDWN = .TRUE.
      ELSE
         TOPDWN = .FALSE.
      END IF
*
      IF( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) THEN
*
*        Scale by DMAX
*
         TEMP = ABS( D( 1 ) )
         DO 20 I = 2, MNMINNEW
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   20    CONTINUE
*
         IF( TEMP.GT.ZERO ) THEN
            ALPHA = DMAX / TEMP
         ELSE
            INFO = 2
            RETURN
         END IF
*
         CALL DSCAL( MNMINNEW, ALPHA, D, 1 )
*
      END IF
*
      IF( ISYM.EQ.3 ) THEN
         DO I = MNMIN, 2*MNMINNEW + 1, -1
            D(I) = ZERO
         END DO
         DO I = MNMINNEW, 1, -1
            D(2*I - 1) = D(I)
            D(2*I) = ZERO
         END DO
      END IF
*
*     3)      Generate Banded Matrix using Givens rotations.
*             Also the special case of UUB=LLB=0
*
*               Compute Addressing constants to cover all
*               storage formats.  Whether GE, SY, GB, or SB,
*               upper or lower triangle or both,
*               the (i,j)-th element is in
*               A( i - ISKEW*j + IOFFST, j )
*
      IF( IPACK.GT.4 ) THEN
         ILDA = LDA - 1
         ISKEW = 1
         IF( IPACK.GT.5 ) THEN
            IOFFST = UUB + 1
         ELSE
            IOFFST = 1
         END IF
      ELSE
         ILDA = LDA
         ISKEW = 0
         IOFFST = 0
      END IF
*
*     IPACKG is the format that the matrix is generated in. If this is
*     different from IPACK, then the matrix must be repacked at the
*     end.  It also signals how to compute the norm, for scaling.
*
      IPACKG = 0
      CALL DLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
*
*     Diagonal Matrix -- We are done, unless it
*     is to be stored SP/PP/TP (PACK='R' or 'C')
*
      IF( LLB.EQ.0 .AND. UUB.EQ.0 ) THEN
         CALL DCOPY( MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 )
         IF( IPACK.LE.2 .OR. IPACK.GE.5 )
     $      IPACKG = IPACK
*
      ELSE IF( GIVENS ) THEN
*
*        Check whether to use Givens rotations,
*        Householder transformations, or nothing.
*
         IF( ISYM.EQ.1 ) THEN
*
*           Non-symmetric -- A = U D V
*
            IF( IPACK.GT.4 ) THEN
               IPACKG = IPACK
            ELSE
               IPACKG = 0
            END IF
*
            CALL DCOPY( MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 )
*
            IF( TOPDWN ) THEN
               JKL = 0
               DO 50 JKU = 1, UUB
*
*                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
*
*                 Last row actually rotated is M
*                 Last column actually rotated is MIN( M+JKU, N )
*
                  DO 40 JR = 1, MIN( M+JKU, N ) + JKL - 1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     ICOL = MAX( 1, JR-JKL )
                     IF( JR.LT.M ) THEN
                        IL = MIN( N, JR+JKU ) + 1 - ICOL
                        CALL DLAROT( .TRUE., JR.GT.JKL, .FALSE., IL,
     $                               C,
     $                               S, A( JR-ISKEW*ICOL+IOFFST, ICOL ),
     $                               ILDA, EXTRA, DUMMY )
                     END IF
*
*                    Chase "EXTRA" back up
*
                     IR = JR
                     IC = ICOL
                     DO 30 JCH = JR - JKL, 1, -JKL - JKU
                        IF( IR.LT.M ) THEN
                           CALL DLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST,
     $                                  IC+1 ), EXTRA, C, S, DUMMY )
                        END IF
                        IROW = MAX( 1, JCH-JKU )
                        IL = IR + 2 - IROW
                        TEMP = ZERO
                        ILTEMP = JCH.GT.JKU
                        CALL DLAROT( .FALSE., ILTEMP, .TRUE., IL, C,
     $                               -S,
     $                               A( IROW-ISKEW*IC+IOFFST, IC ),
     $                               ILDA, TEMP, EXTRA )
                        IF( ILTEMP ) THEN
                           CALL DLARTG( A( IROW+1-ISKEW*( IC+1 )+IOFFST,
     $                                  IC+1 ), TEMP, C, S, DUMMY )
                           ICOL = MAX( 1, JCH-JKU-JKL )
                           IL = IC + 2 - ICOL
                           EXTRA = ZERO
                           CALL DLAROT( .TRUE., JCH.GT.JKU+JKL,
     $                                  .TRUE.,
     $                                  IL, C, -S, A( IROW-ISKEW*ICOL+
     $                                  IOFFST, ICOL ), ILDA, EXTRA,
     $                                  TEMP )
                           IC = ICOL
                           IR = IROW
                        END IF
   30                CONTINUE
   40             CONTINUE
   50          CONTINUE
*
               JKU = UUB
               DO 80 JKL = 1, LLB
*
*                 Transform from bandwidth JKL-1, JKU to JKL, JKU
*
                  DO 70 JC = 1, MIN( N+JKL, M ) + JKU - 1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     IROW = MAX( 1, JC-JKU )
                     IF( JC.LT.N ) THEN
                        IL = MIN( M, JC+JKL ) + 1 - IROW
                        CALL DLAROT( .FALSE., JC.GT.JKU, .FALSE., IL,
     $                               C,
     $                               S, A( IROW-ISKEW*JC+IOFFST, JC ),
     $                               ILDA, EXTRA, DUMMY )
                     END IF
*
*                    Chase "EXTRA" back up
*
                     IC = JC
                     IR = IROW
                     DO 60 JCH = JC - JKU, 1, -JKL - JKU
                        IF( IC.LT.N ) THEN
                           CALL DLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST,
     $                                  IC+1 ), EXTRA, C, S, DUMMY )
                        END IF
                        ICOL = MAX( 1, JCH-JKL )
                        IL = IC + 2 - ICOL
                        TEMP = ZERO
                        ILTEMP = JCH.GT.JKL
                        CALL DLAROT( .TRUE., ILTEMP, .TRUE., IL, C,
     $                               -S,
     $                               A( IR-ISKEW*ICOL+IOFFST, ICOL ),
     $                               ILDA, TEMP, EXTRA )
                        IF( ILTEMP ) THEN
                           CALL DLARTG( A( IR+1-ISKEW*( ICOL+1 )+IOFFST,
     $                                  ICOL+1 ), TEMP, C, S, DUMMY )
                           IROW = MAX( 1, JCH-JKL-JKU )
                           IL = IR + 2 - IROW
                           EXTRA = ZERO
                           CALL DLAROT( .FALSE., JCH.GT.JKL+JKU,
     $                                  .TRUE.,
     $                                  IL, C, -S, A( IROW-ISKEW*ICOL+
     $                                  IOFFST, ICOL ), ILDA, EXTRA,
     $                                  TEMP )
                           IC = ICOL
                           IR = IROW
                        END IF
   60                CONTINUE
   70             CONTINUE
   80          CONTINUE
*
            ELSE
*
*              Bottom-Up -- Start at the bottom right.
*
               JKL = 0
               DO 110 JKU = 1, UUB
*
*                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
*
*                 First row actually rotated is M
*                 First column actually rotated is MIN( M+JKU, N )
*
                  IENDCH = MIN( M, N+JKL ) - 1
                  DO 100 JC = MIN( M+JKU, N ) - 1, 1 - JKL, -1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     IROW = MAX( 1, JC-JKU+1 )
                     IF( JC.GT.0 ) THEN
                        IL = MIN( M, JC+JKL+1 ) + 1 - IROW
                        CALL DLAROT( .FALSE., .FALSE., JC+JKL.LT.M,
     $                               IL,
     $                               C, S, A( IROW-ISKEW*JC+IOFFST,
     $                               JC ), ILDA, DUMMY, EXTRA )
                     END IF
*
*                    Chase "EXTRA" back down
*
                     IC = JC
                     DO 90 JCH = JC + JKL, IENDCH, JKL + JKU
                        ILEXTR = IC.GT.0
                        IF( ILEXTR ) THEN
                           CALL DLARTG( A( JCH-ISKEW*IC+IOFFST, IC ),
     $                                  EXTRA, C, S, DUMMY )
                        END IF
                        IC = MAX( 1, IC )
                        ICOL = MIN( N-1, JCH+JKU )
                        ILTEMP = JCH + JKU.LT.N
                        TEMP = ZERO
                        CALL DLAROT( .TRUE., ILEXTR, ILTEMP,
     $                               ICOL+2-IC,
     $                               C, S, A( JCH-ISKEW*IC+IOFFST, IC ),
     $                               ILDA, EXTRA, TEMP )
                        IF( ILTEMP ) THEN
                           CALL DLARTG( A( JCH-ISKEW*ICOL+IOFFST,
     $                                  ICOL ), TEMP, C, S, DUMMY )
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = ZERO
                           CALL DLAROT( .FALSE., .TRUE.,
     $                                  JCH+JKL+JKU.LE.IENDCH, IL, C, S,
     $                                  A( JCH-ISKEW*ICOL+IOFFST,
     $                                  ICOL ), ILDA, TEMP, EXTRA )
                           IC = ICOL
                        END IF
   90                CONTINUE
  100             CONTINUE
  110          CONTINUE
*
               JKU = UUB
               DO 140 JKL = 1, LLB
*
*                 Transform from bandwidth JKL-1, JKU to JKL, JKU
*
*                 First row actually rotated is MIN( N+JKL, M )
*                 First column actually rotated is N
*
                  IENDCH = MIN( N, M+JKU ) - 1
                  DO 130 JR = MIN( N+JKL, M ) - 1, 1 - JKU, -1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     ICOL = MAX( 1, JR-JKL+1 )
                     IF( JR.GT.0 ) THEN
                        IL = MIN( N, JR+JKU+1 ) + 1 - ICOL
                        CALL DLAROT( .TRUE., .FALSE., JR+JKU.LT.N,
     $                               IL,
     $                               C, S, A( JR-ISKEW*ICOL+IOFFST,
     $                               ICOL ), ILDA, DUMMY, EXTRA )
                     END IF
*
*                    Chase "EXTRA" back down
*
                     IR = JR
                     DO 120 JCH = JR + JKU, IENDCH, JKL + JKU
                        ILEXTR = IR.GT.0
                        IF( ILEXTR ) THEN
                           CALL DLARTG( A( IR-ISKEW*JCH+IOFFST, JCH ),
     $                                  EXTRA, C, S, DUMMY )
                        END IF
                        IR = MAX( 1, IR )
                        IROW = MIN( M-1, JCH+JKL )
                        ILTEMP = JCH + JKL.LT.M
                        TEMP = ZERO
                        CALL DLAROT( .FALSE., ILEXTR, ILTEMP,
     $                               IROW+2-IR,
     $                               C, S, A( IR-ISKEW*JCH+IOFFST,
     $                               JCH ), ILDA, EXTRA, TEMP )
                        IF( ILTEMP ) THEN
                           CALL DLARTG( A( IROW-ISKEW*JCH+IOFFST, JCH ),
     $                                  TEMP, C, S, DUMMY )
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = ZERO
                           CALL DLAROT( .TRUE., .TRUE.,
     $                                  JCH+JKL+JKU.LE.IENDCH, IL, C, S,
     $                                  A( IROW-ISKEW*JCH+IOFFST, JCH ),
     $                                  ILDA, TEMP, EXTRA )
                           IR = IROW
                        END IF
  120                CONTINUE
  130             CONTINUE
  140          CONTINUE
            END IF
*
         ELSE
*
*           Symmetric -- A = U D U'
*
            IPACKG = IPACK
            IOFFG = IOFFST
*
            IF( TOPDWN ) THEN
*
*              Top-Down -- Generate Upper triangle only
*
               IF( IPACK.GE.5 ) THEN
                  IPACKG = 6
                  IOFFG = UUB + 1
               ELSE
                  IPACKG = 1
               END IF
               CALL DCOPY( MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 )
*
               DO 170 K = 1, UUB
                  DO 160 JC = 1, N - 1
                     IROW = MAX( 1, JC-K )
                     IL = MIN( JC+1, K+2 )
                     EXTRA = ZERO
                     TEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 )
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     CALL DLAROT( .FALSE., JC.GT.K, .TRUE., IL, C, S,
     $                            A( IROW-ISKEW*JC+IOFFG, JC ), ILDA,
     $                            EXTRA, TEMP )
                     CALL DLAROT( .TRUE., .TRUE., .FALSE.,
     $                            MIN( K, N-JC )+1, C, S,
     $                            A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA,
     $                            TEMP, DUMMY )
*
*                    Chase EXTRA back up the matrix
*
                     ICOL = JC
                     DO 150 JCH = JC - K, 1, -K
                        CALL DLARTG( A( JCH+1-ISKEW*( ICOL+1 )+IOFFG,
     $                               ICOL+1 ), EXTRA, C, S, DUMMY )
                        TEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 )
                        CALL DLAROT( .TRUE., .TRUE., .TRUE., K+2, C,
     $                               -S,
     $                               A( ( 1-ISKEW )*JCH+IOFFG, JCH ),
     $                               ILDA, TEMP, EXTRA )
                        IROW = MAX( 1, JCH-K )
                        IL = MIN( JCH+1, K+2 )
                        EXTRA = ZERO
                        CALL DLAROT( .FALSE., JCH.GT.K, .TRUE., IL,
     $                               C,
     $                               -S, A( IROW-ISKEW*JCH+IOFFG, JCH ),
     $                               ILDA, EXTRA, TEMP )
                        ICOL = JCH
  150                CONTINUE
  160             CONTINUE
  170          CONTINUE
*
*              If we need lower triangle, copy from upper. Note that
*              the order of copying is chosen to work for 'q' -> 'b'
*
               IF( IPACK.NE.IPACKG .AND. IPACK.NE.3 ) THEN
                  DO 190 JC = 1, N
                     IROW = IOFFST - ISKEW*JC
                     DO 180 JR = JC, MIN( N, JC+UUB )
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
  180                CONTINUE
  190             CONTINUE
                  IF( IPACK.EQ.5 ) THEN
                     DO 210 JC = N - UUB + 1, N
                        DO 200 JR = N + 2 - JC, UUB + 1
                           A( JR, JC ) = ZERO
  200                   CONTINUE
  210                CONTINUE
                  END IF
                  IF( IPACKG.EQ.6 ) THEN
                     IPACKG = IPACK
                  ELSE
                     IPACKG = 0
                  END IF
               END IF
            ELSE
*
*              Bottom-Up -- Generate Lower triangle only
*
               IF( IPACK.GE.5 ) THEN
                  IPACKG = 5
                  IF( IPACK.EQ.6 )
     $               IOFFG = 1
               ELSE
                  IPACKG = 2
               END IF
               CALL DCOPY( MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 )
*
               DO 240 K = 1, UUB
                  DO 230 JC = N - 1, 1, -1
                     IL = MIN( N+1-JC, K+2 )
                     EXTRA = ZERO
                     TEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC )
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = -SIN( ANGLE )
                     CALL DLAROT( .FALSE., .TRUE., N-JC.GT.K, IL, C,
     $                            S,
     $                            A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA,
     $                            TEMP, EXTRA )
                     ICOL = MAX( 1, JC-K+1 )
                     CALL DLAROT( .TRUE., .FALSE., .TRUE., JC+2-ICOL,
     $                            C,
     $                            S, A( JC-ISKEW*ICOL+IOFFG, ICOL ),
     $                            ILDA, DUMMY, TEMP )
*
*                    Chase EXTRA back down the matrix
*
                     ICOL = JC
                     DO 220 JCH = JC + K, N - 1, K
                        CALL DLARTG( A( JCH-ISKEW*ICOL+IOFFG, ICOL ),
     $                               EXTRA, C, S, DUMMY )
                        TEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH )
                        CALL DLAROT( .TRUE., .TRUE., .TRUE., K+2, C,
     $                               S,
     $                               A( JCH-ISKEW*ICOL+IOFFG, ICOL ),
     $                               ILDA, EXTRA, TEMP )
                        IL = MIN( N+1-JCH, K+2 )
                        EXTRA = ZERO
                        CALL DLAROT( .FALSE., .TRUE., N-JCH.GT.K, IL,
     $                               C,
     $                               S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ),
     $                               ILDA, TEMP, EXTRA )
                        ICOL = JCH
  220                CONTINUE
  230             CONTINUE
  240          CONTINUE
*
*              If we need upper triangle, copy from lower. Note that
*              the order of copying is chosen to work for 'b' -> 'q'
*
               IF( IPACK.NE.IPACKG .AND. IPACK.NE.4 ) THEN
                  DO 260 JC = N, 1, -1
                     IROW = IOFFST - ISKEW*JC
                     DO 250 JR = JC, MAX( 1, JC-UUB ), -1
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
  250                CONTINUE
  260             CONTINUE
                  IF( IPACK.EQ.6 ) THEN
                     DO 280 JC = 1, UUB
                        DO 270 JR = 1, UUB + 1 - JC
                           A( JR, JC ) = ZERO
  270                   CONTINUE
  280                CONTINUE
                  END IF
                  IF( IPACKG.EQ.5 ) THEN
                     IPACKG = IPACK
                  ELSE
                     IPACKG = 0
                  END IF
               END IF
            END IF
         END IF
*
      ELSE
*
*        4)      Generate Banded Matrix by first
*                Rotating by random Unitary matrices,
*                then reducing the bandwidth using Householder
*                transformations.
*
*                Note: we should get here only if LDA .ge. N
*
         IF( ISYM.EQ.1 ) THEN
*
*           Non-symmetric -- A = U D V
*
            CALL DLAGGE( MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK,
     $                   IINFO )
         ELSEIF( ISYM.EQ.2 ) THEN
*
*           Symmetric -- A = U D U'
*
            CALL DLAGSY( M, LLB, D, A, LDA, ISEED, WORK, IINFO )
*
         ELSE
*
*           Skew-symmetric -- A = U D U'
*
            CALL DLAGKY( M, LLB, D, A, LDA, ISEED, WORK, IINFO )
*
         END IF
         IF( IINFO.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
      END IF
*
*     5)      Pack the matrix
*
      IF( IPACK.NE.IPACKG ) THEN
         IF( IPACK.EQ.1 ) THEN
*
*           'U' -- Upper triangular, not packed
*
            DO 300 J = 1, M
               DO 290 I = J + 1, M
                  A( I, J ) = ZERO
  290          CONTINUE
  300       CONTINUE
*
         ELSE IF( IPACK.EQ.2 ) THEN
*
*           'L' -- Lower triangular, not packed
*
            DO 320 J = 2, M
               DO 310 I = 1, J - 1
                  A( I, J ) = ZERO
  310          CONTINUE
  320       CONTINUE
*
         ELSE IF( IPACK.EQ.3 ) THEN
*
*           'C' -- Upper triangle packed Columnwise.
*
            ICOL = 1
            IROW = 0
            DO 340 J = 1, M
               DO 330 I = 1, J
                  IROW = IROW + 1
                  IF( IROW.GT.LDA ) THEN
                     IROW = 1
                     ICOL = ICOL + 1
                  END IF
                  A( IROW, ICOL ) = A( I, J )
  330          CONTINUE
  340       CONTINUE
*
         ELSE IF( IPACK.EQ.4 ) THEN
*
*           'R' -- Lower triangle packed Columnwise.
*
            ICOL = 1
            IROW = 0
            DO 360 J = 1, M
               DO 350 I = J, M
                  IROW = IROW + 1
                  IF( IROW.GT.LDA ) THEN
                     IROW = 1
                     ICOL = ICOL + 1
                  END IF
                  A( IROW, ICOL ) = A( I, J )
  350          CONTINUE
  360       CONTINUE
*
         ELSE IF( IPACK.GE.5 ) THEN
*
*           'B' -- The lower triangle is packed as a band matrix.
*           'Q' -- The upper triangle is packed as a band matrix.
*           'Z' -- The whole matrix is packed as a band matrix.
*
            IF( IPACK.EQ.5 )
     $         UUB = 0
            IF( IPACK.EQ.6 )
     $         LLB = 0
*
            DO 380 J = 1, UUB
               DO 370 I = MIN( J+LLB, M ), 1, -1
                  A( I-J+UUB+1, J ) = A( I, J )
  370          CONTINUE
  380       CONTINUE
*
            DO 400 J = UUB + 2, N
               DO 390 I = J - UUB, MIN( J+LLB, M )
                  A( I-J+UUB+1, J ) = A( I, J )
  390          CONTINUE
  400       CONTINUE
         END IF
*
*        If packed, zero out extraneous elements.
*
*        Symmetric/Triangular Packed --
*        zero out everything after A(IROW,ICOL)
*
         IF( IPACK.EQ.3 .OR. IPACK.EQ.4 ) THEN
            DO 420 JC = ICOL, M
               DO 410 JR = IROW + 1, LDA
                  A( JR, JC ) = ZERO
  410          CONTINUE
               IROW = 0
  420       CONTINUE
*
         ELSE IF( IPACK.GE.5 ) THEN
*
*           Packed Band --
*              1st row is now in A( UUB+2-j, j), zero above it
*              m-th row is now in A( M+UUB-j,j), zero below it
*              last non-zero diagonal is now in A( UUB+LLB+1,j ),
*                 zero below it, too.
*
            IR1 = UUB + LLB + 2
            IR2 = UUB + M + 2
            DO 450 JC = 1, N
               DO 430 JR = 1, UUB + 1 - JC
                  A( JR, JC ) = ZERO
  430          CONTINUE
               DO 440 JR = MAX( 1, MIN( IR1, IR2-JC ) ), LDA
                  A( JR, JC ) = ZERO
  440          CONTINUE
  450       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DLATMS
*
      END
