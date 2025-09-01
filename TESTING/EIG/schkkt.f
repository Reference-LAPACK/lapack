*> \brief \b SCHKKT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SCHKKT( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
*                          NOUNIT, A, LDA, AP, SD, SE, D1, D2, D3, D4, D5,
*                          WA1, WA2, WA3, WR, U, LDU, V, VP, TAU, Z, WORK,
*                          LWORK, IWORK, LIWORK, RESULT, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES,
*      $                   NTYPES
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
*       REAL               A( LDA, * ), AP( * ), D1( * ), D2( * ),
*      $                   D3( * ), D4( * ), D5( * ), RESULT( * ),
*      $                   SD( * ), SE( * ), TAU( * ), U( LDU, * ),
*      $                   V( LDU, * ), VP( * ), WA1( * ), WA2( * ),
*      $                   WA3( * ), WORK( * ), WR( * ), Z( LDU, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SCHKKT  checks the skew-symmetric eigenvalue problem routines.
*>
*>    SKYTRD factors A as  U S U' , where ' means transpose,
*>    S is skew-symmetric tridiagonal, and U is orthogonal.
*>    SKYTRD can use either just the lower or just the upper triangle
*>    of A; SCHKKT checks both cases.
*>    U is represented as a product of Householder
*>    transformations, whose vectors are stored in the first
*>    n-1 columns of V, and whose scale factors are in TAU.
*>
*>    SKTEQR factors S as  Z D1 Z' , where Z is the orthogonal
*>    matrix of eigenvectors and D1 is a diagonal matrix with
*>    the eigenvalues on the diagonal.  D2 is the matrix of
*>    eigenvalues computed when Z is not computed.
*>
*> When SCHKKT is called, a number of matrix "sizes" ("n's") and a
*> number of matrix "types" are specified.  For each size ("n")
*> and each type of matrix, one matrix will be generated and used
*> to test the skew-symmetric eigenroutines.  For each matrix, a
*> number of tests will be performed:
*>
*> (1)     | A - V S V' | / ( |A| n ulp ) SKYTRD( UPLO='U', ... )
*>
*> (2)     | I - UV' | / ( n ulp )        SORGTR( UPLO='U', ... )
*>
*> (3)     | A - V S V' | / ( |A| n ulp ) SKYTRD( UPLO='L', ... )
*>
*> (4)     | I - UV' | / ( n ulp )        SORGTR( UPLO='L', ... )
*>
*> (5-8)   Same as 1-4, but for SSPTRD and SOPGTR.
*>
*> (9)     | S - Z D Z' | / ( |S| n ulp ) SKTEQR('V',...)
*>
*> (10)    | I - ZZ' | / ( n ulp )        SKTEQR('V',...)
*>
*> (11)    | D1 - D2 | / ( |D1| ulp )        SKTEQR('N',...)
*>
*> (12)    | D1 - D3 | / ( |D1| ulp )        SSTERF
*>
*> (13)    0 if the true eigenvalues (computed by sturm count)
*>         of S are within THRESH of
*>         those in D1.  2*THRESH if they are not.  (Tested using
*>         SSTECH)
*>
*> For S positive definite,
*>
*> (14)    | S - Z4 D4 Z4' | / ( |S| n ulp ) SPTEQR('V',...)
*>
*> (15)    | I - Z4 Z4' | / ( n ulp )        SPTEQR('V',...)
*>
*> (16)    | D4 - D5 | / ( 100 |D4| ulp )       SPTEQR('N',...)
*>
*> When S is also diagonally dominant by the factor gamma < 1,
*>
*> (17)    max | D4(i) - WR(i) | / ( |D4(i)| omega ) ,
*>          i
*>         omega = 2 (2n-1) ULP (1 + 8 gamma**2) / (1 - gamma)**4
*>                                              SSTEBZ( 'A', 'E', ...)
*>
*> (18)    | WA1 - D3 | / ( |D3| ulp )          SSTEBZ( 'A', 'E', ...)
*>
*> (19)    ( max { min | WA2(i)-WA3(j) | } +
*>            i     j
*>           max { min | WA3(i)-WA2(j) | } ) / ( |D3| ulp )
*>            i     j
*>                                              SSTEBZ( 'I', 'E', ...)
*>
*> (20)    | S - Y WA1 Y' | / ( |S| n ulp )  SSTEBZ, SSTEIN
*>
*> (21)    | I - Y Y' | / ( n ulp )          SSTEBZ, SSTEIN
*>
*> (22)    | S - Z D Z' | / ( |S| n ulp )    SSTEDC('I')
*>
*> (23)    | I - ZZ' | / ( n ulp )           SSTEDC('I')
*>
*> (24)    | S - Z D Z' | / ( |S| n ulp )    SSTEDC('V')
*>
*> (25)    | I - ZZ' | / ( n ulp )           SSTEDC('V')
*>
*> (26)    | D1 - D2 | / ( |D1| ulp )           SSTEDC('V') and
*>                                              SSTEDC('N')
*>
*> Test 27 is disabled at the moment because SSTEMR does not
*> guarantee high relatvie accuracy.
*>
*> (27)    max | D6(i) - WR(i) | / ( |D6(i)| omega ) ,
*>          i
*>         omega = 2 (2n-1) ULP (1 + 8 gamma**2) / (1 - gamma)**4
*>                                              SSTEMR('V', 'A')
*>
*> (28)    max | D6(i) - WR(i) | / ( |D6(i)| omega ) ,
*>          i
*>         omega = 2 (2n-1) ULP (1 + 8 gamma**2) / (1 - gamma)**4
*>                                              SSTEMR('V', 'I')
*>
*> Tests 29 through 34 are disable at present because SSTEMR
*> does not handle partial spectrum requests.
*>
*> (29)    | S - Z D Z' | / ( |S| n ulp )    SSTEMR('V', 'I')
*>
*> (30)    | I - ZZ' | / ( n ulp )           SSTEMR('V', 'I')
*>
*> (31)    ( max { min | WA2(i)-WA3(j) | } +
*>            i     j
*>           max { min | WA3(i)-WA2(j) | } ) / ( |D3| ulp )
*>            i     j
*>         SSTEMR('N', 'I') vs. SSTEMR('V', 'I')
*>
*> (32)    | S - Z D Z' | / ( |S| n ulp )    SSTEMR('V', 'V')
*>
*> (33)    | I - ZZ' | / ( n ulp )           SSTEMR('V', 'V')
*>
*> (34)    ( max { min | WA2(i)-WA3(j) | } +
*>            i     j
*>           max { min | WA3(i)-WA2(j) | } ) / ( |D3| ulp )
*>            i     j
*>         SSTEMR('N', 'V') vs. SSTEMR('V', 'V')
*>
*> (35)    | S - Z D Z' | / ( |S| n ulp )    SSTEMR('V', 'A')
*>
*> (36)    | I - ZZ' | / ( n ulp )           SSTEMR('V', 'A')
*>
*> (37)    ( max { min | WA2(i)-WA3(j) | } +
*>            i     j
*>           max { min | WA3(i)-WA2(j) | } ) / ( |D3| ulp )
*>            i     j
*>         SSTEMR('N', 'A') vs. SSTEMR('V', 'A')
*>
*> The "sizes" are specified by an array NN(1:NSIZES); the value of
*> each element NN(j) specifies one size.
*> The "types" are specified by a logical array DOTYPE( 1:NTYPES );
*> if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
*> Currently, the list of possible types is:
*>
*> (1)  The zero matrix.
*> (2)  The identity matrix.
*>
*> (3)  A diagonal matrix with evenly spaced entries
*>      1, ..., ULP  and random signs.
*>      (ULP = (first number larger than 1) - 1 )
*> (4)  A diagonal matrix with geometrically spaced entries
*>      1, ..., ULP  and random signs.
*> (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
*>      and random signs.
*>
*> (6)  Same as (4), but multiplied by SQRT( overflow threshold )
*> (7)  Same as (4), but multiplied by SQRT( underflow threshold )
*>
*> (8)  A matrix of the form  U' D U, where U is orthogonal and
*>      D has evenly spaced entries 1, ..., ULP with random signs
*>      on the diagonal.
*>
*> (9)  A matrix of the form  U' D U, where U is orthogonal and
*>      D has geometrically spaced entries 1, ..., ULP with random
*>      signs on the diagonal.
*>
*> (10) A matrix of the form  U' D U, where U is orthogonal and
*>      D has "clustered" entries 1, ULP,..., ULP with random
*>      signs on the diagonal.
*>
*> (11) Same as (8), but multiplied by SQRT( overflow threshold )
*> (12) Same as (8), but multiplied by SQRT( underflow threshold )
*>
*> (13) Symmetric matrix with random entries chosen from (-1,1).
*> (14) Same as (13), but multiplied by SQRT( overflow threshold )
*> (15) Same as (13), but multiplied by SQRT( underflow threshold )
*> (16) Same as (8), but diagonal elements are all positive.
*> (17) Same as (9), but diagonal elements are all positive.
*> (18) Same as (10), but diagonal elements are all positive.
*> (19) Same as (16), but multiplied by SQRT( overflow threshold )
*> (20) Same as (16), but multiplied by SQRT( underflow threshold )
*> (21) A diagonally dominant tridiagonal matrix with geometrically
*>      spaced diagonal entries 1, ..., ULP.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NSIZES
*> \verbatim
*>          NSIZES is INTEGER
*>          The number of sizes of matrices to use.  If it is zero,
*>          SCHKKT does nothing.  It must be at least zero.
*> \endverbatim
*>
*> \param[in] NN
*> \verbatim
*>          NN is INTEGER array, dimension (NSIZES)
*>          An array containing the sizes to be used for the matrices.
*>          Zero values will be skipped.  The values must be at least
*>          zero.
*> \endverbatim
*>
*> \param[in] NTYPES
*> \verbatim
*>          NTYPES is INTEGER
*>          The number of elements in DOTYPE.   If it is zero, SCHKKT
*>          does nothing.  It must be at least zero.  If it is MAXTYP+1
*>          and NSIZES is 1, then an additional type, MAXTYP+1 is
*>          defined, which is to use whatever matrix is in A.  This
*>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
*>          DOTYPE(MAXTYP+1) is .TRUE. .
*> \endverbatim
*>
*> \param[in] DOTYPE
*> \verbatim
*>          DOTYPE is LOGICAL array, dimension (NTYPES)
*>          If DOTYPE(j) is .TRUE., then for each size in NN a
*>          matrix of that size and of type j will be generated.
*>          If NTYPES is smaller than the maximum number of types
*>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
*>          MAXTYP will not be generated.  If NTYPES is larger
*>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
*>          will be ignored.
*> \endverbatim
*>
*> \param[in,out] ISEED
*> \verbatim
*>          ISEED is INTEGER array, dimension (4)
*>          On entry ISEED specifies the seed of the random number
*>          generator. The array elements should be between 0 and 4095;
*>          if not they will be reduced mod 4096.  Also, ISEED(4) must
*>          be odd.  The random number generator uses a linear
*>          congruential sequence limited to small integers, and so
*>          should produce machine independent random numbers. The
*>          values of ISEED are changed on exit, and can be used in the
*>          next call to SCHKKT to continue the same random number
*>          sequence.
*> \endverbatim
*>
*> \param[in] THRESH
*> \verbatim
*>          THRESH is REAL
*>          A test will count as "failed" if the "error", computed as
*>          described above, exceeds THRESH.  Note that the error
*>          is scaled to be O(1), so THRESH should be a reasonably
*>          small multiple of 1, e.g., 10 or 100.  In particular,
*>          it should not depend on the precision (single vs. double)
*>          or the size of the matrix.  It must be at least zero.
*> \endverbatim
*>
*> \param[in] NOUNIT
*> \verbatim
*>          NOUNIT is INTEGER
*>          The FORTRAN unit number for printing out error messages
*>          (e.g., if a routine returns IINFO not equal to 0.)
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array of
*>                                  dimension ( LDA , max(NN) )
*>          Used to hold the matrix whose eigenvalues are to be
*>          computed.  On exit, A contains the last matrix actually
*>          used.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  It must be at
*>          least 1 and at least max( NN ).
*> \endverbatim
*>
*> \param[out] AP
*> \verbatim
*>          AP is REAL array of
*>                      dimension( max(NN)*max(NN+1)/2 )
*>          The matrix A stored in packed format.
*> \endverbatim
*>
*> \param[out] SD
*> \verbatim
*>          SD is REAL array of
*>                             dimension( max(NN) )
*>          The diagonal of the tridiagonal matrix computed by SKYTRD.
*>          On exit, SD and SE contain the tridiagonal form of the
*>          matrix in A.
*> \endverbatim
*>
*> \param[out] SE
*> \verbatim
*>          SE is REAL array of
*>                             dimension( max(NN) )
*>          The off-diagonal of the tridiagonal matrix computed by
*>          SKYTRD.  On exit, SD and SE contain the tridiagonal form of
*>          the matrix in A.
*> \endverbatim
*>
*> \param[out] D1
*> \verbatim
*>          D1 is REAL array of
*>                             dimension( max(NN) )
*>          The eigenvalues of A, as computed by SKTEQR simlutaneously
*>          with Z.  On exit, the eigenvalues in D1 correspond with the
*>          matrix in A.
*> \endverbatim
*>
*> \param[out] D2
*> \verbatim
*>          D2 is REAL array of
*>                             dimension( max(NN) )
*>          The eigenvalues of A, as computed by SKTEQR if Z is not
*>          computed.  On exit, the eigenvalues in D2 correspond with
*>          the matrix in A.
*> \endverbatim
*>
*> \param[out] D3
*> \verbatim
*>          D3 is REAL array of
*>                             dimension( max(NN) )
*>          The eigenvalues of A, as computed by SSTERF.  On exit, the
*>          eigenvalues in D3 correspond with the matrix in A.
*> \endverbatim
*>
*> \param[out] D4
*> \verbatim
*>          D4 is REAL array of
*>                             dimension( max(NN) )
*>          The eigenvalues of A, as computed by SPTEQR(V).
*>          ZPTEQR factors S as  Z4 D4 Z4*
*>          On exit, the eigenvalues in D4 correspond with the matrix in A.
*> \endverbatim
*>
*> \param[out] D5
*> \verbatim
*>          D5 is REAL array of
*>                             dimension( max(NN) )
*>          The eigenvalues of A, as computed by SPTEQR(N)
*>          when Z is not computed. On exit, the
*>          eigenvalues in D4 correspond with the matrix in A.
*> \endverbatim
*>
*> \param[out] WA1
*> \verbatim
*>          WA1 is REAL array of
*>                             dimension( max(NN) )
*>          All eigenvalues of A, computed to high
*>          absolute accuracy, with different range options.
*>          as computed by SSTEBZ.
*> \endverbatim
*>
*> \param[out] WA2
*> \verbatim
*>          WA2 is REAL array of
*>                             dimension( max(NN) )
*>          Selected eigenvalues of A, computed to high
*>          absolute accuracy, with different range options.
*>          as computed by SSTEBZ.
*>          Choose random values for IL and IU, and ask for the
*>          IL-th through IU-th eigenvalues.
*> \endverbatim
*>
*> \param[out] WA3
*> \verbatim
*>          WA3 is REAL array of
*>                             dimension( max(NN) )
*>          Selected eigenvalues of A, computed to high
*>          absolute accuracy, with different range options.
*>          as computed by SSTEBZ.
*>          Determine the values VL and VU of the IL-th and IU-th
*>          eigenvalues and ask for all eigenvalues in this range.
*> \endverbatim
*>
*> \param[out] WR
*> \verbatim
*>          WR is REAL array of
*>                             dimension( max(NN) )
*>          All eigenvalues of A, computed to high
*>          absolute accuracy, with different options.
*>          as computed by SSTEBZ.
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is REAL array of
*>                             dimension( LDU, max(NN) ).
*>          The orthogonal matrix computed by SKYTRD + SORGTR.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of U, Z, and V.  It must be at least 1
*>          and at least max( NN ).
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is REAL array of
*>                             dimension( LDU, max(NN) ).
*>          The Housholder vectors computed by SKYTRD in reducing A to
*>          tridiagonal form.  The vectors computed with UPLO='U' are
*>          in the upper triangle, and the vectors computed with UPLO='L'
*>          are in the lower triangle.  (As described in SKYTRD, the
*>          sub- and superdiagonal are not set to 1, although the
*>          true Householder vector has a 1 in that position.  The
*>          routines that use V, such as SORGTR, set those entries to
*>          1 before using them, and then restore them later.)
*> \endverbatim
*>
*> \param[out] VP
*> \verbatim
*>          VP is REAL array of
*>                      dimension( max(NN)*max(NN+1)/2 )
*>          The matrix V stored in packed format.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is REAL array of
*>                             dimension( max(NN) )
*>          The Householder factors computed by SKYTRD in reducing A
*>          to tridiagonal form.
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is REAL array of
*>                             dimension( LDU, max(NN) ).
*>          The orthogonal matrix of eigenvectors computed by SKTEQR,
*>          SPTEQR, and SSTEIN.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array of
*>                      dimension( LWORK )
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The number of entries in WORK.  This must be at least
*>          1 + 4 * Nmax + 2 * Nmax * lg Nmax + 3 * Nmax**2
*>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array,
*>          Workspace.
*> \endverbatim
*>
*> \param[out] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          The number of entries in IWORK.  This must be at least
*>                  6 + 6*Nmax + 5 * Nmax * lg Nmax
*>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is REAL array, dimension (26)
*>          The values computed by the tests described above.
*>          The values are currently limited to 1/ulp, to avoid
*>          overflow.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          If 0, then everything ran OK.
*>           -1: NSIZES < 0
*>           -2: Some NN(j) < 0
*>           -3: NTYPES < 0
*>           -5: THRESH < 0
*>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
*>          -23: LDU < 1 or LDU < NMAX.
*>          -29: LWORK too small.
*>          If  SLATMR, SLATMS, SKYTRD, SORGTR, SKTEQR, SSTERF,
*>              or SORMC2 returns an error code, the
*>              absolute value of it is returned.
*>
*>-----------------------------------------------------------------------
*>
*>       Some Local Variables and Parameters:
*>       ---- ----- --------- --- ----------
*>       ZERO, ONE       Real 0 and 1.
*>       MAXTYP          The number of types defined.
*>       NTEST           The number of tests performed, or which can
*>                       be performed so far, for the current matrix.
*>       NTESTT          The total number of tests performed so far.
*>       NBLOCK          Blocksize as returned by ENVIR.
*>       NMAX            Largest value in NN.
*>       NMATS           The number of matrices generated so far.
*>       NERRS           The number of tests which have exceeded THRESH
*>                       so far.
*>       COND, IMODE     Values to be passed to the matrix generators.
*>       ANORM           Norm of A; passed to matrix generators.
*>
*>       OVFL, UNFL      Overflow and underflow thresholds.
*>       ULP, ULPINV     Finest relative precision and its inverse.
*>       RTOVFL, RTUNFL  Square roots of the previous 2 values.
*>               The following four arrays decode JTYPE:
*>       KTYPE(j)        The general type (1-10) for type "j".
*>       KMODE(j)        The MODE value to be passed to the matrix
*>                       generator for type "j".
*>       KMAGN(j)        The order of magnitude ( O(1),
*>                       O(overflow^(1/2) ), O(underflow^(1/2) )
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
*> \ingroup single_eig
*
*  =====================================================================
      SUBROUTINE SCHKKT( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $                   NOUNIT, A, LDA, AP, SD, SE, D1, D2, D3, D4, D5,
     $                   WA1, WA2, WA3, WR, U, LDU, V, VP, TAU, Z, WORK,
     $                   LWORK, IWORK, LIWORK, RESULT, INFO )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES,
     $                   NTYPES
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
      REAL               A( LDA, * ), AP( * ), D1( * ), D2( * ),
     $                   D3( * ), D4( * ), D5( * ), RESULT( * ),
     $                   SD( * ), SE( * ), TAU( * ), U( LDU, * ),
     $                   V( LDU, * ), VP( * ), WA1( * ), WA2( * ),
     $                   WA3( * ), WORK( * ), WR( * ), Z( LDU, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO, EIGHT, TEN, HUN
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   EIGHT = 8.0E0, TEN = 10.0E0, HUN = 100.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = ONE / TWO )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 21 )
      LOGICAL            SRANGE
      PARAMETER          ( SRANGE = .FALSE. )
      LOGICAL            SREL
      PARAMETER          ( SREL = .FALSE. )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADNN, TRYRAC
      INTEGER            I, IINFO, IL, IMODE, ITEMP, ITYPE, IU, J, JC,
     $                   JR, JSIZE, JTYPE, LGN, LIWEDC, LOG2UI, LWEDC,
     $                   M, M2, M3, MTYPES, N, NAP, NBLOCK, NERRS,
     $                   NMATS, NMAX, NSPLIT, NTEST, NTESTT
      REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL,
     $                   RTUNFL, TEMP1, TEMP2, TEMP3, TEMP4, ULP,
     $                   ULPINV, UNFL, VL, VU
*     ..
*     .. Local Arrays ..
      INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ),
     $                   KMAGN( MAXTYP ), KMODE( MAXTYP ),
     $                   KTYPE( MAXTYP )
      REAL               DUMMA( 1 )
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      REAL               SLAMCH
      EXTERNAL           ILAENV, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLABAD, SLACPY, SLASET, SLASUM, SLATMR,
     $                   SLATMS, SORGTR, SKTEQR, SKTT21, SKYT21,
     $                   SKYTRD, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, REAL, SQRT
*     ..
*     .. Data statements ..
      DATA               KTYPE / 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8,
     $                   8, 8, 9, 9, 9, 9, 9, 10 /
      DATA               KMAGN / 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1,
     $                   2, 3, 1, 1, 1, 2, 3, 1 /
      DATA               KMODE / 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0,
     $                   0, 0, 4, 3, 1, 4, 4, 3 /
*     ..
*     .. Executable Statements ..
*
*     Keep ftnchek happy
      IDUMMA( 1 ) = 1
*
*     Check for errors
*
      NTESTT = 0
      INFO = 0
*
*     Important constants
*
      BADNN = .FALSE.
      TRYRAC = .TRUE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 )
     $      BADNN = .TRUE.
   10 CONTINUE
*
      NBLOCK = ILAENV( 1, 'SKYTRD', 'L', NMAX, -1, -1, -1 )
      NBLOCK = MIN( NMAX, MAX( 1, NBLOCK ) )
*
*     Check for errors
*
      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADNN ) THEN
         INFO = -2
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.NMAX ) THEN
         INFO = -9
      ELSE IF( LDU.LT.NMAX ) THEN
         INFO = -23
      ELSE IF( 2*MAX( 2, NMAX )**2.GT.LWORK ) THEN
         INFO = -29
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SCHKKT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 )
     $   RETURN
*
*     More Important constants
*
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL SLABAD( UNFL, OVFL )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      LOG2UI = INT( LOG( ULPINV ) / LOG( TWO ) )
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
*
*     Loop over sizes, types
*
      DO 20 I = 1, 4
         ISEED2( I ) = ISEED( I )
   20 CONTINUE
      NERRS = 0
      NMATS = 0
*
      DO 310 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         IF( N.GT.0 ) THEN
            LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
            IF( 2**LGN.LT.N )
     $         LGN = LGN + 1
            IF( 2**LGN.LT.N )
     $         LGN = LGN + 1
            LWEDC = 1 + 4*N + 2*N*LGN + 4*N**2
            LIWEDC = 6 + 6*N + 5*N*LGN
         ELSE
            LWEDC = 8
            LIWEDC = 12
         END IF
         NAP = ( N*( N+1 ) ) / 2
         ANINV = ONE / REAL( MAX( 1, N ) )
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 300 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) )
     $         GO TO 300
            NMATS = NMATS + 1
            NTEST = 0
*
            DO 30 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   30       CONTINUE
*
*           Compute "A"
*
*           Control parameters:
*
*               KMAGN  KMODE        KTYPE
*           =1  O(1)   clustered 1  zero
*           =2  large  clustered 2  identity
*           =3  small  exponential  (none)
*           =4         arithmetic   diagonal, (w/ eigenvalues)
*           =5         random log   symmetric, w/ eigenvalues
*           =6         random       (none)
*           =7                      random diagonal
*           =8                      random symmetric
*           =9                      positive definite
*           =10                     diagonally dominant tridiagonal
*
            IF( MTYPES.GT.MAXTYP )
     $         GO TO 100
*
            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )
*
*           Compute norm
*
            GO TO ( 40, 50, 60 )KMAGN( JTYPE )
*
   40       CONTINUE
            ANORM = ONE
            GO TO 70
*
   50       CONTINUE
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 70
*
   60       CONTINUE
            ANORM = RTUNFL*N*ULPINV
            GO TO 70
*
   70       CONTINUE
*
            CALL SLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
            IINFO = 0
            IF( JTYPE.LE.15 ) THEN
               COND = ULPINV
            ELSE
               COND = ULPINV*ANINV / TEN
            END IF
*
*           Special Matrices -- Identity & Jordan block
*
*              Zero
*
            IF( ITYPE.EQ.1 ) THEN
               IINFO = 0
*
            ELSE IF( ITYPE.EQ.2 ) THEN
*
*              Identity
*
               DO 80 JC = 1, N
                  A( JC, JC ) = ANORM
   80          CONTINUE
*
            ELSE IF( ITYPE.EQ.4 ) THEN
*
*              tridiagonal Matrix, [Eigen]values Specified
*
               CALL SLATMS( N, N, 'S', ISEED, 'K', WORK, IMODE, COND,
     $                      ANORM, 1, 1, 'N', A, LDA, WORK( N+1 ),
     $                      IINFO )
*
*
            ELSE IF( ITYPE.EQ.5 ) THEN
*
*              skew-ymmetric, eigenvalues specified
*
               CALL SLATMS( N, N, 'S', ISEED, 'K', WORK, IMODE, COND,
     $                      ANORM, N, N, 'N', A, LDA, WORK( N+1 ),
     $                      IINFO )
*
            ELSE IF( ITYPE.EQ.7 ) THEN
*
*              tridiagonal, random eigenvalues
*
               CALL SLATMR( N, N, 'S', ISEED, 'K', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 1, 1,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.8 ) THEN
*
*              skew-ymmetric, random eigenvalues
*
               CALL SLATMR( N, N, 'S', ISEED, 'K', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.9 ) THEN
*
*              skew-ymmetric, eigenvalues specified.
*
               CALL SLATMS( N, N, 'S', ISEED, 'K', WORK, IMODE, COND,
     $                      ANORM, N, N, 'N', A, LDA, WORK( N+1 ),
     $                      IINFO )
*
            ELSE IF( ITYPE.EQ.10 ) THEN
*
*              skew-ymmetric tridiagonal, eigenvalues specified.
*
               CALL SLATMS( N, N, 'S', ISEED, 'K', WORK, IMODE, COND,
     $                      ANORM, 1, 1, 'N', A, LDA, WORK( N+1 ),
     $                      IINFO )
               DO 90 I = 2, N
                  TEMP1 = ABS( A( I-1, I ) ) /
     $                    SQRT( ABS( A( I-1, I-1 )*A( I, I ) ) )
                  IF( TEMP1.GT.HALF ) THEN
                     A( I-1, I ) = HALF*SQRT( ABS( A( I-1, I-1 )*A( I,
     $                             I ) ) )
                     A( I, I-1 ) = A( I-1, I )
                  END IF
   90          CONTINUE
*
            ELSE
*
               IINFO = 1
            END IF
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
*
  100       CONTINUE
*
*           Call SKYTRD and SORGTR to compute S and U from
*           upper triangle.
*
            CALL SLACPY( 'U', N, N, A, LDA, V, LDU )
*
            NTEST = 1
            CALL SKYTRD( 'U', N, V, LDU, SE, TAU, WORK, LWORK,
     $                   IINFO )
            CALL SLASET( 'N', N, 1, ZERO, ZERO, SD, N)
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SKYTRD(U)', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 1 ) = ULPINV
                  GO TO 280
               END IF
            END IF
*
            CALL SLACPY( 'U', N, N, V, LDU, U, LDU )
*
            NTEST = 2
            CALL SORGTR( 'U', N, U, LDU, TAU, WORK, LWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SORGTR(U)', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 2 ) = ULPINV
                  GO TO 280
               END IF
            END IF
*
*           Do tests 1 and 2
*
            CALL SKYT21( 2, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V,
     $                   LDU, TAU, WORK, RESULT( 1 ) )
            CALL SKYT21( 3, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V,
     $                   LDU, TAU, WORK, RESULT( 2 ) )
*
*           Call SKYTRD and SORGTR to compute S and U from
*           lower triangle, do tests.
*
            CALL SLACPY( 'L', N, N, A, LDA, V, LDU )
*
            NTEST = 3
            CALL SKYTRD( 'L', N, V, LDU, SE, TAU, WORK, LWORK,
     $                   IINFO )
            CALL SLASET( 'N', N, 1, ZERO, ZERO, SD, N)
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SKYTRD(L)', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 3 ) = ULPINV
                  GO TO 280
               END IF
            END IF
*
            CALL SLACPY( 'L', N, N, V, LDU, U, LDU )
*
            NTEST = 4
            CALL SORGTR( 'L', N, U, LDU, TAU, WORK, LWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SORGTR(L)', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 4 ) = ULPINV
                  GO TO 280
               END IF
            END IF
*
*           Do tests 3 and 4
*
            CALL SKYT21( 2, 'Lower', N, 1, A, LDA, SD, SE, U, LDU, V,
     $                   LDU, TAU, WORK, RESULT( 3 ) )
            CALL SKYT21( 3, 'Lower', N, 1, A, LDA, SD, SE, U, LDU, V,
     $                   LDU, TAU, WORK, RESULT( 4 ) )
*
*           Call SKTEQR to compute D1, D2, and Z, do tests.
*
*           Compute D1 and Z
*
            CALL SCOPY( N, SD, 1, D1, 1 )
            IF( N.GT.0 )
     $         CALL SCOPY( N-1, SE, 1, WORK, 1 )
            CALL SLASET( 'Full', N, N, ZERO, ONE, Z, LDU )
*
            NTEST = 5
            CALL SKTEQR( 'V', N, WORK, Z, LDU, WORK( N+1 ), IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SKTEQR(V)', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 5 ) = ULPINV
                  GO TO 280
               END IF
            END IF
            IF( N.GT.0 )
     $         CALL SCOPY( N-1, WORK, 1, D1, 1 )
*
*           Compute D2
*
            CALL SCOPY( N, SD, 1, D2, 1 )
            IF( N.GT.0 )
     $         CALL SCOPY( N-1, SE, 1, WORK, 1 )
*
            NTEST = 7
            CALL SKTEQR( 'N', N, WORK, WORK( N+1 ), LDU,
     $                   WORK( N+1 ), IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SKTEQR(N)', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 6 ) = ULPINV
                  GO TO 280
               END IF
            END IF
            IF( N.GT.0 )
     $         CALL SCOPY( N-1, WORK, 1, D2, 1 )
*
*           Do Tests 5 and 6
*
            CALL SKTT21( N, 1, DUMMA, SE, DUMMA, D1, Z, LDU, WORK,
     $                   RESULT( 5 ) )
*
*           Do Tests 7
*
            TEMP1 = ZERO
            TEMP2 = ZERO
*
            DO 150 J = 1, N-1
               TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
  150       CONTINUE
*
            RESULT( 7 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
  280       CONTINUE
            NTESTT = NTESTT + NTEST
*
*           End of Loop -- Check for RESULT(j) > THRESH
*
*
*           Print out tests which fail.
*
            DO 290 JR = 1, NTEST
               IF( RESULT( JR ).GE.THRESH ) THEN
*
*                 If this is the first test to fail,
*                 print a header to the data file.
*
                  IF( NERRS.EQ.0 ) THEN
                     WRITE( NOUNIT, FMT = 9998 )'SKT'
                     WRITE( NOUNIT, FMT = 9997 )
                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )'Skew-symmetric'
                     WRITE( NOUNIT, FMT = 9994 )
*
*                    Tests performed
*
                     WRITE( NOUNIT, FMT = 9988 )
                  END IF
                  NERRS = NERRS + 1
                  WRITE( NOUNIT, FMT = 9990 )N, IOLDSD, JTYPE, JR,
     $               RESULT( JR )
               END IF
  290       CONTINUE
  300    CONTINUE
  310 CONTINUE
*
*     Summary
*
      CALL SLASUM( 'SKT', NOUNIT, NERRS, NTESTT )
      RETURN
*
 9999 FORMAT( ' SCHKKT: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
 9998 FORMAT( / 1X, A3, ' -- Real Skew-symmetric eigenvalue problem' )
 9997 FORMAT( ' Matrix types (see SCHKKT for details): ' )
*
 9996 FORMAT( / ' Special Matrices:',
     $      / '  1=Zero matrix.                        ',
     $      '  5=Diagonal: clustered entries.',
     $      / '  2=Identity matrix.                    ',
     $      '  6=Diagonal: large, evenly spaced.',
     $      / '  3=Diagonal: evenly spaced entries.    ',
     $      '  7=Diagonal: small, evenly spaced.',
     $      / '  4=Diagonal: geometr. spaced entries.' )
 9995 FORMAT( ' Dense ', A, ' Matrices:',
     $      / '  8=Evenly spaced eigenvals.            ',
     $      ' 12=Small, evenly spaced eigenvals.',
     $      / '  9=Geometrically spaced eigenvals.     ',
     $      ' 13=Matrix with random O(1) entries.',
     $      / ' 10=Clustered eigenvalues.              ',
     $      ' 14=Matrix with large random entries.',
     $      / ' 11=Large, evenly spaced eigenvals.     ',
     $      ' 15=Matrix with small random entries.' )
 9994 FORMAT( ' 16=Positive definite, evenly spaced eigenvalues',
     $      / ' 17=Positive definite, geometrically spaced eigenvlaues',
     $      / ' 18=Positive definite, clustered eigenvalues',
     $      / ' 19=Positive definite, small evenly spaced eigenvalues',
     $      / ' 20=Positive definite, large evenly spaced eigenvalues',
     $      / ' 21=Diagonally dominant tridiagonal, geometrically',
     $      ' spaced eigenvalues' )
*
 9990 FORMAT( ' N=', I5, ', seed=', 4( I4, ',' ), ' type ', I2,
     $      ', test(', I2, ')=', G10.3 )
*
 9988 FORMAT( / 'Test performed:  see SCHKKT for details.', / )
*     End of SCHKKT
*
      END
