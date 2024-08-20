*> \brief \b SDRVKT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SDRVKT( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
*                          NOUNIT, A, LDA, D1, D2, D3, D4, EVEIGS, WA1,
*                          WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK,
*                          IWORK, LIWORK, RESULT, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES,
*      $                   NTYPES
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
*       REAL               A( LDA, * ), D1( * ), D2( * ), D3( * ),
*      $                   D4( * ), EVEIGS( * ), RESULT( * ), TAU( * ),
*      $                   U( LDU, * ), V( LDU, * ), WA1( * ), WA2( * ),
*      $                   WA3( * ), WORK( * ), Z( LDU, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>      SDRVKT  checks the skew-symmetric eigenvalue problem drivers.
*>
*>              SKTEV computes all eigenvalues and, optionally,
*>              eigenvectors of a real skew-symmetric tridiagonal matrix.
*>
*>              SKYEV computes all eigenvalues and, optionally,
*>              eigenvectors of a real skew-symmetric matrix.
*>
*>      When SDRVKT is called, a number of matrix "sizes" ("n's") and a
*>      number of matrix "types" are specified.  For each size ("n")
*>      and each type of matrix, one matrix will be generated and used
*>      to test the appropriate drivers.  For each matrix and each
*>      driver routine called, the following tests will be performed:
*>
*>      (1)     | A - Z D Z' | / ( |A| n ulp )
*>
*>      (2)     | I - Z Z' | / ( n ulp )
*>
*>      (3)     | D1 - D2 | / ( |D1| ulp )
*>
*>      where Z is the matrix of eigenvectors returned when the
*>      eigenvector option is given and D1 and D2 are the eigenvalues
*>      returned with and without the eigenvector option.
*>
*>      The "sizes" are specified by an array NN(1:NSIZES); the value of
*>      each element NN(j) specifies one size.
*>      The "types" are specified by a logical array DOTYPE( 1:NTYPES );
*>      if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
*>      Currently, the list of possible types is:
*>
*>      (1)  The zero matrix.
*>      (2)  The identity matrix.
*>
*>      (3)  A diagonal matrix with evenly spaced eigenvalues
*>           1, ..., ULP  and random signs.
*>           (ULP = (first number larger than 1) - 1 )
*>      (4)  A diagonal matrix with geometrically spaced eigenvalues
*>           1, ..., ULP  and random signs.
*>      (5)  A diagonal matrix with "clustered" eigenvalues
*>           1, ULP, ..., ULP and random signs.
*>
*>      (6)  Same as (4), but multiplied by SQRT( overflow threshold )
*>      (7)  Same as (4), but multiplied by SQRT( underflow threshold )
*>
*>      (8)  A matrix of the form  U' D U, where U is orthogonal and
*>           D has evenly spaced entries 1, ..., ULP with random signs
*>           on the diagonal.
*>
*>      (9)  A matrix of the form  U' D U, where U is orthogonal and
*>           D has geometrically spaced entries 1, ..., ULP with random
*>           signs on the diagonal.
*>
*>      (10) A matrix of the form  U' D U, where U is orthogonal and
*>           D has "clustered" entries 1, ULP,..., ULP with random
*>           signs on the diagonal.
*>
*>      (11) Same as (8), but multiplied by SQRT( overflow threshold )
*>      (12) Same as (8), but multiplied by SQRT( underflow threshold )
*>
*>      (13) skew-symmetric matrix with random entries chosen from (-1,1).
*>      (14) Same as (13), but multiplied by SQRT( overflow threshold )
*>      (15) Same as (13), but multiplied by SQRT( underflow threshold )
*>      (16) A band matrix with half bandwidth randomly chosen between
*>           0 and N-1, with evenly spaced eigenvalues 1, ..., ULP
*>           with random signs.
*>      (17) Same as (16), but multiplied by SQRT( overflow threshold )
*>      (18) Same as (16), but multiplied by SQRT( underflow threshold )
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \verbatim
*>  NSIZES  INTEGER
*>          The number of sizes of matrices to use.  If it is zero,
*>          SDRVKT does nothing.  It must be at least zero.
*>          Not modified.
*>
*>  NN      INTEGER array, dimension (NSIZES)
*>          An array containing the sizes to be used for the matrices.
*>          Zero values will be skipped.  The values must be at least
*>          zero.
*>          Not modified.
*>
*>  NTYPES  INTEGER
*>          The number of elements in DOTYPE.   If it is zero, SDRVKT
*>          does nothing.  It must be at least zero.  If it is MAXTYP+1
*>          and NSIZES is 1, then an additional type, MAXTYP+1 is
*>          defined, which is to use whatever matrix is in A.  This
*>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
*>          DOTYPE(MAXTYP+1) is .TRUE. .
*>          Not modified.
*>
*>  DOTYPE  LOGICAL array, dimension (NTYPES)
*>          If DOTYPE(j) is .TRUE., then for each size in NN a
*>          matrix of that size and of type j will be generated.
*>          If NTYPES is smaller than the maximum number of types
*>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
*>          MAXTYP will not be generated.  If NTYPES is larger
*>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
*>          will be ignored.
*>          Not modified.
*>
*>  ISEED   INTEGER array, dimension (4)
*>          On entry ISEED specifies the seed of the random number
*>          generator. The array elements should be between 0 and 4095;
*>          if not they will be reduced mod 4096.  Also, ISEED(4) must
*>          be odd.  The random number generator uses a linear
*>          congruential sequence limited to small integers, and so
*>          should produce machine independent random numbers. The
*>          values of ISEED are changed on exit, and can be used in the
*>          next call to SDRVKT to continue the same random number
*>          sequence.
*>          Modified.
*>
*>  THRESH  REAL
*>          A test will count as "failed" if the "error", computed as
*>          described above, exceeds THRESH.  Note that the error
*>          is scaled to be O(1), so THRESH should be a reasonably
*>          small multiple of 1, e.g., 10 or 100.  In particular,
*>          it should not depend on the precision (single vs. double)
*>          or the size of the matrix.  It must be at least zero.
*>          Not modified.
*>
*>  NOUNIT  INTEGER
*>          The FORTRAN unit number for printing out error messages
*>          (e.g., if a routine returns IINFO not equal to 0.)
*>          Not modified.
*>
*>  A       REAL array, dimension (LDA , max(NN))
*>          Used to hold the matrix whose eigenvalues are to be
*>          computed.  On exit, A contains the last matrix actually
*>          used.
*>          Modified.
*>
*>  LDA     INTEGER
*>          The leading dimension of A.  It must be at
*>          least 1 and at least max( NN ).
*>          Not modified.
*>
*>  D1      REAL array, dimension (max(NN))
*>          The eigenvalues of A, as computed by SSTEQR simlutaneously
*>          with Z.  On exit, the eigenvalues in D1 correspond with the
*>          matrix in A.
*>          Modified.
*>
*>  D2      REAL array, dimension (max(NN))
*>          The eigenvalues of A, as computed by SSTEQR if Z is not
*>          computed.  On exit, the eigenvalues in D2 correspond with
*>          the matrix in A.
*>          Modified.
*>
*>  D3      REAL array, dimension (max(NN))
*>          The eigenvalues of A, as computed by SSTERF.  On exit, the
*>          eigenvalues in D3 correspond with the matrix in A.
*>          Modified.
*>
*>  D4      REAL array, dimension
*>
*>  EVEIGS  REAL array, dimension (max(NN))
*>          The eigenvalues as computed by SKTEV('N', ... )
*>          (I reserve the right to change this to the output of
*>          whichever algorithm computes the most accurate eigenvalues).
*>
*>  WA1     REAL array, dimension
*>
*>  WA2     REAL array, dimension
*>
*>  WA3     REAL array, dimension
*>
*>  U       REAL array, dimension (LDU, max(NN))
*>          The orthogonal matrix computed by SSYTRD + SORGTR.
*>          Modified.
*>
*>  LDU     INTEGER
*>          The leading dimension of U, Z, and V.  It must be at
*>          least 1 and at least max( NN ).
*>          Not modified.
*>
*>  V       REAL array, dimension (LDU, max(NN))
*>          The Housholder vectors computed by SSYTRD in reducing A to
*>          tridiagonal form.
*>          Modified.
*>
*>  TAU     REAL array, dimension (max(NN))
*>          The Householder factors computed by SSYTRD in reducing A
*>          to tridiagonal form.
*>          Modified.
*>
*>  Z       REAL array, dimension (LDU, max(NN))
*>          The orthogonal matrix of eigenvectors computed by SSTEQR,
*>          SPTEQR, and SSTEIN.
*>          Modified.
*>
*>  WORK    REAL array, dimension (LWORK)
*>          Workspace.
*>          Modified.
*>
*>  LWORK   INTEGER
*>          The number of entries in WORK.  This must be at least
*>          1 + 4 * Nmax + 2 * Nmax * lg Nmax + 4 * Nmax**2
*>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
*>          Not modified.
*>
*>  IWORK   INTEGER array,
*>             dimension (6 + 6*Nmax + 5 * Nmax * lg Nmax )
*>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
*>          Workspace.
*>          Modified.
*>
*>  RESULT  REAL array, dimension (105)
*>          The values computed by the tests described above.
*>          The values are currently limited to 1/ulp, to avoid
*>          overflow.
*>          Modified.
*>
*>  INFO    INTEGER
*>          If 0, then everything ran OK.
*>           -1: NSIZES < 0
*>           -2: Some NN(j) < 0
*>           -3: NTYPES < 0
*>           -5: THRESH < 0
*>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
*>          -16: LDU < 1 or LDU < NMAX.
*>          -21: LWORK too small.
*>          If  SLATMR, SLATMS, SSYTRD, SORGTR, SSTEQR, SSTERF,
*>              or SORMTR returns an error code, the
*>              absolute value of it is returned.
*>          Modified.
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
*>       NMAX            Largest value in NN.
*>       NMATS           The number of matrices generated so far.
*>       NERRS           The number of tests which have exceeded THRESH
*>                       so far (computed by SLAFTS).
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
*>
*>     The tests performed are:                 Routine tested
*>    1= | A - U S U' | / ( |A| n ulp )         SKTEV('V', ... )
*>    2= | I - U U' | / ( n ulp )               SKTEV('V', ... )
*>    3= |D(with Z) - D(w/o Z)| / (|D| ulp)     SKTEV('N', ... )
*>    4= | A - U S U' | / ( |A| n ulp )         SSTEVX('V','A', ... )
*>    5= | I - U U' | / ( n ulp )               SSTEVX('V','A', ... )
*>    6= |D(with Z) - EVEIGS| / (|D| ulp)       SSTEVX('N','A', ... )
*>    7= | A - U S U' | / ( |A| n ulp )         SSTEVR('V','A', ... )
*>    8= | I - U U' | / ( n ulp )               SSTEVR('V','A', ... )
*>    9= |D(with Z) - EVEIGS| / (|D| ulp)       SSTEVR('N','A', ... )
*>    10= | A - U S U' | / ( |A| n ulp )        SSTEVX('V','I', ... )
*>    11= | I - U U' | / ( n ulp )              SSTEVX('V','I', ... )
*>    12= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSTEVX('N','I', ... )
*>    13= | A - U S U' | / ( |A| n ulp )        SSTEVX('V','V', ... )
*>    14= | I - U U' | / ( n ulp )              SSTEVX('V','V', ... )
*>    15= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSTEVX('N','V', ... )
*>    16= | A - U S U' | / ( |A| n ulp )        SSTEVD('V', ... )
*>    17= | I - U U' | / ( n ulp )              SSTEVD('V', ... )
*>    18= |D(with Z) - EVEIGS| / (|D| ulp)      SSTEVD('N', ... )
*>    19= | A - U S U' | / ( |A| n ulp )        SSTEVR('V','I', ... )
*>    20= | I - U U' | / ( n ulp )              SSTEVR('V','I', ... )
*>    21= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSTEVR('N','I', ... )
*>    22= | A - U S U' | / ( |A| n ulp )        SSTEVR('V','V', ... )
*>    23= | I - U U' | / ( n ulp )              SSTEVR('V','V', ... )
*>    24= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSTEVR('N','V', ... )
*>
*>    25= | A - U S U' | / ( |A| n ulp )        SKYEV('L','V', ... )
*>    26= | I - U U' | / ( n ulp )              SKYEV('L','V', ... )
*>    27= |D(with Z) - D(w/o Z)| / (|D| ulp)    SKYEV('L','N', ... )
*>    28= | A - U S U' | / ( |A| n ulp )        SSYEVX('L','V','A', ... )
*>    29= | I - U U' | / ( n ulp )              SSYEVX('L','V','A', ... )
*>    30= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVX('L','N','A', ... )
*>    31= | A - U S U' | / ( |A| n ulp )        SSYEVX('L','V','I', ... )
*>    32= | I - U U' | / ( n ulp )              SSYEVX('L','V','I', ... )
*>    33= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVX('L','N','I', ... )
*>    34= | A - U S U' | / ( |A| n ulp )        SSYEVX('L','V','V', ... )
*>    35= | I - U U' | / ( n ulp )              SSYEVX('L','V','V', ... )
*>    36= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVX('L','N','V', ... )
*>    37= | A - U S U' | / ( |A| n ulp )        SSPEV('L','V', ... )
*>    38= | I - U U' | / ( n ulp )              SSPEV('L','V', ... )
*>    39= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEV('L','N', ... )
*>    40= | A - U S U' | / ( |A| n ulp )        SSPEVX('L','V','A', ... )
*>    41= | I - U U' | / ( n ulp )              SSPEVX('L','V','A', ... )
*>    42= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVX('L','N','A', ... )
*>    43= | A - U S U' | / ( |A| n ulp )        SSPEVX('L','V','I', ... )
*>    44= | I - U U' | / ( n ulp )              SSPEVX('L','V','I', ... )
*>    45= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVX('L','N','I', ... )
*>    46= | A - U S U' | / ( |A| n ulp )        SSPEVX('L','V','V', ... )
*>    47= | I - U U' | / ( n ulp )              SSPEVX('L','V','V', ... )
*>    48= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVX('L','N','V', ... )
*>    49= | A - U S U' | / ( |A| n ulp )        SSBEV('L','V', ... )
*>    50= | I - U U' | / ( n ulp )              SSBEV('L','V', ... )
*>    51= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEV('L','N', ... )
*>    52= | A - U S U' | / ( |A| n ulp )        SSBEVX('L','V','A', ... )
*>    53= | I - U U' | / ( n ulp )              SSBEVX('L','V','A', ... )
*>    54= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVX('L','N','A', ... )
*>    55= | A - U S U' | / ( |A| n ulp )        SSBEVX('L','V','I', ... )
*>    56= | I - U U' | / ( n ulp )              SSBEVX('L','V','I', ... )
*>    57= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVX('L','N','I', ... )
*>    58= | A - U S U' | / ( |A| n ulp )        SSBEVX('L','V','V', ... )
*>    59= | I - U U' | / ( n ulp )              SSBEVX('L','V','V', ... )
*>    60= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVX('L','N','V', ... )
*>    61= | A - U S U' | / ( |A| n ulp )        SSYEVD('L','V', ... )
*>    62= | I - U U' | / ( n ulp )              SSYEVD('L','V', ... )
*>    63= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVD('L','N', ... )
*>    64= | A - U S U' | / ( |A| n ulp )        SSPEVD('L','V', ... )
*>    65= | I - U U' | / ( n ulp )              SSPEVD('L','V', ... )
*>    66= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVD('L','N', ... )
*>    67= | A - U S U' | / ( |A| n ulp )        SSBEVD('L','V', ... )
*>    68= | I - U U' | / ( n ulp )              SSBEVD('L','V', ... )
*>    69= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVD('L','N', ... )
*>    70= | A - U S U' | / ( |A| n ulp )        SSYEVR('L','V','A', ... )
*>    71= | I - U U' | / ( n ulp )              SSYEVR('L','V','A', ... )
*>    72= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVR('L','N','A', ... )
*>    73= | A - U S U' | / ( |A| n ulp )        SSYEVR('L','V','I', ... )
*>    74= | I - U U' | / ( n ulp )              SSYEVR('L','V','I', ... )
*>    75= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVR('L','N','I', ... )
*>    76= | A - U S U' | / ( |A| n ulp )        SSYEVR('L','V','V', ... )
*>    77= | I - U U' | / ( n ulp )              SSYEVR('L','V','V', ... )
*>    78= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSYEVR('L','N','V', ... )
*>
*>    Tests 25 through 78 are repeated (as tests 79 through 132)
*>    with UPLO='U'
*>
*>    To be added in 1999
*>
*>    79= | A - U S U' | / ( |A| n ulp )        SSPEVR('L','V','A', ... )
*>    80= | I - U U' | / ( n ulp )              SSPEVR('L','V','A', ... )
*>    81= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVR('L','N','A', ... )
*>    82= | A - U S U' | / ( |A| n ulp )        SSPEVR('L','V','I', ... )
*>    83= | I - U U' | / ( n ulp )              SSPEVR('L','V','I', ... )
*>    84= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVR('L','N','I', ... )
*>    85= | A - U S U' | / ( |A| n ulp )        SSPEVR('L','V','V', ... )
*>    86= | I - U U' | / ( n ulp )              SSPEVR('L','V','V', ... )
*>    87= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSPEVR('L','N','V', ... )
*>    88= | A - U S U' | / ( |A| n ulp )        SSBEVR('L','V','A', ... )
*>    89= | I - U U' | / ( n ulp )              SSBEVR('L','V','A', ... )
*>    90= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVR('L','N','A', ... )
*>    91= | A - U S U' | / ( |A| n ulp )        SSBEVR('L','V','I', ... )
*>    92= | I - U U' | / ( n ulp )              SSBEVR('L','V','I', ... )
*>    93= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVR('L','N','I', ... )
*>    94= | A - U S U' | / ( |A| n ulp )        SSBEVR('L','V','V', ... )
*>    95= | I - U U' | / ( n ulp )              SSBEVR('L','V','V', ... )
*>    96= |D(with Z) - D(w/o Z)| / (|D| ulp)    SSBEVR('L','N','V', ... )
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
      SUBROUTINE SDRVKT( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $                   NOUNIT, A, LDA, D1, D2, D3, D4, EVEIGS, WA1,
     $                   WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK,
     $                   IWORK, LIWORK, RESULT, INFO )
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
      REAL               A( LDA, * ), D1( * ), D2( * ), D3( * ),
     $                   D4( * ), EVEIGS( * ), RESULT( * ), TAU( * ),
     $                   U( LDU, * ), V( LDU, * ), WA1( * ), WA2( * ),
     $                   WA3( * ), WORK( * ), Z( LDU, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO, TEN
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   TEN = 10.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 18 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADNN
      CHARACTER          UPLO
      INTEGER            I, IDIAG, IHBW, IINFO, IL, IMODE, IROW,
     $                   ITEMP, ITYPE, IU, IUPLO, J, J1, J2, JCOL,
     $                   JSIZE, JTYPE, LGN, LIWEDC, LWEDC,
     $                   MTYPES, N, NERRS, NMATS, NMAX, NTEST,
     $                   NTESTT
      REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL,
     $                   RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL,
     $                   VL, VU
*     ..
*     .. Local Arrays ..
      INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ),
     $                   ISEED3( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ),
     $                   KTYPE( MAXTYP )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLARND, SSXT1
      EXTERNAL           SLAMCH, SLARND, SSXT1
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALASVM, SLABAD, SLACPY, SLAFTS, SLASET, SLATMR,
     $                   SLATMS, SKTEV, SKTT21, SKYEV, SKYT21, XERBLA
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, REAL, SQRT
*     ..
*     .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8, 3*9 /
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1,
     $                   2, 3, 1, 2, 3 /
      DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0,
     $                   0, 0, 4, 4, 4 /
*     ..
*     .. Executable Statements ..
*
*     Keep ftrnchek happy
*
      VL = ZERO
      VU = ZERO
*
*     1)      Check for errors
*
      NTESTT = 0
      INFO = 0
*
      BADNN = .FALSE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 )
     $      BADNN = .TRUE.
   10 CONTINUE
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
         INFO = -16
      ELSE IF( 2*MAX( 2, NMAX )**2.GT.LWORK ) THEN
         INFO = -21
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SDRVKT', -INFO )
         RETURN
      END IF
*
*     Quick return if nothing to do
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 )
     $   RETURN
*
*     More Important constants
*
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = SLAMCH( 'Overflow' )
      CALL SLABAD( UNFL, OVFL )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
*
*     Loop over sizes, types
*
      DO 20 I = 1, 4
         ISEED2( I ) = ISEED( I )
         ISEED3( I ) = ISEED( I )
   20 CONTINUE
*
      NERRS = 0
      NMATS = 0
*
*
      DO 1740 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         IF( N.GT.0 ) THEN
            LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
            IF( 2**LGN.LT.N )
     $         LGN = LGN + 1
            IF( 2**LGN.LT.N )
     $         LGN = LGN + 1
            LWEDC = 1 + 4*N + 2*N*LGN + 4*N**2
c           LIWEDC = 6 + 6*N + 5*N*LGN
            LIWEDC = 3 + 5*N
         ELSE
            LWEDC = 9
c           LIWEDC = 12
            LIWEDC = 8
         END IF
         ANINV = ONE / REAL( MAX( 1, N ) )
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 1730 JTYPE = 1, MTYPES
*
            IF( .NOT.DOTYPE( JTYPE ) )
     $         GO TO 1730
            NMATS = NMATS + 1
            NTEST = 0
*
            DO 30 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   30       CONTINUE
*
*           2)      Compute "A"
*
*                   Control parameters:
*
*               KMAGN  KMODE        KTYPE
*           =1  O(1)   clustered 1  zero
*           =2  large  clustered 2  identity
*           =3  small  exponential  (none)
*           =4         arithmetic   diagonal, (w/ eigenvalues)
*           =5         random log   skew-symmetric, w/ eigenvalues
*           =6         random       (none)
*           =7                      random diagonal
*           =8                      random skew-symmetric
*           =9                      band skew-symmetric, w/ eigenvalues
*
            IF( MTYPES.GT.MAXTYP )
     $         GO TO 110
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
            COND = ULPINV
*
*           Special Matrices -- Identity & Jordan block
*
*                   Zero
*
            IF( ITYPE.EQ.1 ) THEN
               IINFO = 0
*
            ELSE IF( ITYPE.EQ.2 ) THEN
*
*              Identity
*
               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
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
            ELSE IF( ITYPE.EQ.5 ) THEN
*
*              skew-symmetric, eigenvalues specified
*
               CALL SLATMS( N, N, 'S', ISEED, 'K', WORK, IMODE, COND,
     $                      ANORM, N, N, 'N', A, LDA, WORK( N+1 ),
     $                      IINFO )
*
            ELSE IF( ITYPE.EQ.7 ) THEN
*
*              tridiagonal, random eigenvalues
*
               IDUMMA( 1 ) = 1
               CALL SLATMR( N, N, 'S', ISEED, 'K', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 1, 1,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.8 ) THEN
*
*              skew-symmetric, random eigenvalues
*
               IDUMMA( 1 ) = 1
               CALL SLATMR( N, N, 'S', ISEED, 'K', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.9 ) THEN
*
*              skew-symmetric banded, eigenvalues specified
*
               IHBW = INT( ( N-1 )*SLARND( 1, ISEED3 ) )
               CALL SLATMS( N, N, 'S', ISEED, 'K', WORK, IMODE, COND,
     $                      ANORM, IHBW, IHBW, 'Z', U, LDU, WORK( N+1 ),
     $                      IINFO )
*
*              Store as dense matrix for most routines.
*
               CALL SLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
               DO 100 IDIAG = -IHBW, IHBW
                  IROW = IHBW - IDIAG + 1
                  J1 = MAX( 1, IDIAG+1 )
                  J2 = MIN( N, N+IDIAG )
                  DO 90 J = J1, J2
                     I = J - IDIAG
                     A( I, J ) = U( IROW, J )
   90             CONTINUE
  100          CONTINUE
            ELSE
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
  110       CONTINUE
*
            ABSTOL = UNFL + UNFL
            IF( N.LE.1 ) THEN
               IL = 1
               IU = N
            ELSE
               IL = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               IF( IL.GT.IU ) THEN
                  ITEMP = IL
                  IL = IU
                  IU = ITEMP
               END IF
            END IF
*
*           3)      If matrix is tridiagonal, call SKTEV and SSTEVX.
*
            IF( JTYPE.LE.7 ) THEN
               NTEST = 1
               DO 120 I = 1, N
                  D1( I ) = REAL( A( I, I ) )
  120          CONTINUE
               DO 130 I = 1, N - 1
                  D2( I ) = REAL( A( I+1, I ) )
  130          CONTINUE
               SRNAMT = 'SKTEV'
               CALL SKTEV( 'V', N, D1, D2, Z, LDU, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SKTEV(V)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 1 ) = ULPINV
                     RESULT( 2 ) = ULPINV
                     RESULT( 3 ) = ULPINV
                     GO TO 180
                  END IF
               END IF
*
*              Do tests 1 and 2.
*
               DO 140 I = 1, N
                  D3( I ) = REAL( A( I, I ) )
  140          CONTINUE
               DO 150 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  150          CONTINUE
               CALL SKTT21( N, 1, D3, D4, D2, D1, Z, LDU, WORK,
     $                      RESULT( 1 ) )
*
               NTEST = 3
               DO 160 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  160          CONTINUE
               SRNAMT = 'SKTEV'
               CALL SKTEV( 'N', N, D3, D4, Z, LDU, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SKTEV(N)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 3 ) = ULPINV
                     GO TO 180
                  END IF
               END IF
*
*              Do test 3.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 170 J = 1, N-1
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  170          CONTINUE
               RESULT( 3 ) = TEMP2 / MAX( UNFL,
     $                       ULP*MAX( TEMP1, TEMP2 ) )
*
  180          CONTINUE
*
            ELSE
*
               DO 640 I = 1, 3
                  RESULT( I ) = ZERO
  640          CONTINUE
               NTEST = 3
            END IF
*
*           Perform remaining tests storing upper or lower triangular
*           part of matrix.
*
            DO 1720 IUPLO = 0, 1
               IF( IUPLO.EQ.0 ) THEN
                  UPLO = 'L'
               ELSE
                  UPLO = 'U'
               END IF
*
*              4)      Call SKYEV and SSYEVX.
*
               CALL SLACPY( ' ', N, N, A, LDA, V, LDU )
*
               NTEST = NTEST + 1
               SRNAMT = 'SKYEV'
               CALL SKYEV( 'V', UPLO, N, A, LDU, D1, WORK, LWORK,
     $                     IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SKYEV(V,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 660
                  END IF
               END IF
*
*              Do tests 25 and 26 (or +54)
*
               CALL SKYT21( 1, UPLO, N, 1, V, LDU, D2, D1, A, LDU, Z,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 2
               SRNAMT = 'SKYEV'
               CALL SKYEV( 'N', UPLO, N, A, LDU, D3, WORK, LWORK,
     $                     IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SKYEV(N,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 660
                  END IF
               END IF
*
*              Do test 27 (or +54)
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 650 J = 1, N-1
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  650          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
  660          CONTINUE
*
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
*
 1720       CONTINUE
*
*           End of Loop -- Check for RESULT(j) > THRESH
*
            NTESTT = NTESTT + NTEST
*
            CALL SLAFTS( 'SKT', N, N, JTYPE, NTEST, RESULT, IOLDSD,
     $                   THRESH, NOUNIT, NERRS )
*
 1730    CONTINUE
 1740 CONTINUE
*
*     Summary
*
      CALL ALASVM( 'SKT', NOUNIT, NERRS, NTESTT, 0 )
*
 9999 FORMAT( ' SDRVKT: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of SDRVKT
*
      END
