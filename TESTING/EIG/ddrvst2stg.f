*> \brief \b DDRVST2STG
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DDRVST2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
*                          NOUNIT, A, LDA, D1, D2, D3, D4, EVEIGS, WA1,
*                          WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK,
*                          IWORK, LIWORK, RESULT, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES,
*      $                   NTYPES
*       DOUBLE PRECISION   THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
*       DOUBLE PRECISION   A( LDA, * ), D1( * ), D2( * ), D3( * ),
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
*>      DDRVST2STG  checks the symmetric eigenvalue problem drivers.
*>
*>              DSTEV computes all eigenvalues and, optionally,
*>              eigenvectors of a real symmetric tridiagonal matrix.
*>
*>              DSTEVX computes selected eigenvalues and, optionally,
*>              eigenvectors of a real symmetric tridiagonal matrix.
*>
*>              DSTEVR computes selected eigenvalues and, optionally,
*>              eigenvectors of a real symmetric tridiagonal matrix
*>              using the Relatively Robust Representation where it can.
*>
*>              DSYEV computes all eigenvalues and, optionally,
*>              eigenvectors of a real symmetric matrix.
*>
*>              DSYEVX computes selected eigenvalues and, optionally,
*>              eigenvectors of a real symmetric matrix.
*>
*>              DSYEVR computes selected eigenvalues and, optionally,
*>              eigenvectors of a real symmetric matrix
*>              using the Relatively Robust Representation where it can.
*>
*>              DSPEV computes all eigenvalues and, optionally,
*>              eigenvectors of a real symmetric matrix in packed
*>              storage.
*>
*>              DSPEVX computes selected eigenvalues and, optionally,
*>              eigenvectors of a real symmetric matrix in packed
*>              storage.
*>
*>              DSBEV computes all eigenvalues and, optionally,
*>              eigenvectors of a real symmetric band matrix.
*>
*>              DSBEVX computes selected eigenvalues and, optionally,
*>              eigenvectors of a real symmetric band matrix.
*>
*>              DSYEVD computes all eigenvalues and, optionally,
*>              eigenvectors of a real symmetric matrix using
*>              a divide and conquer algorithm.
*>
*>              DSPEVD computes all eigenvalues and, optionally,
*>              eigenvectors of a real symmetric matrix in packed
*>              storage, using a divide and conquer algorithm.
*>
*>              DSBEVD computes all eigenvalues and, optionally,
*>              eigenvectors of a real symmetric band matrix,
*>              using a divide and conquer algorithm.
*>
*>      When DDRVST2STG is called, a number of matrix "sizes" ("n's") and a
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
*>      (13) Symmetric matrix with random entries chosen from (-1,1).
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
*>          DDRVST2STG does nothing.  It must be at least zero.
*>          Not modified.
*>
*>  NN      INTEGER array, dimension (NSIZES)
*>          An array containing the sizes to be used for the matrices.
*>          Zero values will be skipped.  The values must be at least
*>          zero.
*>          Not modified.
*>
*>  NTYPES  INTEGER
*>          The number of elements in DOTYPE.   If it is zero, DDRVST2STG
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
*>          next call to DDRVST2STG to continue the same random number
*>          sequence.
*>          Modified.
*>
*>  THRESH  DOUBLE PRECISION
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
*>  A       DOUBLE PRECISION array, dimension (LDA , max(NN))
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
*>  D1      DOUBLE PRECISION array, dimension (max(NN))
*>          The eigenvalues of A, as computed by DSTEQR simultaneously
*>          with Z.  On exit, the eigenvalues in D1 correspond with the
*>          matrix in A.
*>          Modified.
*>
*>  D2      DOUBLE PRECISION array, dimension (max(NN))
*>          The eigenvalues of A, as computed by DSTEQR if Z is not
*>          computed.  On exit, the eigenvalues in D2 correspond with
*>          the matrix in A.
*>          Modified.
*>
*>  D3      DOUBLE PRECISION array, dimension (max(NN))
*>          The eigenvalues of A, as computed by DSTERF.  On exit, the
*>          eigenvalues in D3 correspond with the matrix in A.
*>          Modified.
*>
*>  D4      DOUBLE PRECISION array, dimension
*>
*>  EVEIGS  DOUBLE PRECISION array, dimension (max(NN))
*>          The eigenvalues as computed by DSTEV('N', ... )
*>          (I reserve the right to change this to the output of
*>          whichever algorithm computes the most accurate eigenvalues).
*>
*>  WA1     DOUBLE PRECISION array, dimension
*>
*>  WA2     DOUBLE PRECISION array, dimension
*>
*>  WA3     DOUBLE PRECISION array, dimension
*>
*>  U       DOUBLE PRECISION array, dimension (LDU, max(NN))
*>          The orthogonal matrix computed by DSYTRD + DORGTR.
*>          Modified.
*>
*>  LDU     INTEGER
*>          The leading dimension of U, Z, and V.  It must be at
*>          least 1 and at least max( NN ).
*>          Not modified.
*>
*>  V       DOUBLE PRECISION array, dimension (LDU, max(NN))
*>          The Housholder vectors computed by DSYTRD in reducing A to
*>          tridiagonal form.
*>          Modified.
*>
*>  TAU     DOUBLE PRECISION array, dimension (max(NN))
*>          The Householder factors computed by DSYTRD in reducing A
*>          to tridiagonal form.
*>          Modified.
*>
*>  Z       DOUBLE PRECISION array, dimension (LDU, max(NN))
*>          The orthogonal matrix of eigenvectors computed by DSTEQR,
*>          DPTEQR, and DSTEIN.
*>          Modified.
*>
*>  WORK    DOUBLE PRECISION array, dimension (LWORK)
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
*>  RESULT  DOUBLE PRECISION array, dimension (105)
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
*>          If  DLATMR, DLATMS, DSYTRD, DORGTR, DSTEQR, DSTERF,
*>              or DORMTR returns an error code, the
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
*>                       so far (computed by DLAFTS).
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
*>    1= | A - U S U' | / ( |A| n ulp )         DSTEV('V', ... )
*>    2= | I - U U' | / ( n ulp )               DSTEV('V', ... )
*>    3= |D(with Z) - D(w/o Z)| / (|D| ulp)     DSTEV('N', ... )
*>    4= | A - U S U' | / ( |A| n ulp )         DSTEVX('V','A', ... )
*>    5= | I - U U' | / ( n ulp )               DSTEVX('V','A', ... )
*>    6= |D(with Z) - EVEIGS| / (|D| ulp)       DSTEVX('N','A', ... )
*>    7= | A - U S U' | / ( |A| n ulp )         DSTEVR('V','A', ... )
*>    8= | I - U U' | / ( n ulp )               DSTEVR('V','A', ... )
*>    9= |D(with Z) - EVEIGS| / (|D| ulp)       DSTEVR('N','A', ... )
*>    10= | A - U S U' | / ( |A| n ulp )        DSTEVX('V','I', ... )
*>    11= | I - U U' | / ( n ulp )              DSTEVX('V','I', ... )
*>    12= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSTEVX('N','I', ... )
*>    13= | A - U S U' | / ( |A| n ulp )        DSTEVX('V','V', ... )
*>    14= | I - U U' | / ( n ulp )              DSTEVX('V','V', ... )
*>    15= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSTEVX('N','V', ... )
*>    16= | A - U S U' | / ( |A| n ulp )        DSTEVD('V', ... )
*>    17= | I - U U' | / ( n ulp )              DSTEVD('V', ... )
*>    18= |D(with Z) - EVEIGS| / (|D| ulp)      DSTEVD('N', ... )
*>    19= | A - U S U' | / ( |A| n ulp )        DSTEVR('V','I', ... )
*>    20= | I - U U' | / ( n ulp )              DSTEVR('V','I', ... )
*>    21= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSTEVR('N','I', ... )
*>    22= | A - U S U' | / ( |A| n ulp )        DSTEVR('V','V', ... )
*>    23= | I - U U' | / ( n ulp )              DSTEVR('V','V', ... )
*>    24= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSTEVR('N','V', ... )
*>
*>    25= | A - U S U' | / ( |A| n ulp )        DSYEV('L','V', ... )
*>    26= | I - U U' | / ( n ulp )              DSYEV('L','V', ... )
*>    27= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEV_2STAGE('L','N', ... )
*>    28= | A - U S U' | / ( |A| n ulp )        DSYEVX('L','V','A', ... )
*>    29= | I - U U' | / ( n ulp )              DSYEVX('L','V','A', ... )
*>    30= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVX_2STAGE('L','N','A', ... )
*>    31= | A - U S U' | / ( |A| n ulp )        DSYEVX('L','V','I', ... )
*>    32= | I - U U' | / ( n ulp )              DSYEVX('L','V','I', ... )
*>    33= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVX_2STAGE('L','N','I', ... )
*>    34= | A - U S U' | / ( |A| n ulp )        DSYEVX('L','V','V', ... )
*>    35= | I - U U' | / ( n ulp )              DSYEVX('L','V','V', ... )
*>    36= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVX_2STAGE('L','N','V', ... )
*>    37= | A - U S U' | / ( |A| n ulp )        DSPEV('L','V', ... )
*>    38= | I - U U' | / ( n ulp )              DSPEV('L','V', ... )
*>    39= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEV('L','N', ... )
*>    40= | A - U S U' | / ( |A| n ulp )        DSPEVX('L','V','A', ... )
*>    41= | I - U U' | / ( n ulp )              DSPEVX('L','V','A', ... )
*>    42= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVX('L','N','A', ... )
*>    43= | A - U S U' | / ( |A| n ulp )        DSPEVX('L','V','I', ... )
*>    44= | I - U U' | / ( n ulp )              DSPEVX('L','V','I', ... )
*>    45= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVX('L','N','I', ... )
*>    46= | A - U S U' | / ( |A| n ulp )        DSPEVX('L','V','V', ... )
*>    47= | I - U U' | / ( n ulp )              DSPEVX('L','V','V', ... )
*>    48= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVX('L','N','V', ... )
*>    49= | A - U S U' | / ( |A| n ulp )        DSBEV('L','V', ... )
*>    50= | I - U U' | / ( n ulp )              DSBEV('L','V', ... )
*>    51= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEV_2STAGE('L','N', ... )
*>    52= | A - U S U' | / ( |A| n ulp )        DSBEVX('L','V','A', ... )
*>    53= | I - U U' | / ( n ulp )              DSBEVX('L','V','A', ... )
*>    54= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVX_2STAGE('L','N','A', ... )
*>    55= | A - U S U' | / ( |A| n ulp )        DSBEVX('L','V','I', ... )
*>    56= | I - U U' | / ( n ulp )              DSBEVX('L','V','I', ... )
*>    57= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVX_2STAGE('L','N','I', ... )
*>    58= | A - U S U' | / ( |A| n ulp )        DSBEVX('L','V','V', ... )
*>    59= | I - U U' | / ( n ulp )              DSBEVX('L','V','V', ... )
*>    60= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVX_2STAGE('L','N','V', ... )
*>    61= | A - U S U' | / ( |A| n ulp )        DSYEVD('L','V', ... )
*>    62= | I - U U' | / ( n ulp )              DSYEVD('L','V', ... )
*>    63= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVD_2STAGE('L','N', ... )
*>    64= | A - U S U' | / ( |A| n ulp )        DSPEVD('L','V', ... )
*>    65= | I - U U' | / ( n ulp )              DSPEVD('L','V', ... )
*>    66= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVD('L','N', ... )
*>    67= | A - U S U' | / ( |A| n ulp )        DSBEVD('L','V', ... )
*>    68= | I - U U' | / ( n ulp )              DSBEVD('L','V', ... )
*>    69= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVD_2STAGE('L','N', ... )
*>    70= | A - U S U' | / ( |A| n ulp )        DSYEVR('L','V','A', ... )
*>    71= | I - U U' | / ( n ulp )              DSYEVR('L','V','A', ... )
*>    72= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVR_2STAGE('L','N','A', ... )
*>    73= | A - U S U' | / ( |A| n ulp )        DSYEVR('L','V','I', ... )
*>    74= | I - U U' | / ( n ulp )              DSYEVR('L','V','I', ... )
*>    75= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVR_2STAGE('L','N','I', ... )
*>    76= | A - U S U' | / ( |A| n ulp )        DSYEVR('L','V','V', ... )
*>    77= | I - U U' | / ( n ulp )              DSYEVR('L','V','V', ... )
*>    78= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVR_2STAGE('L','N','V', ... )
*>
*>    Tests 25 through 78 are repeated (as tests 79 through 132)
*>    with UPLO='U'
*>
*>    To be added in 1999
*>
*>    79= | A - U S U' | / ( |A| n ulp )        DSPEVR('L','V','A', ... )
*>    80= | I - U U' | / ( n ulp )              DSPEVR('L','V','A', ... )
*>    81= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVR('L','N','A', ... )
*>    82= | A - U S U' | / ( |A| n ulp )        DSPEVR('L','V','I', ... )
*>    83= | I - U U' | / ( n ulp )              DSPEVR('L','V','I', ... )
*>    84= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVR('L','N','I', ... )
*>    85= | A - U S U' | / ( |A| n ulp )        DSPEVR('L','V','V', ... )
*>    86= | I - U U' | / ( n ulp )              DSPEVR('L','V','V', ... )
*>    87= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVR('L','N','V', ... )
*>    88= | A - U S U' | / ( |A| n ulp )        DSBEVR('L','V','A', ... )
*>    89= | I - U U' | / ( n ulp )              DSBEVR('L','V','A', ... )
*>    90= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVR('L','N','A', ... )
*>    91= | A - U S U' | / ( |A| n ulp )        DSBEVR('L','V','I', ... )
*>    92= | I - U U' | / ( n ulp )              DSBEVR('L','V','I', ... )
*>    93= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVR('L','N','I', ... )
*>    94= | A - U S U' | / ( |A| n ulp )        DSBEVR('L','V','V', ... )
*>    95= | I - U U' | / ( n ulp )              DSBEVR('L','V','V', ... )
*>    96= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVR('L','N','V', ... )
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
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE DDRVST2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $                   NOUNIT, A, LDA, D1, D2, D3, D4, EVEIGS, WA1,
     $                   WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK,
     $                   IWORK, LIWORK, RESULT, INFO )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES,
     $                   NTYPES
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
      DOUBLE PRECISION   A( LDA, * ), D1( * ), D2( * ), D3( * ),
     $                   D4( * ), EVEIGS( * ), RESULT( * ), TAU( * ),
     $                   U( LDU, * ), V( LDU, * ), WA1( * ), WA2( * ),
     $                   WA3( * ), WORK( * ), Z( LDU, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, TEN
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   TEN = 10.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 18 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADNN
      CHARACTER          UPLO
      INTEGER            I, IDIAG, IHBW, IINFO, IL, IMODE, INDX, IROW,
     $                   ITEMP, ITYPE, IU, IUPLO, J, J1, J2, JCOL,
     $                   JSIZE, JTYPE, KD, LGN, LIWEDC, LWEDC, M, M2,
     $                   M3, MTYPES, N, NERRS, NMATS, NMAX, NTEST,
     $                   NTESTT
      DOUBLE PRECISION   ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL,
     $                   RTUNFL, TEMP1, TEMP2, TEMP3, ULP, ULPINV, UNFL,
     $                   VL, VU
*     ..
*     .. Local Arrays ..
      INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ),
     $                   ISEED3( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ),
     $                   KTYPE( MAXTYP )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLARND, DSXT1
      EXTERNAL           DLAMCH, DLARND, DSXT1
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALASVM, DLACPY, DLAFTS, DLASET, DLATMR,
     $                   DLATMS, DSBEV, DSBEVD, DSBEVX, DSPEV, DSPEVD,
     $                   DSPEVX, DSTEV, DSTEVD, DSTEVR, DSTEVX, DSTT21,
     $                   DSTT22, DSYEV, DSYEVD, DSYEVR, DSYEVX, DSYT21,
     $                   DSYEVD_2STAGE, DSYEVR_2STAGE, DSYEVX_2STAGE,
     $                   DSYEV_2STAGE, DSBEV_2STAGE, DSBEVD_2STAGE,
     $                   DSBEVX_2STAGE, DSYTRD_2STAGE, DSYTRD_SY2SB, 
     $                   DSYTRD_SB2ST, DSYT22, XERBLA
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX, MIN, SQRT
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
         CALL XERBLA( 'DDRVST2STG', -INFO )
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
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
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
            LGN = INT( LOG( DBLE( N ) ) / LOG( TWO ) )
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
         ANINV = ONE / DBLE( MAX( 1, N ) )
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
*           =5         random log   symmetric, w/ eigenvalues
*           =6         random       (none)
*           =7                      random diagonal
*           =8                      random symmetric
*           =9                      band symmetric, w/ eigenvalues
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
            CALL DLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
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
*              Diagonal Matrix, [Eigen]values Specified
*
               CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND,
     $                      ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ),
     $                      IINFO )
*
            ELSE IF( ITYPE.EQ.5 ) THEN
*
*              Symmetric, eigenvalues specified
*
               CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND,
     $                      ANORM, N, N, 'N', A, LDA, WORK( N+1 ),
     $                      IINFO )
*
            ELSE IF( ITYPE.EQ.7 ) THEN
*
*              Diagonal, random eigenvalues
*
               IDUMMA( 1 ) = 1
               CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.8 ) THEN
*
*              Symmetric, random eigenvalues
*
               IDUMMA( 1 ) = 1
               CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.9 ) THEN
*
*              Symmetric banded, eigenvalues specified
*
               IHBW = INT( ( N-1 )*DLARND( 1, ISEED3 ) )
               CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND,
     $                      ANORM, IHBW, IHBW, 'Z', U, LDU, WORK( N+1 ),
     $                      IINFO )
*
*              Store as dense matrix for most routines.
*
               CALL DLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
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
               IL = 1 + INT( ( N-1 )*DLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( N-1 )*DLARND( 1, ISEED2 ) )
               IF( IL.GT.IU ) THEN
                  ITEMP = IL
                  IL = IU
                  IU = ITEMP
               END IF
            END IF
*
*           3)      If matrix is tridiagonal, call DSTEV and DSTEVX.
*
            IF( JTYPE.LE.7 ) THEN
               NTEST = 1
               DO 120 I = 1, N
                  D1( I ) = DBLE( A( I, I ) )
  120          CONTINUE
               DO 130 I = 1, N - 1
                  D2( I ) = DBLE( A( I+1, I ) )
  130          CONTINUE
               SRNAMT = 'DSTEV'
               CALL DSTEV( 'V', N, D1, D2, Z, LDU, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEV(V)', IINFO, N,
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
                  D3( I ) = DBLE( A( I, I ) )
  140          CONTINUE
               DO 150 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  150          CONTINUE
               CALL DSTT21( N, 0, D3, D4, D1, D2, Z, LDU, WORK,
     $                      RESULT( 1 ) )
*
               NTEST = 3
               DO 160 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  160          CONTINUE
               SRNAMT = 'DSTEV'
               CALL DSTEV( 'N', N, D3, D4, Z, LDU, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEV(N)', IINFO, N,
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
               DO 170 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  170          CONTINUE
               RESULT( 3 ) = TEMP2 / MAX( UNFL,
     $                       ULP*MAX( TEMP1, TEMP2 ) )
*
  180          CONTINUE
*
               NTEST = 4
               DO 190 I = 1, N
                  EVEIGS( I ) = D3( I )
                  D1( I ) = DBLE( A( I, I ) )
  190          CONTINUE
               DO 200 I = 1, N - 1
                  D2( I ) = DBLE( A( I+1, I ) )
  200          CONTINUE
               SRNAMT = 'DSTEVX'
               CALL DSTEVX( 'V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL,
     $                      M, WA1, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ),
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVX(V,A)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 4 ) = ULPINV
                     RESULT( 5 ) = ULPINV
                     RESULT( 6 ) = ULPINV
                     GO TO 250
                  END IF
               END IF
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
*
*              Do tests 4 and 5.
*
               DO 210 I = 1, N
                  D3( I ) = DBLE( A( I, I ) )
  210          CONTINUE
               DO 220 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  220          CONTINUE
               CALL DSTT21( N, 0, D3, D4, WA1, D2, Z, LDU, WORK,
     $                      RESULT( 4 ) )
*
               NTEST = 6
               DO 230 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  230          CONTINUE
               SRNAMT = 'DSTEVX'
               CALL DSTEVX( 'N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL,
     $                      M2, WA2, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVX(N,A)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 6 ) = ULPINV
                     GO TO 250
                  END IF
               END IF
*
*              Do test 6.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 240 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA2( J ) ),
     $                    ABS( EVEIGS( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA2( J )-EVEIGS( J ) ) )
  240          CONTINUE
               RESULT( 6 ) = TEMP2 / MAX( UNFL,
     $                       ULP*MAX( TEMP1, TEMP2 ) )
*
  250          CONTINUE
*
               NTEST = 7
               DO 260 I = 1, N
                  D1( I ) = DBLE( A( I, I ) )
  260          CONTINUE
               DO 270 I = 1, N - 1
                  D2( I ) = DBLE( A( I+1, I ) )
  270          CONTINUE
               SRNAMT = 'DSTEVR'
               CALL DSTEVR( 'V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL,
     $                      M, WA1, Z, LDU, IWORK, WORK, LWORK,
     $                      IWORK(2*N+1), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVR(V,A)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 7 ) = ULPINV
                     RESULT( 8 ) = ULPINV
                     GO TO 320
                  END IF
               END IF
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
*
*              Do tests 7 and 8.
*
               DO 280 I = 1, N
                  D3( I ) = DBLE( A( I, I ) )
  280          CONTINUE
               DO 290 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  290          CONTINUE
               CALL DSTT21( N, 0, D3, D4, WA1, D2, Z, LDU, WORK,
     $                      RESULT( 7 ) )
*
               NTEST = 9
               DO 300 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  300          CONTINUE
               SRNAMT = 'DSTEVR'
               CALL DSTEVR( 'N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL,
     $                      M2, WA2, Z, LDU, IWORK, WORK, LWORK,
     $                      IWORK(2*N+1), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVR(N,A)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 9 ) = ULPINV
                     GO TO 320
                  END IF
               END IF
*
*              Do test 9.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 310 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA2( J ) ),
     $                    ABS( EVEIGS( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA2( J )-EVEIGS( J ) ) )
  310          CONTINUE
               RESULT( 9 ) = TEMP2 / MAX( UNFL,
     $                       ULP*MAX( TEMP1, TEMP2 ) )
*
  320          CONTINUE
*
*
               NTEST = 10
               DO 330 I = 1, N
                  D1( I ) = DBLE( A( I, I ) )
  330          CONTINUE
               DO 340 I = 1, N - 1
                  D2( I ) = DBLE( A( I+1, I ) )
  340          CONTINUE
               SRNAMT = 'DSTEVX'
               CALL DSTEVX( 'V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL,
     $                      M2, WA2, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVX(V,I)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 10 ) = ULPINV
                     RESULT( 11 ) = ULPINV
                     RESULT( 12 ) = ULPINV
                     GO TO 380
                  END IF
               END IF
*
*              Do tests 10 and 11.
*
               DO 350 I = 1, N
                  D3( I ) = DBLE( A( I, I ) )
  350          CONTINUE
               DO 360 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  360          CONTINUE
               CALL DSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
     $                      MAX( 1, M2 ), RESULT( 10 ) )
*
*
               NTEST = 12
               DO 370 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  370          CONTINUE
               SRNAMT = 'DSTEVX'
               CALL DSTEVX( 'N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL,
     $                      M3, WA3, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVX(N,I)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 12 ) = ULPINV
                     GO TO 380
                  END IF
               END IF
*
*              Do test 12.
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 12 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )
*
  380          CONTINUE
*
               NTEST = 12
               IF( N.GT.0 ) THEN
                  IF( IL.NE.1 ) THEN
                     VL = WA1( IL ) - MAX( HALF*
     $                    ( WA1( IL )-WA1( IL-1 ) ), TEN*ULP*TEMP3,
     $                    TEN*RTUNFL )
                  ELSE
                     VL = WA1( 1 ) - MAX( HALF*( WA1( N )-WA1( 1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
                  IF( IU.NE.N ) THEN
                     VU = WA1( IU ) + MAX( HALF*
     $                    ( WA1( IU+1 )-WA1( IU ) ), TEN*ULP*TEMP3,
     $                    TEN*RTUNFL )
                  ELSE
                     VU = WA1( N ) + MAX( HALF*( WA1( N )-WA1( 1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
               ELSE
                  VL = ZERO
                  VU = ONE
               END IF
*
               DO 390 I = 1, N
                  D1( I ) = DBLE( A( I, I ) )
  390          CONTINUE
               DO 400 I = 1, N - 1
                  D2( I ) = DBLE( A( I+1, I ) )
  400          CONTINUE
               SRNAMT = 'DSTEVX'
               CALL DSTEVX( 'V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL,
     $                      M2, WA2, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVX(V,V)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 13 ) = ULPINV
                     RESULT( 14 ) = ULPINV
                     RESULT( 15 ) = ULPINV
                     GO TO 440
                  END IF
               END IF
*
               IF( M2.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( 13 ) = ULPINV
                  RESULT( 14 ) = ULPINV
                  RESULT( 15 ) = ULPINV
                  GO TO 440
               END IF
*
*              Do tests 13 and 14.
*
               DO 410 I = 1, N
                  D3( I ) = DBLE( A( I, I ) )
  410          CONTINUE
               DO 420 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  420          CONTINUE
               CALL DSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
     $                      MAX( 1, M2 ), RESULT( 13 ) )
*
               NTEST = 15
               DO 430 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  430          CONTINUE
               SRNAMT = 'DSTEVX'
               CALL DSTEVX( 'N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL,
     $                      M3, WA3, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVX(N,V)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 15 ) = ULPINV
                     GO TO 440
                  END IF
               END IF
*
*              Do test 15.
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 15 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
  440          CONTINUE
*
               NTEST = 16
               DO 450 I = 1, N
                  D1( I ) = DBLE( A( I, I ) )
  450          CONTINUE
               DO 460 I = 1, N - 1
                  D2( I ) = DBLE( A( I+1, I ) )
  460          CONTINUE
               SRNAMT = 'DSTEVD'
               CALL DSTEVD( 'V', N, D1, D2, Z, LDU, WORK, LWEDC, IWORK,
     $                      LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVD(V)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 16 ) = ULPINV
                     RESULT( 17 ) = ULPINV
                     RESULT( 18 ) = ULPINV
                     GO TO 510
                  END IF
               END IF
*
*              Do tests 16 and 17.
*
               DO 470 I = 1, N
                  D3( I ) = DBLE( A( I, I ) )
  470          CONTINUE
               DO 480 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  480          CONTINUE
               CALL DSTT21( N, 0, D3, D4, D1, D2, Z, LDU, WORK,
     $                      RESULT( 16 ) )
*
               NTEST = 18
               DO 490 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  490          CONTINUE
               SRNAMT = 'DSTEVD'
               CALL DSTEVD( 'N', N, D3, D4, Z, LDU, WORK, LWEDC, IWORK,
     $                      LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVD(N)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 18 ) = ULPINV
                     GO TO 510
                  END IF
               END IF
*
*              Do test 18.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 500 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( EVEIGS( J ) ),
     $                    ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( EVEIGS( J )-D3( J ) ) )
  500          CONTINUE
               RESULT( 18 ) = TEMP2 / MAX( UNFL,
     $                        ULP*MAX( TEMP1, TEMP2 ) )
*
  510          CONTINUE
*
               NTEST = 19
               DO 520 I = 1, N
                  D1( I ) = DBLE( A( I, I ) )
  520          CONTINUE
               DO 530 I = 1, N - 1
                  D2( I ) = DBLE( A( I+1, I ) )
  530          CONTINUE
               SRNAMT = 'DSTEVR'
               CALL DSTEVR( 'V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL,
     $                      M2, WA2, Z, LDU, IWORK, WORK, LWORK,
     $                      IWORK(2*N+1), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVR(V,I)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 19 ) = ULPINV
                     RESULT( 20 ) = ULPINV
                     RESULT( 21 ) = ULPINV
                     GO TO 570
                  END IF
               END IF
*
*              DO tests 19 and 20.
*
               DO 540 I = 1, N
                  D3( I ) = DBLE( A( I, I ) )
  540          CONTINUE
               DO 550 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  550          CONTINUE
               CALL DSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
     $                      MAX( 1, M2 ), RESULT( 19 ) )
*
*
               NTEST = 21
               DO 560 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  560          CONTINUE
               SRNAMT = 'DSTEVR'
               CALL DSTEVR( 'N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL,
     $                      M3, WA3, Z, LDU, IWORK, WORK, LWORK,
     $                      IWORK(2*N+1), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVR(N,I)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 21 ) = ULPINV
                     GO TO 570
                  END IF
               END IF
*
*              Do test 21.
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 21 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )
*
  570          CONTINUE
*
               NTEST = 21
               IF( N.GT.0 ) THEN
                  IF( IL.NE.1 ) THEN
                     VL = WA1( IL ) - MAX( HALF*
     $                    ( WA1( IL )-WA1( IL-1 ) ), TEN*ULP*TEMP3,
     $                    TEN*RTUNFL )
                  ELSE
                     VL = WA1( 1 ) - MAX( HALF*( WA1( N )-WA1( 1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
                  IF( IU.NE.N ) THEN
                     VU = WA1( IU ) + MAX( HALF*
     $                    ( WA1( IU+1 )-WA1( IU ) ), TEN*ULP*TEMP3,
     $                    TEN*RTUNFL )
                  ELSE
                     VU = WA1( N ) + MAX( HALF*( WA1( N )-WA1( 1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
               ELSE
                  VL = ZERO
                  VU = ONE
               END IF
*
               DO 580 I = 1, N
                  D1( I ) = DBLE( A( I, I ) )
  580          CONTINUE
               DO 590 I = 1, N - 1
                  D2( I ) = DBLE( A( I+1, I ) )
  590          CONTINUE
               SRNAMT = 'DSTEVR'
               CALL DSTEVR( 'V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL,
     $                      M2, WA2, Z, LDU, IWORK, WORK, LWORK,
     $                      IWORK(2*N+1), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVR(V,V)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 22 ) = ULPINV
                     RESULT( 23 ) = ULPINV
                     RESULT( 24 ) = ULPINV
                     GO TO 630
                  END IF
               END IF
*
               IF( M2.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( 22 ) = ULPINV
                  RESULT( 23 ) = ULPINV
                  RESULT( 24 ) = ULPINV
                  GO TO 630
               END IF
*
*              Do tests 22 and 23.
*
               DO 600 I = 1, N
                  D3( I ) = DBLE( A( I, I ) )
  600          CONTINUE
               DO 610 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  610          CONTINUE
               CALL DSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
     $                      MAX( 1, M2 ), RESULT( 22 ) )
*
               NTEST = 24
               DO 620 I = 1, N - 1
                  D4( I ) = DBLE( A( I+1, I ) )
  620          CONTINUE
               SRNAMT = 'DSTEVR'
               CALL DSTEVR( 'N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL,
     $                      M3, WA3, Z, LDU, IWORK, WORK, LWORK,
     $                      IWORK(2*N+1), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSTEVR(N,V)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 24 ) = ULPINV
                     GO TO 630
                  END IF
               END IF
*
*              Do test 24.
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 24 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
  630          CONTINUE
*
*
*
            ELSE
*
               DO 640 I = 1, 24
                  RESULT( I ) = ZERO
  640          CONTINUE
               NTEST = 24
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
*              4)      Call DSYEV and DSYEVX.
*
               CALL DLACPY( ' ', N, N, A, LDA, V, LDU )
*
               NTEST = NTEST + 1
               SRNAMT = 'DSYEV'
               CALL DSYEV( 'V', UPLO, N, A, LDU, D1, WORK, LWORK,
     $                     IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSYEV(V,' // UPLO // ')',
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
               CALL DSYT21( 1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 2
               SRNAMT = 'DSYEV_2STAGE'
               CALL DSYEV_2STAGE( 'N', UPLO, N, A, LDU, D3, WORK, LWORK,
     $                     IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSYEV_2STAGE(N,' // UPLO // ')',
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
               DO 650 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  650          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
  660          CONTINUE
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 1
*
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
                  IF( IL.NE.1 ) THEN
                     VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  ELSE IF( N.GT.0 ) THEN
                     VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
                  IF( IU.NE.N ) THEN
                     VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  ELSE IF( N.GT.0 ) THEN
                     VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
               ELSE
                  TEMP3 = ZERO
                  VL = ZERO
                  VU = ONE
               END IF
*
               SRNAMT = 'DSYEVX'
               CALL DSYEVX( 'V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M, WA1, Z, LDU, WORK, LWORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSYEVX(V,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 680
                  END IF
               END IF
*
*              Do tests 28 and 29 (or +54)
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL DSYT21( 1, UPLO, N, 0, A, LDU, D1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               SRNAMT = 'DSYEVX_2STAGE'
               CALL DSYEVX_2STAGE( 'N', 'A', UPLO, N, A, LDU, VL, VU,
     $                      IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK,
     $                      LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSYEVX_2STAGE(N,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 680
                  END IF
               END IF
*
*              Do test 30 (or +54)
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 670 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  670          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
  680          CONTINUE
*
               NTEST = NTEST + 1
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'DSYEVX'
               CALL DSYEVX( 'V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSYEVX(V,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 690
                  END IF
               END IF
*
*              Do tests 31 and 32 (or +54)
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL DSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'DSYEVX_2STAGE'
               CALL DSYEVX_2STAGE( 'N', 'I', UPLO, N, A, LDU, VL, VU,
     $                      IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK,
     $                      LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSYEVX_2STAGE(N,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 690
                  END IF
               END IF
*
*              Do test 33 (or +54)
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, ULP*TEMP3 )
  690          CONTINUE
*
               NTEST = NTEST + 1
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'DSYEVX'
               CALL DSYEVX( 'V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSYEVX(V,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 700
                  END IF
               END IF
*
*              Do tests 34 and 35 (or +54)
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL DSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'DSYEVX_2STAGE'
               CALL DSYEVX_2STAGE( 'N', 'V', UPLO, N, A, LDU, VL, VU,
     $                      IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK,
     $                      LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSYEVX_2STAGE(N,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 700
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 700
               END IF
*
*              Do test 36 (or +54)
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
  700          CONTINUE
*
*              5)      Call DSPEV and DSPEVX.
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
*              Load array WORK with the upper or lower triangular
*              part of the matrix in packed form.
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 720 J = 1, N
                     DO 710 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  710                CONTINUE
  720             CONTINUE
               ELSE
                  INDX = 1
                  DO 740 J = 1, N
                     DO 730 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  730                CONTINUE
  740             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               SRNAMT = 'DSPEV'
               CALL DSPEV( 'V', UPLO, N, WORK, D1, Z, LDU, V, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEV(V,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 800
                  END IF
               END IF
*
*              Do tests 37 and 38 (or +54)
*
               CALL DSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 760 J = 1, N
                     DO 750 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  750                CONTINUE
  760             CONTINUE
               ELSE
                  INDX = 1
                  DO 780 J = 1, N
                     DO 770 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  770                CONTINUE
  780             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               SRNAMT = 'DSPEV'
               CALL DSPEV( 'N', UPLO, N, WORK, D3, Z, LDU, V, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEV(N,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 800
                  END IF
               END IF
*
*              Do test 39 (or +54)
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 790 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  790          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
*              Load array WORK with the upper or lower triangular part
*              of the matrix in packed form.
*
  800          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 820 J = 1, N
                     DO 810 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  810                CONTINUE
  820             CONTINUE
               ELSE
                  INDX = 1
                  DO 840 J = 1, N
                     DO 830 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  830                CONTINUE
  840             CONTINUE
               END IF
*
               NTEST = NTEST + 1
*
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
                  IF( IL.NE.1 ) THEN
                     VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  ELSE IF( N.GT.0 ) THEN
                     VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
                  IF( IU.NE.N ) THEN
                     VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  ELSE IF( N.GT.0 ) THEN
                     VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ),
     $                    TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
               ELSE
                  TEMP3 = ZERO
                  VL = ZERO
                  VU = ONE
               END IF
*
               SRNAMT = 'DSPEVX'
               CALL DSPEVX( 'V', 'A', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M, WA1, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEVX(V,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 900
                  END IF
               END IF
*
*              Do tests 40 and 41 (or +54)
*
               CALL DSYT21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 860 J = 1, N
                     DO 850 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  850                CONTINUE
  860             CONTINUE
               ELSE
                  INDX = 1
                  DO 880 J = 1, N
                     DO 870 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  870                CONTINUE
  880             CONTINUE
               END IF
*
               SRNAMT = 'DSPEVX'
               CALL DSPEVX( 'N', 'A', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEVX(N,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 900
                  END IF
               END IF
*
*              Do test 42 (or +54)
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 890 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  890          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
  900          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 920 J = 1, N
                     DO 910 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  910                CONTINUE
  920             CONTINUE
               ELSE
                  INDX = 1
                  DO 940 J = 1, N
                     DO 930 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  930                CONTINUE
  940             CONTINUE
               END IF
*
               NTEST = NTEST + 1
*
               SRNAMT = 'DSPEVX'
               CALL DSPEVX( 'V', 'I', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEVX(V,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 990
                  END IF
               END IF
*
*              Do tests 43 and 44 (or +54)
*
               CALL DSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 960 J = 1, N
                     DO 950 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  950                CONTINUE
  960             CONTINUE
               ELSE
                  INDX = 1
                  DO 980 J = 1, N
                     DO 970 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  970                CONTINUE
  980             CONTINUE
               END IF
*
               SRNAMT = 'DSPEVX'
               CALL DSPEVX( 'N', 'I', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M3, WA3, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEVX(N,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 990
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 990
               END IF
*
*              Do test 45 (or +54)
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
  990          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 1010 J = 1, N
                     DO 1000 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1000                CONTINUE
 1010             CONTINUE
               ELSE
                  INDX = 1
                  DO 1030 J = 1, N
                     DO 1020 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1020                CONTINUE
 1030             CONTINUE
               END IF
*
               NTEST = NTEST + 1
*
               SRNAMT = 'DSPEVX'
               CALL DSPEVX( 'V', 'V', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEVX(V,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1080
                  END IF
               END IF
*
*              Do tests 46 and 47 (or +54)
*
               CALL DSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 1050 J = 1, N
                     DO 1040 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1040                CONTINUE
 1050             CONTINUE
               ELSE
                  INDX = 1
                  DO 1070 J = 1, N
                     DO 1060 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1060                CONTINUE
 1070             CONTINUE
               END IF
*
               SRNAMT = 'DSPEVX'
               CALL DSPEVX( 'N', 'V', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M3, WA3, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEVX(N,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1080
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 1080
               END IF
*
*              Do test 48 (or +54)
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
 1080          CONTINUE
*
*              6)      Call DSBEV and DSBEVX.
*
               IF( JTYPE.LE.7 ) THEN
                  KD = 1
               ELSE IF( JTYPE.GE.8 .AND. JTYPE.LE.15 ) THEN
                  KD = MAX( N-1, 0 )
               ELSE
                  KD = IHBW
               END IF
*
*              Load array V with the upper or lower triangular part
*              of the matrix in band form.
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1100 J = 1, N
                     DO 1090 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1090                CONTINUE
 1100             CONTINUE
               ELSE
                  DO 1120 J = 1, N
                     DO 1110 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1110                CONTINUE
 1120             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               SRNAMT = 'DSBEV'
               CALL DSBEV( 'V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK,
     $                     IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSBEV(V,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1180
                  END IF
               END IF
*
*              Do tests 49 and 50 (or ... )
*
               CALL DSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1140 J = 1, N
                     DO 1130 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1130                CONTINUE
 1140             CONTINUE
               ELSE
                  DO 1160 J = 1, N
                     DO 1150 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1150                CONTINUE
 1160             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               SRNAMT = 'DSBEV_2STAGE'
               CALL DSBEV_2STAGE( 'N', UPLO, N, KD, V, LDU, D3, Z, LDU,
     $                     WORK, LWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSBEV_2STAGE(N,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1180
                  END IF
               END IF
*
*              Do test 51 (or +54)
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1170 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1170          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
*              Load array V with the upper or lower triangular part
*              of the matrix in band form.
*
 1180          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  DO 1200 J = 1, N
                     DO 1190 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1190                CONTINUE
 1200             CONTINUE
               ELSE
                  DO 1220 J = 1, N
                     DO 1210 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1210                CONTINUE
 1220             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               SRNAMT = 'DSBEVX'
               CALL DSBEVX( 'V', 'A', UPLO, N, KD, V, LDU, U, LDU, VL,
     $                      VU, IL, IU, ABSTOL, M, WA2, Z, LDU, WORK,
     $                      IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSBEVX(V,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1280
                  END IF
               END IF
*
*              Do tests 52 and 53 (or +54)
*
               CALL DSYT21( 1, UPLO, N, 0, A, LDU, WA2, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1240 J = 1, N
                     DO 1230 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1230                CONTINUE
 1240             CONTINUE
               ELSE
                  DO 1260 J = 1, N
                     DO 1250 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1250                CONTINUE
 1260             CONTINUE
               END IF
*
               SRNAMT = 'DSBEVX_2STAGE'
               CALL DSBEVX_2STAGE( 'N', 'A', UPLO, N, KD, V, LDU,
     $                      U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3,
     $                      Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ),
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSBEVX_2STAGE(N,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1280
                  END IF
               END IF
*
*              Do test 54 (or +54)
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1270 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA2( J ) ), ABS( WA3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA2( J )-WA3( J ) ) )
 1270          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
 1280          CONTINUE
               NTEST = NTEST + 1
               IF( IUPLO.EQ.1 ) THEN
                  DO 1300 J = 1, N
                     DO 1290 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1290                CONTINUE
 1300             CONTINUE
               ELSE
                  DO 1320 J = 1, N
                     DO 1310 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1310                CONTINUE
 1320             CONTINUE
               END IF
*
               SRNAMT = 'DSBEVX'
               CALL DSBEVX( 'V', 'I', UPLO, N, KD, V, LDU, U, LDU, VL,
     $                      VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK,
     $                      IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSBEVX(V,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1370
                  END IF
               END IF
*
*              Do tests 55 and 56 (or +54)
*
               CALL DSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1340 J = 1, N
                     DO 1330 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1330                CONTINUE
 1340             CONTINUE
               ELSE
                  DO 1360 J = 1, N
                     DO 1350 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1350                CONTINUE
 1360             CONTINUE
               END IF
*
               SRNAMT = 'DSBEVX_2STAGE'
               CALL DSBEVX_2STAGE( 'N', 'I', UPLO, N, KD, V, LDU,
     $                      U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3,
     $                      Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ),
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSBEVX_2STAGE(N,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1370
                  END IF
               END IF
*
*              Do test 57 (or +54)
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
 1370          CONTINUE
               NTEST = NTEST + 1
               IF( IUPLO.EQ.1 ) THEN
                  DO 1390 J = 1, N
                     DO 1380 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1380                CONTINUE
 1390             CONTINUE
               ELSE
                  DO 1410 J = 1, N
                     DO 1400 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1400                CONTINUE
 1410             CONTINUE
               END IF
*
               SRNAMT = 'DSBEVX'
               CALL DSBEVX( 'V', 'V', UPLO, N, KD, V, LDU, U, LDU, VL,
     $                      VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK,
     $                      IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSBEVX(V,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1460
                  END IF
               END IF
*
*              Do tests 58 and 59 (or +54)
*
               CALL DSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1430 J = 1, N
                     DO 1420 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1420                CONTINUE
 1430             CONTINUE
               ELSE
                  DO 1450 J = 1, N
                     DO 1440 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1440                CONTINUE
 1450             CONTINUE
               END IF
*
               SRNAMT = 'DSBEVX_2STAGE'
               CALL DSBEVX_2STAGE( 'N', 'V', UPLO, N, KD, V, LDU,
     $                      U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3,
     $                      Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ),
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSBEVX_2STAGE(N,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1460
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 1460
               END IF
*
*              Do test 60 (or +54)
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
 1460          CONTINUE
*
*              7)      Call DSYEVD
*
               CALL DLACPY( ' ', N, N, A, LDA, V, LDU )
*
               NTEST = NTEST + 1
               SRNAMT = 'DSYEVD'
               CALL DSYEVD( 'V', UPLO, N, A, LDU, D1, WORK, LWEDC,
     $                      IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSYEVD(V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1480
                  END IF
               END IF
*
*              Do tests 61 and 62 (or +54)
*
               CALL DSYT21( 1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 2
               SRNAMT = 'DSYEVD_2STAGE'
               CALL DSYEVD_2STAGE( 'N', UPLO, N, A, LDU, D3, WORK, 
     $                              LWORK, IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSYEVD_2STAGE(N,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1480
                  END IF
               END IF
*
*              Do test 63 (or +54)
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1470 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1470          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
 1480          CONTINUE
*
*              8)      Call DSPEVD.
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
*              Load array WORK with the upper or lower triangular
*              part of the matrix in packed form.
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 1500 J = 1, N
                     DO 1490 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1490                CONTINUE
 1500             CONTINUE
               ELSE
                  INDX = 1
                  DO 1520 J = 1, N
                     DO 1510 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1510                CONTINUE
 1520             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               SRNAMT = 'DSPEVD'
               CALL DSPEVD( 'V', UPLO, N, WORK, D1, Z, LDU,
     $                      WORK( INDX ), LWEDC-INDX+1, IWORK, LIWEDC,
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEVD(V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1580
                  END IF
               END IF
*
*              Do tests 64 and 65 (or +54)
*
               CALL DSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 1540 J = 1, N
                     DO 1530 I = 1, J
*
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1530                CONTINUE
 1540             CONTINUE
               ELSE
                  INDX = 1
                  DO 1560 J = 1, N
                     DO 1550 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1550                CONTINUE
 1560             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               SRNAMT = 'DSPEVD'
               CALL DSPEVD( 'N', UPLO, N, WORK, D3, Z, LDU,
     $                      WORK( INDX ), LWEDC-INDX+1, IWORK, LIWEDC,
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSPEVD(N,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1580
                  END IF
               END IF
*
*              Do test 66 (or +54)
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1570 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1570          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
 1580          CONTINUE
*
*              9)      Call DSBEVD.
*
               IF( JTYPE.LE.7 ) THEN
                  KD = 1
               ELSE IF( JTYPE.GE.8 .AND. JTYPE.LE.15 ) THEN
                  KD = MAX( N-1, 0 )
               ELSE
                  KD = IHBW
               END IF
*
*              Load array V with the upper or lower triangular part
*              of the matrix in band form.
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1600 J = 1, N
                     DO 1590 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1590                CONTINUE
 1600             CONTINUE
               ELSE
                  DO 1620 J = 1, N
                     DO 1610 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1610                CONTINUE
 1620             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               SRNAMT = 'DSBEVD'
               CALL DSBEVD( 'V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK,
     $                      LWEDC, IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSBEVD(V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1680
                  END IF
               END IF
*
*              Do tests 67 and 68 (or +54)
*
               CALL DSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1640 J = 1, N
                     DO 1630 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1630                CONTINUE
 1640             CONTINUE
               ELSE
                  DO 1660 J = 1, N
                     DO 1650 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1650                CONTINUE
 1660             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               SRNAMT = 'DSBEVD_2STAGE'
               CALL DSBEVD_2STAGE( 'N', UPLO, N, KD, V, LDU, D3, Z, LDU,
     $                             WORK, LWORK, IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSBEVD_2STAGE(N,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1680
                  END IF
               END IF
*
*              Do test 69 (or +54)
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1670 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1670          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
 1680          CONTINUE
*
*
               CALL DLACPY( ' ', N, N, A, LDA, V, LDU )
               NTEST = NTEST + 1
               SRNAMT = 'DSYEVR'
               CALL DSYEVR( 'V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M, WA1, Z, LDU, IWORK, WORK, LWORK,
     $                      IWORK(2*N+1), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSYEVR(V,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1700
                  END IF
               END IF
*
*              Do tests 70 and 71 (or ... )
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL DSYT21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               SRNAMT = 'DSYEVR_2STAGE'
               CALL DSYEVR_2STAGE( 'N', 'A', UPLO, N, A, LDU, VL, VU,
     $                      IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK,
     $                      WORK, LWORK, IWORK(2*N+1), LIWORK-2*N,
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSYEVR_2STAGE(N,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1700
                  END IF
               END IF
*
*              Do test 72 (or ... )
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1690 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
 1690          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
 1700          CONTINUE
*
               NTEST = NTEST + 1
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'DSYEVR'
               CALL DSYEVR( 'V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK,
     $                      IWORK(2*N+1), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSYEVR(V,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1710
                  END IF
               END IF
*
*              Do tests 73 and 74 (or +54)
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL DSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'DSYEVR_2STAGE'
               CALL DSYEVR_2STAGE( 'N', 'I', UPLO, N, A, LDU, VL, VU,
     $                      IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK,
     $                      WORK, LWORK, IWORK(2*N+1), LIWORK-2*N,
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSYEVR_2STAGE(N,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1710
                  END IF
               END IF
*
*              Do test 75 (or +54)
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, ULP*TEMP3 )
 1710          CONTINUE
*
               NTEST = NTEST + 1
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'DSYEVR'
               CALL DSYEVR( 'V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK,
     $                      IWORK(2*N+1), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSYEVR(V,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1720
                  END IF
               END IF
*
*              Do tests 76 and 77 (or +54)
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL DSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'DSYEVR_2STAGE'
               CALL DSYEVR_2STAGE( 'N', 'V', UPLO, N, A, LDU, VL, VU,
     $                      IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK,
     $                      WORK, LWORK, IWORK(2*N+1), LIWORK-2*N,
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
     $               'DSYEVR_2STAGE(N,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1720
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 1720
               END IF
*
*              Do test 78 (or +54)
*
               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
               CALL DLACPY( ' ', N, N, V, LDU, A, LDA )
*
 1720       CONTINUE
*
*           End of Loop -- Check for RESULT(j) > THRESH
*
            NTESTT = NTESTT + NTEST
*
            CALL DLAFTS( 'DST', N, N, JTYPE, NTEST, RESULT, IOLDSD,
     $                   THRESH, NOUNIT, NERRS )
*
 1730    CONTINUE
 1740 CONTINUE
*
*     Summary
*
      CALL ALASVM( 'DST', NOUNIT, NERRS, NTESTT, 0 )
*
 9999 FORMAT( ' DDRVST2STG: ', A, ' returned INFO=', I6, '.', / 9X,
     $    'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of DDRVST2STG
*
      END
