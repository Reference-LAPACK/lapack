      SUBROUTINE CSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL,
     $                   M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $                   LIWORK, INFO )
*
*  -- LAPACK computational routine (instru to count ops, version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N
      REAL               ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            ISUPPZ( * ), IWORK( * )
      REAL               D( * ), E( * ), W( * ), WORK( * )
      COMPLEX            Z( LDZ, * )
*     ..
*     Common block to return operation count ..
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      REAL               ITCNT, OPS
*     ..
*
*  Purpose
*  =======
*
*  CSTEGR computes eigenvalues by the dqds algorithm, while
*  orthogonal eigenvectors are computed from various "good" L D L^T
*  representations (also known as Relatively Robust Representations).
*  Gram-Schmidt orthogonalization is avoided as far as possible. More
*  specifically, the various steps of the algorithm are as follows.
*  For the i-th unreduced block of T,
*     (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
*         is a relatively robust representation,
*     (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
*         relative accuracy by the dqds algorithm,
*     (c) If there is a cluster of close eigenvalues, "choose" sigma_i
*         close to the cluster, and go to step (a),
*     (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
*         compute the corresponding eigenvector by forming a
*         rank-revealing twisted factorization.
*  The desired accuracy of the output can be specified by the input
*  parameter ABSTOL.
*
*  For more details, see "A new O(n^2) algorithm for the symmetric
*  tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
*  Computer Science Division Technical Report No. UCB/CSD-97-971,
*  UC Berkeley, May 1997.
*
*  Note 1 : Currently CSTEGR is only set up to find ALL the n
*  eigenvalues and eigenvectors of T in O(n^2) time
*  Note 2 : Currently the routine CSTEIN is called when an appropriate
*  sigma_i cannot be chosen in step (c) above. CSTEIN invokes modified
*  Gram-Schmidt when eigenvalues are close.
*  Note 3 : CSTEGR works only on machines which follow ieee-754
*  floating-point standard in their handling of infinities and NaNs.
*  Normal execution of CSTEGR may create NaNs and infinities and hence
*  may abort due to a floating point exception in environments which
*  do not conform to the ieee standard.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
********** Only RANGE = 'A' is currently supported *********************
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix
*          T. On exit, D is overwritten.
*
*  E       (input/output) REAL array, dimension (N)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix T in elements 1 to N-1 of E; E(N) need not be set.
*          On exit, E is overwritten.
*
*  VL      (input) REAL
*  VU      (input) REAL
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) REAL
*          The absolute error tolerance for the
*          eigenvalues/eigenvectors. IF JOBZ = 'V', the eigenvalues and
*          eigenvectors output have residual norms bounded by ABSTOL,
*          and the dot products between different eigenvectors are
*          bounded by ABSTOL. If ABSTOL is less than N*EPS*|T|, then
*          N*EPS*|T| will be used in its place, where EPS is the
*          machine precision and |T| is the 1-norm of the tridiagonal
*          matrix. The eigenvalues are computed to an accuracy of
*          EPS*|T| irrespective of ABSTOL. If high relative accuracy
*          is important, set ABSTOL to DLAMCH( 'Safe minimum' ).
*          See Barlow and Demmel "Computing Accurate Eigensystems of
*          Scaled Diagonally Dominant Matrices", LAPACK Working Note #7
*          for a discussion of which matrices define their eigenvalues
*          to high relative accuracy.
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (output) REAL array, dimension (N)
*          The first M elements contain the selected eigenvalues in
*          ascending order.
*
*  Z       (output) COMPLEX array, dimension (LDZ, max(1,M) )
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix T
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,M) )
*          The support of the eigenvectors in Z, i.e., the indices
*          indicating the nonzero elements in Z. The i-th eigenvector
*          is nonzero only in elements ISUPPZ( 2*i-1 ) through
*          ISUPPZ( 2*i ).
*
*  WORK    (workspace/output) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal
*          (and minimal) LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,18*N)
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.  LIWORK >= max(1,10*N)
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = 1, internal error in SLARRE,
*                if INFO = 2, internal error in CLARRV.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Inderjit Dhillon, IBM Almaden, USA
*     Osni Marques, LBNL/NERSC, USA
*     Ken Stanley, Computer Science Division, University of
*       California at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LQUERY, VALEIG, WANTZ
      INTEGER            I, IBEGIN, IEND, IINDBL, IINDWK, IINFO, IINSPL, 
     $                   INDGRS, INDWOF, INDWRK, ITMP, J, JJ, LIWMIN, 
     $                   LWMIN, NSPLIT
      REAL               BIGNUM, EPS, RMAX, RMIN, SAFMIN, SCALE, SMLNUM,
     $                   THRESH, TMP, TNRM, TOL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANST
      EXTERNAL           LSAME, SLAMCH, SLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARRV, CLASET, CSWAP, SLARRE, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )
      LWMIN = 18*N
      LIWMIN = 10*N
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
*
*     The following two lines need to be removed once the 
*     RANGE = 'V' and RANGE = 'I' options are provided.
*
      ELSE IF( VALEIG .OR. INDEIG ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
         INFO = -7
      ELSE IF( INDEIG .AND. IL.LT.1 ) THEN
         INFO = -8
*     The following change should be made in DSTEVX also, otherwise
*     IL can be specified as N+1 and IU as N.
*     ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) THEN
      ELSE IF( INDEIG .AND. ( IU.LT.IL .OR. IU.GT.N ) ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -17
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -19
      END IF
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSTEGR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = D( 1 )
         ELSE
            IF( VL.LT.D( 1 ) .AND. VU.GE.D( 1 ) ) THEN
               M = 1
               W( 1 ) = D( 1 )
            END IF
         END IF
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      OPS = OPS + REAL( 7 )
      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Scale matrix to allowable range, if necessary.
*
      SCALE = ONE
      TNRM = SLANST( 'M', N, D, E )
      IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
         OPS = OPS + REAL( 1 )
         SCALE = RMIN / TNRM
      ELSE IF( TNRM.GT.RMAX ) THEN
         OPS = OPS + REAL( 1 )
         SCALE = RMAX / TNRM
      END IF
      IF( SCALE.NE.ONE ) THEN
         OPS = OPS + REAL( 2*N )
         CALL SSCAL( N, SCALE, D, 1 )
         CALL SSCAL( N-1, SCALE, E, 1 )
         TNRM = TNRM*SCALE
      END IF
      INDGRS = 1
      INDWOF = 2*N + 1
      INDWRK = 3*N + 1
*
      IINSPL = 1
      IINDBL = N + 1
      IINDWK = 2*N + 1
*
      CALL CLASET( 'Full', N, N, CZERO, CZERO, Z, LDZ )
*
*     Compute the desired eigenvalues of the tridiagonal after splitting
*     into smaller subblocks if the corresponding of-diagonal elements
*     are small
*
      OPS = OPS + REAL( 1 )
      THRESH = EPS*TNRM
      CALL SLARRE( N, D, E, THRESH, NSPLIT, IWORK( IINSPL ), M, W,
     $             WORK( INDWOF ), WORK( INDGRS ), WORK( INDWRK ),
     $             IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
*
      IF( WANTZ ) THEN
*
*        Compute the desired eigenvectors corresponding to the computed
*        eigenvalues
*
         OPS = OPS + REAL( 1 )
         TOL = MAX( ABSTOL, REAL( N )*THRESH )
         IBEGIN = 1
         DO 20 I = 1, NSPLIT
            IEND = IWORK( IINSPL+I-1 )
            DO 10 J = IBEGIN, IEND
               IWORK( IINDBL+J-1 ) = I
   10       CONTINUE
            IBEGIN = IEND + 1
   20    CONTINUE
*
         CALL CLARRV( N, D, E, IWORK( IINSPL ), M, W, IWORK( IINDBL ),
     $                WORK( INDGRS ), TOL, Z, LDZ, ISUPPZ,
     $                WORK( INDWRK ), IWORK( IINDWK ), IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 2
            RETURN
         END IF
*
      END IF
*
      IBEGIN = 1
      DO 40 I = 1, NSPLIT
         IEND = IWORK( IINSPL+I-1 )
         DO 30 J = IBEGIN, IEND
            OPS = OPS + REAL( 1 )
            W( J ) = W( J ) + WORK( INDWOF+I-1 )
   30    CONTINUE
         IBEGIN = IEND + 1
   40 CONTINUE
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( SCALE.NE.ONE ) THEN
         CALL SSCAL( M, ONE / SCALE, W, 1 )
      END IF
*
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.
*
      IF( NSPLIT.GT.1 ) THEN
         DO 60 J = 1, M - 1
            I = 0
            TMP = W( J )
            DO 50 JJ = J + 1, M
               IF( W( JJ ).LT.TMP ) THEN
                  I = JJ
                  TMP = W( JJ )
               END IF
   50       CONTINUE
            IF( I.NE.0 ) THEN
               W( I ) = W( J )
               W( J ) = TMP
               IF( WANTZ ) THEN
                  CALL CSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
                  ITMP = ISUPPZ( 2*I-1 )
                  ISUPPZ( 2*I-1 ) = ISUPPZ( 2*J-1 )
                  ISUPPZ( 2*J-1 ) = ITMP
                  ITMP = ISUPPZ( 2*I )
                  ISUPPZ( 2*I ) = ISUPPZ( 2*J )
                  ISUPPZ( 2*J ) = ITMP
               END IF
            END IF
   60    CONTINUE
      END IF
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of CSTEGR
*
      END
