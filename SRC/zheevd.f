*> \brief <b> ZHEEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZHEEVD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevd.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevd.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevd.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
*                          LRWORK, IWORK, LIWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   RWORK( * ), W( * )
*       COMPLEX*16         A( LDA, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a
*> complex Hermitian matrix A.  If eigenvectors are desired, it uses a
*> divide and conquer algorithm.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBZ
*> \verbatim
*>          JOBZ is CHARACTER*1
*>          = 'N':  Compute eigenvalues only;
*>          = 'V':  Compute eigenvalues and eigenvectors.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, N)
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
*>          leading N-by-N upper triangular part of A contains the
*>          upper triangular part of the matrix A.  If UPLO = 'L',
*>          the leading N-by-N lower triangular part of A contains
*>          the lower triangular part of the matrix A.
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*>          orthonormal eigenvectors of the matrix A.
*>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*>          or the upper triangle (if UPLO='U') of A, including the
*>          diagonal, is destroyed.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is DOUBLE PRECISION array, dimension (N)
*>          If INFO = 0, the eigenvalues in ascending order.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of the array WORK.
*>          If N <= 1,                LWORK must be at least 1.
*>          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.
*>          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal sizes of the WORK, RWORK and
*>          IWORK arrays, returns these values as the first entries of
*>          the WORK, RWORK and IWORK arrays, and no error message
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
*>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
*> \endverbatim
*>
*> \param[in] LRWORK
*> \verbatim
*>          LRWORK is INTEGER
*>          The dimension of the array RWORK.
*>          If N <= 1,                LRWORK must be at least 1.
*>          If JOBZ  = 'N' and N > 1, LRWORK must be at least N.
*>          If JOBZ  = 'V' and N > 1, LRWORK must be at least
*>                         1 + 5*N + 2*N**2.
*>
*>          If LRWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal sizes of the WORK, RWORK
*>          and IWORK arrays, returns these values as the first entries
*>          of the WORK, RWORK and IWORK arrays, and no error message
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*> \endverbatim
*>
*> \param[in] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          The dimension of the array IWORK.
*>          If N <= 1,                LIWORK must be at least 1.
*>          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
*>          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
*>
*>          If LIWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal sizes of the WORK, RWORK
*>          and IWORK arrays, returns these values as the first entries
*>          of the WORK, RWORK and IWORK arrays, and no error message
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
*>                to converge; i off-diagonal elements of an intermediate
*>                tridiagonal form did not converge to zero;
*>                if INFO = i and JOBZ = 'V', then the algorithm failed
*>                to compute an eigenvalue while working on the submatrix
*>                lying in rows and columns INFO/(N+1) through
*>                mod(INFO,N+1).
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
*> \ingroup heevd
*
*> \par Further Details:
*  =====================
*>
*>  Modified description of INFO. Sven, 16 Feb 05.
*
*> \par Contributors:
*  ==================
*>
*> Jeff Rutter, Computer Science Division, University of California
*> at Berkeley, USA
*>
*  =====================================================================
      SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK,
     $                   RWORK,
     $                   LRWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDRWK, INDTAU, INDWK2,
     $                   INDWRK, ISCALE, LIOPT, LIWMIN, LLRWK, LLWORK,
     $                   LLWRK2, LOPT, LROPT, LRWMIN, LWMIN
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANHE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTERF, XERBLA, ZHETRD, ZLACPY,
     $                   ZLASCL,
     $                   ZSTEDC, ZUNMTR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.LE.1 ) THEN
            LWMIN = 1
            LRWMIN = 1
            LIWMIN = 1
            LOPT = LWMIN
            LROPT = LRWMIN
            LIOPT = LIWMIN
         ELSE
            IF( WANTZ ) THEN
               LWMIN = 2*N + N*N
               LRWMIN = 1 + 5*N + 2*N**2
               LIWMIN = 3 + 5*N
            ELSE
               LWMIN = N + 1
               LRWMIN = N
               LIWMIN = 1
            END IF
            LOPT = MAX( LWMIN, N +
     $                  N*ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1,
     $                            -1 ) )
            LROPT = LRWMIN
            LIOPT = LIWMIN
         END IF
         WORK( 1 ) = LOPT
         RWORK( 1 ) = REAL( LROPT )
         IWORK( 1 ) = LIOPT
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -8
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -10
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = DBLE( A( 1, 1 ) )
         IF( WANTZ )
     $      A( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL ZLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
*
*     Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = 1
      INDWRK = INDTAU + N
      INDRWK = INDE + N
      INDWK2 = INDWRK + N*N
      LLWORK = LWORK - INDWRK + 1
      LLWRK2 = LWORK - INDWK2 + 1
      LLRWK = LRWORK - INDRWK + 1
      CALL ZHETRD( UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, first call
*     ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
*     tridiagonal matrix, then call ZUNMTR to multiply it to the
*     Householder transformations represented as Householder vectors in
*     A.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL ZSTEDC( 'I', N, W, RWORK( INDE ), WORK( INDWRK ), N,
     $                WORK( INDWK2 ), LLWRK2, RWORK( INDRWK ), LLRWK,
     $                IWORK, LIWORK, INFO )
         CALL ZUNMTR( 'L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ),
     $                WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO )
         CALL ZLACPY( 'A', N, N, WORK( INDWRK ), N, A, LDA )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
      WORK( 1 ) = LOPT
      RWORK( 1 ) = REAL( LROPT )
      IWORK( 1 ) = LIOPT
*
      RETURN
*
*     End of ZHEEVD
*
      END
