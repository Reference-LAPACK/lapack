*> \brief <b> SSYEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b>
*
*  @generated from dsyev_2stage.f, fortran d -> s, Sat Nov  5 23:55:51 2016
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SSYEV_2STAGE + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyev_2stage.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyev_2stage.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyev_2stage.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SSYEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK,
*                                INFO )
*
*       IMPLICIT NONE
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), W( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SSYEV_2STAGE computes all eigenvalues and, optionally, eigenvectors of a
*> real symmetric matrix A using the 2stage technique for
*> the reduction to tridiagonal.
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
*>                  Not available in this release.
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
*>          A is REAL array, dimension (LDA, N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the
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
*>          W is REAL array, dimension (N)
*>          If INFO = 0, the eigenvalues in ascending order.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension LWORK
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of the array WORK. LWORK >= 1, when N <= 1;
*>          otherwise
*>          If JOBZ = 'N' and N > 1, LWORK must be queried.
*>                                   LWORK = MAX(1, dimension) where
*>                                   dimension = max(stage1,stage2) + (KD+1)*N + 2*N
*>                                             = N*KD + N*max(KD+1,FACTOPTNB)
*>                                               + max(2*KD*KD, KD*NTHREADS)
*>                                               + (KD+1)*N + 2*N
*>                                   where KD is the blocking size of the reduction,
*>                                   FACTOPTNB is the blocking used by the QR or LQ
*>                                   algorithm, usually FACTOPTNB=128 is a good choice
*>                                   NTHREADS is the number of threads used when
*>                                   openMP compilation is enabled, otherwise =1.
*>          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the algorithm failed to converge; i
*>                off-diagonal elements of an intermediate tridiagonal
*>                form did not converge to zero.
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
*> \ingroup heev_2stage
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  All details about the 2stage techniques are available in:
*>
*>  Azzam Haidar, Hatem Ltaief, and Jack Dongarra.
*>  Parallel reduction to condensed forms for symmetric eigenvalue problems
*>  using aggregated fine-grained and memory-aware kernels. In Proceedings
*>  of 2011 International Conference for High Performance Computing,
*>  Networking, Storage and Analysis (SC '11), New York, NY, USA,
*>  Article 8 , 11 pages.
*>  http://doi.acm.org/10.1145/2063384.2063394
*>
*>  A. Haidar, J. Kurzak, P. Luszczek, 2013.
*>  An improved parallel singular value algorithm and its implementation
*>  for multicore hardware, In Proceedings of 2013 International Conference
*>  for High Performance Computing, Networking, Storage and Analysis (SC '13).
*>  Denver, Colorado, USA, 2013.
*>  Article 90, 12 pages.
*>  http://doi.acm.org/10.1145/2503210.2503292
*>
*>  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra.
*>  A novel hybrid CPU-GPU generalized eigensolver for electronic structure
*>  calculations based on fine-grained memory aware tasks.
*>  International Journal of High Performance Computing Applications.
*>  Volume 28 Issue 2, Pages 196-209, May 2014.
*>  http://hpc.sagepub.com/content/28/2/196
*>
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE SSYEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK,
     $                         INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), W( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV2STAGE
      REAL               SLAMCH, SLANSY, SROUNDUP_LWORK
      EXTERNAL           LSAME, SLAMCH, SLANSY, ILAENV2STAGE,
     $                   SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASCL, SORGTR, SSCAL, SSTEQR,
     $                   SSTERF,
     $                   XERBLA, SSYTRD_2STAGE
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
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( LSAME( JOBZ, 'N' ) ) ) THEN
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
         KD    = ILAENV2STAGE( 1, 'SSYTRD_2STAGE', JOBZ, N, -1, -1,
     $                         -1 )
         IB    = ILAENV2STAGE( 2, 'SSYTRD_2STAGE', JOBZ, N, KD, -1,
     $                         -1 )
         LHTRD = ILAENV2STAGE( 3, 'SSYTRD_2STAGE', JOBZ, N, KD, IB,
     $                         -1 )
         LWTRD = ILAENV2STAGE( 4, 'SSYTRD_2STAGE', JOBZ, N, KD, IB,
     $                         -1 )
         LWMIN = 2*N + LHTRD + LWTRD
         WORK( 1 )  = REAL( LWMIN )
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY )
     $      INFO = -8
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYEV_2STAGE ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 2
         IF( WANTZ )
     $      A( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS    = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = SLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL SLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
*
*     Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form.
*
      INDE    = 1
      INDTAU  = INDE + N
      INDHOUS = INDTAU + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1
*
      CALL SSYTRD_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK( INDE ),
     $                    WORK( INDTAU ), WORK( INDHOUS ), LHTRD,
     $                    WORK( INDWRK ), LLWORK, IINFO )
*
*     For eigenvalues only, call SSTERF.  For eigenvectors, first call
*     SORGTR to generate the orthogonal matrix, then call SSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, WORK( INDE ), INFO )
      ELSE
*        Not available in this release, and argument checking should not
*        let it getting here
         RETURN
         CALL SORGTR( UPLO, N, A, LDA, WORK( INDTAU ),
     $                WORK( INDWRK ),
     $                LLWORK, IINFO )
         CALL SSTEQR( JOBZ, N, W, WORK( INDE ), A, LDA,
     $                WORK( INDTAU ),
     $                INFO )
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
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
*
      RETURN
*
*     End of SSYEV_2STAGE
*
      END
