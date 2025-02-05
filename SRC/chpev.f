*> \brief <b> CHPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CHPEV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpev.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpev.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpev.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK,
*                         INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, LDZ, N
*       ..
*       .. Array Arguments ..
*       REAL               RWORK( * ), W( * )
*       COMPLEX            AP( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CHPEV computes all the eigenvalues and, optionally, eigenvectors of a
*> complex Hermitian matrix in packed storage.
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
*> \param[in,out] AP
*> \verbatim
*>          AP is COMPLEX array, dimension (N*(N+1)/2)
*>          On entry, the upper or lower triangle of the Hermitian matrix
*>          A, packed columnwise in a linear array.  The j-th column of A
*>          is stored in the array AP as follows:
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
*>
*>          On exit, AP is overwritten by values generated during the
*>          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
*>          and first superdiagonal of the tridiagonal matrix T overwrite
*>          the corresponding elements of A, and if UPLO = 'L', the
*>          diagonal and first subdiagonal of T overwrite the
*>          corresponding elements of A.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is REAL array, dimension (N)
*>          If INFO = 0, the eigenvalues in ascending order.
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is COMPLEX array, dimension (LDZ, N)
*>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*>          eigenvectors of the matrix A, with the i-th column of Z
*>          holding the eigenvector associated with W(i).
*>          If JOBZ = 'N', then Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1, and if
*>          JOBZ = 'V', LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (max(1, 2*N-1))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (max(1, 3*N-2))
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
*> \ingroup hpev
*
*  =====================================================================
      SUBROUTINE CHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK,
     $                  INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * ), W( * )
      COMPLEX            AP( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTZ
      INTEGER            IINFO, IMAX, INDE, INDRWK, INDTAU, INDWRK,
     $                   ISCALE
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANHP, SLAMCH
      EXTERNAL           LSAME, CLANHP, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHPTRD, CSSCAL, CSTEQR, CUPGTR, SSCAL,
     $                   SSTERF,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( UPLO, 'L' ) .OR.
     $         LSAME( UPLO, 'U' ) ) )
     $          THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHPEV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = REAL( AP( 1 ) )
         RWORK( 1 ) = 1
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = CLANHP( 'M', UPLO, N, AP, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL CSSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
      END IF
*
*     Call CHPTRD to reduce Hermitian packed matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = 1
      CALL CHPTRD( UPLO, N, AP, W, RWORK( INDE ), WORK( INDTAU ),
     $             IINFO )
*
*     For eigenvalues only, call SSTERF.  For eigenvectors, first call
*     CUPGTR to generate the orthogonal matrix, then call CSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         INDWRK = INDTAU + N
         CALL CUPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ,
     $                WORK( INDWRK ), IINFO )
         INDRWK = INDE + N
         CALL CSTEQR( JOBZ, N, W, RWORK( INDE ), Z, LDZ,
     $                RWORK( INDRWK ), INFO )
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
      RETURN
*
*     End of CHPEV
*
      END
