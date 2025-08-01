*> \brief <b> SKTEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SKTEV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sktev.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sktev.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sktev.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ
*       INTEGER            INFO, LDZ, N
*       ..
*       .. Array Arguments ..
*       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKTEV computes all eigenvalues and, optionally, eigenvectors of a
*> real skew-symmetric tridiagonal matrix A.
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
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix.  N >= 0.
*> \endverbatim
*>
*> \param[out] D
*> \verbatim
*>          D is REAL array, dimension (N)
*>          If INFO = 0, the (N-1) lower subdiagonal elements of the
*>          block diagonal matrix at front, and zero at last.
*>          The matrix consists of 2-by-2 skew-symmetric blocks, and zeros.
*>          The values in D, which represent blocks, are always
*>          positive, and sorted in descending order.
*>          The eigenvalues of each blocks can be evaluated directly.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is REAL array, dimension (N-1)
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal
*>          matrix A, stored in elements 1 to N-1 of E.
*>          On exit, the contents of E are destroyed.
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is REAL array, dimension (LDZ, N)
*>          If JOBZ = 'V', then if INFO = 0, Z is the orthogonal matrix
*>          transforming the skew-symmetric tridiagonal matrix to the
*>          block diagonal matrix. The eigenvectors of corresponding matrix
*>          can be evaluated directly.
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
*>          WORK is REAL array, dimension (max(1,2*N-4))
*>          If JOBZ = 'N', WORK is not referenced.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the algorithm failed to converge; i
*>                off-diagonal elements of E did not converge to zero.
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
*> \ingroup stev
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Derived from subroutine sstev.
*>
*> \endverbatim
*
*> \par Contributors:
*  ==================
*>
*>    Shuo Zheng, China, Jul 2025 \n
*>
*  =====================================================================
      SUBROUTINE SKTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
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
      INTEGER            IMAX, ISCALE
      REAL               BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM,
     $                   TNRM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANKT
      EXTERNAL           LSAME, SLAMCH, SLANKT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SCOPY, SKTEQR, XERBLA
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
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -6
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SKTEV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         D(1) = ZERO
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
      ISCALE = 0
      TNRM = SLANKT( 'M', N, E )
      IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / TNRM
      ELSE IF( TNRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / TNRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL SSCAL( N-1, SIGMA, E( 1 ), 1 )
      END IF
*
*     call SKTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL SKTEQR( 'N', N, E, Z, LDZ, WORK, INFO )
      ELSE
         CALL SKTEQR( 'I', N, E, Z, LDZ, WORK, INFO )
      END IF
*
      CALL SCOPY(N-1, E, 1, D, 1)
      D(N) = ZERO
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, D, 1 )
      END IF
*
      RETURN
*
*     End of SKTEV
*
      END
