*> \brief <b> DKTEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DKTEVX + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dktevx.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dktevx.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dktevx.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DKTEVX( JOBZ, RANGE, N, E, VL, VU, IL, IU, ABSTOL,
*                          M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, RANGE
*       INTEGER            IL, INFO, IU, LDZ, M, N
*       DOUBLE PRECISION   ABSTOL, VL, VU
*       ..
*       .. Array Arguments ..
*       INTEGER            IFAIL( * ), IWORK( * )
*       DOUBLE PRECISION   E( * ), W( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DKTEVX computes selected eigenvalues and, optionally, eigenvectors
*> of a real skew-symmetric tridiagonal matrix A.  Eigenvalues and
*> eigenvectors can be selected by specifying either a range of values
*> or a range of indices for the desired eigenvalues.
*>
*> Since the eigenvalues of skew-symmetric matrix are conjugate purely
*> imaginary pairs, to be consistent with output layout of subroutine
*> skyev, the output layout according to VL/VU and IL/IU abides the
*> following rules:
*>
*> If RANGE = 'A', the same as layout of skyev.
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
*> purely imaginary values, and are organized in the same layout as skyev.
*>
*> If RANGE = 'I' and N is even, the same as N being odd, except that the index
*> starts from the first element beside the center in the positive direction,
*> and all the indexes correspond to two eigenvalues.
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
*> \param[in] RANGE
*> \verbatim
*>          RANGE is CHARACTER*1
*>          = 'A': all eigenvalues will be found.
*>          = 'V': all eigenvalues in the half-open interval (VL,VU]
*>                 will be found.
*>          = 'I': the IL-th through IU-th eigenvalues will be found.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (max(1,N-1))
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal
*>          matrix A in elements 1 to N-1 of E.
*>          On exit, E may be multiplied by a constant factor chosen
*>          to avoid over/underflow in computing the eigenvalues.
*> \endverbatim
*>
*> \param[in] VL
*> \verbatim
*>          VL is DOUBLE PRECISION
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
*>          The absolute error tolerance for the eigenvalues.
*>          An approximate eigenvalue is accepted as converged
*>          when it is determined to lie in an interval [a,b]
*>          of width less than or equal to
*>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,
*>
*>          where EPS is the machine precision.  If ABSTOL is less
*>          than or equal to zero, then  EPS*|T|  will be used in
*>          its place, where |T| is the 1-norm of the tridiagonal
*>          matrix.
*>
*>          Eigenvalues will be computed most accurately when ABSTOL is
*>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*>          If this routine returns with INFO>0, indicating that some
*>          eigenvectors did not converge, try setting ABSTOL to
*>          2*DLAMCH('S').
*>
*>          See "Computing Small Singular Values of Bidiagonal Matrices
*>          with Guaranteed High Relative Accuracy," by Demmel and
*>          Kahan, LAPACK Working Note #3.
*> \endverbatim
*>
*> \param[out] M
*> \verbatim
*>          M is INTEGER
*>          The total number of eigenvalues found.  0 <= M <= N.
*>          If RANGE = 'A', M = N. If RANGE = 'I', when N is
*>          odd and IL = 1, M = 2 * (IU-IL) + 1, otherwise M = 
*>          2 * (IU-IL+1).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is DOUBLE PRECISION array, dimension (N)
*>          The first M elements contain the selected eigenvalues in
*>          ascending order.
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M) )
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*>          contain the block form orthonormal eigenvectors of the
*>          matrix A corresponding to the selected eigenvalues, with the
*>          i-th column of Z holding the block form eigenvector associated
*>          with W(i). The block form eigenvectors (v1, v2) corresponding
*>          to a conjugate purely imaginary eigenvalue pair (ri, -ri) means
*>          A*v1 = r*v2 and A*v2 = -r*v1.
*>          If an eigenvector fails to converge (INFO > 0), then that
*>          column of Z contains the latest approximation to the
*>          eigenvector, and the index of the eigenvector is returned
*>          in IFAIL.  If JOBZ = 'N', then Z is not referenced.
*>          Note: the user must ensure that at least max(1,M) columns are
*>          supplied in the array Z; if RANGE = 'V', the exact value of M
*>          is not known in advance and an upper bound must be used.
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
*>          WORK is DOUBLE PRECISION array, dimension (6*N)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (5*N)
*> \endverbatim
*>
*> \param[out] IFAIL
*> \verbatim
*>          IFAIL is INTEGER array, dimension (N)
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*>          indices of the eigenvectors that failed to converge.
*>          If JOBZ = 'N', then IFAIL is not referenced.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, then i eigenvectors failed to converge.
*>                Their indices are stored in array IFAIL.
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
*> \ingroup stevx
*
*  =====================================================================
      SUBROUTINE DKTEVX( JOBZ, RANGE, N, E, VL, VU, IL, IU,
     $                   ABSTOL,
     $                   M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      INTEGER            IL, INFO, IU, LDZ, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   E( * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, TEST, VALEIG, WANTZ
      CHARACTER          ORDER
      INTEGER            I, IMAX, INDISP, INDIWO, INDWRK,
     $                   ISCALE, ITMP1, J, JJ, NSPLIT, K
      DOUBLE PRECISION   BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM,
     $                   TMP1, TNRM, VLL, VUU
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANKT
      EXTERNAL           LSAME, DLAMCH, DLANKT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DSCAL, DKTEBZ, DKTEIN, DKTEQR
     $                   DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
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
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL .OR. VU.LT.ZERO )
     $         INFO = -7
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, (N+1)/2 ) ) THEN
               INFO = -8
            ELSE IF( IU.LT.MIN( (N+1)/2, IL ) .OR. IU.GT.(N+1)/2 )
     $      THEN
               INFO = -9
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) )
     $      INFO = -14
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DKTEVX', -INFO )
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
            W( 1 ) = ZERO
         ELSE
            IF( VL.LT.ZERO ) THEN
               M = 1
               W( 1 ) = ZERO
            END IF
         END IF
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
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
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      IF ( VALEIG ) THEN
         VLL = VL
         VUU = VU
      ELSE
         VLL = ZERO
         VUU = ZERO
      ENDIF
      TNRM = DLANKT( 'M', N, E )
      IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / TNRM
      ELSE IF( TNRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / TNRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N-1, SIGMA, E( 1 ), 1 )
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF
*
*     If all eigenvalues are desired and ABSTOL is less than zero, then
*     call DKTEQR.  If this fails for some eigenvalue, then
*     try DKTEBZ.
*
      TEST = .FALSE.
      IF( INDEIG ) THEN
         IF( IL.EQ.1 .AND. IU.EQ.(N+1)/2 ) THEN
            TEST = .TRUE.
         END IF
      END IF
*      IF( ( ALLEIG .OR. TEST ) .AND. ( ABSTOL.LE.ZERO ) ) THEN
      IF( .FALSE. ) THEN
         CALL DCOPY( N-1, E( 1 ), 1, W, 1 )
         INDWRK = N + 1
         CALL DKTEQR( 'I', N, W, Z, LDZ, WORK( INDWRK ),
     $                INFO )
         IF( INFO.EQ.0 ) THEN
            DO 10 I = 1, N
               IFAIL( I ) = 0
   10         CONTINUE
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 20
         END IF
         INFO = 0
      END IF
*
*     Otherwise, call DKTEBZ and, if eigenvectors are desired, DKTEIN.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDWRK = 1
      INDISP = 1 + N
      INDIWO = INDISP + N
      CALL DKTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTOL, E,
     $             M,
     $             NSPLIT, W, IWORK( 1 ), IWORK( INDISP ),
     $             WORK( INDWRK ), IWORK( INDIWO ), INFO )
*
      IF( WANTZ ) THEN
         CALL DKTEIN( N, E, M, W, IWORK( 1 ), IWORK( INDISP ),
     $                Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL,
     $                INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
   20 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.
*
      IF( WANTZ ) THEN
         JJ = 1
         DO WHILE(JJ.LT.(M-1))
            IF(W(JJ).EQ.ZERO) THEN
               DO K = JJ+1,M-1,2
                  IF(W(K).EQ.ZERO) THEN
                     DO I = JJ, K-2
                        W(I) = W(I+1)
                        CALL DSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ),
     $                              1 )
                        IF( INFO.NE.0 ) THEN
                           ITMP1 = IFAIL( I )
                           IFAIL( I ) = IFAIL( I+1 )
                           IFAIL( I+1 ) = ITMP1
                        END IF
                     END DO
                     W(K-1) = ZERO
                     JJ = K+1
                     EXIT
                  ELSEIF(MOD(M,2).EQ.1 .AND. K.EQ.(M-1)) THEN
                     DO I = JJ, K-1
                        W(I) = W(I+1)
                        CALL DSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ),
     $                              1 )
                        IF( INFO.NE.0 ) THEN
                           ITMP1 = IFAIL( I )
                           IFAIL( I ) = IFAIL( I+1 )
                           IFAIL( I+1 ) = ITMP1
                        END IF
                     END DO
                     CALL DSWAP( N, Z( 1, K ), 1, Z( 1, K+1 ),
     $                           1 )
                     IF( INFO.NE.0 ) THEN
                        ITMP1 = IFAIL( K )
                        IFAIL( K ) = IFAIL( K+1 )
                        IFAIL( K+1 ) = ITMP1
                     END IF
                     W(K) = ZERO
                     JJ = K+1
                     EXIT
                  ELSEIF(MOD(M,2).EQ.0 .AND. K.EQ.(M-2)) THEN
                     DO I = JJ, K-1
                        W(I) = W(I+1)
                        CALL DSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ),
     $                              1 )
                        IF( INFO.NE.0 ) THEN
                           ITMP1 = IFAIL( I )
                           IFAIL( I ) = IFAIL( I+1 )
                           IFAIL( I+1 ) = ITMP1
                        END IF
                     END DO
                     CALL DSWAP( N, Z( 1, K ), 1, Z( 1, K+1 ), 1 )
                     IF( INFO.NE.0 ) THEN
                        ITMP1 = IFAIL( K )
                        IFAIL( K ) = IFAIL( K+1 )
                        IFAIL( K+1 ) = ITMP1
                     END IF
                     W(K) = ZERO
                     JJ = K+1
                     EXIT
                  END IF
               END DO
               IF (JJ.LT.(M-1)) THEN
                  CYCLE
               END IF
            END IF
            JJ = JJ+2
         END DO
*
         DO 40 JJ = 1, M-1, 2
            I = JJ
            TMP1 = ABS(W(JJ))
            DO 30 K = JJ+2, M-1, 2
               IF(ABS(W(K)).GT.TMP1) THEN
                  I = K
                  TMP1 = ABS(W(K))
               END IF
   30       CONTINUE
            IF(I.NE.JJ) THEN
               CALL DSWAP( 1, W( I ), 1, W( JJ ), 1 )
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, JJ ), 1 )
               CALL DSWAP( N, Z( 1, I+1 ), 1, Z( 1, JJ+1 ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( JJ )
                  IFAIL( JJ ) = ITMP1
                  ITMP1 = IFAIL( I+1 )
                  IFAIL( I+1 ) = IFAIL( JJ+1 )
                  IFAIL( JJ+1 ) = ITMP1
               END IF
            END IF
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DKTEVX
*
      END
