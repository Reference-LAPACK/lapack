*> \brief \b DGGHD4
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DGGHD4 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgghd4.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgghd4.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgghd4.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGGHD4( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
*                          LDQ, Z, LDZ, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          COMPQ, COMPZ
*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
*      $                   Z( LDZ, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGGHD4 reduces a pair of real matrices (A,B) to generalized upper
*> Hessenberg form using orthogonal transformations, where A is a
*> general matrix and B is upper triangular.  The form of the
*> generalized eigenvalue problem is
*>    A*x = lambda*B*x,
*> and B is typically made upper triangular by computing its QR
*> factorization and moving the orthogonal matrix Q to the left side
*> of the equation.
*>
*> This subroutine simultaneously reduces A to a Hessenberg matrix H:
*>    Q**T*A*Z = H
*> and transforms B to another upper triangular matrix T:
*>    Q**T*B*Z = T
*> in order to reduce the problem to its standard form
*>    H*y = lambda*T*y
*> where y = Z**T*x.
*>
*> The orthogonal matrices Q and Z are determined as products of Givens
*> rotations.  They may either be formed explicitly, or they may be
*> postmultiplied into input matrices Q1 and Z1, so that
*>
*>      Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T
*>
*>      Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T
*>
*> If Q1 is the orthogonal matrix from the QR factorization of B in the
*> original equation A*x = lambda*B*x, then DGGHD4 reduces the original
*> problem to generalized Hessenberg form.
*>
*> This is a blocked variant of DGGHRD, using matrix-matrix
*> multiplications for parts of the computation to enhance performance.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] COMPQ
*> \verbatim
*>          COMPQ is CHARACTER*1
*>          = 'N': do not compute Q;
*>          = 'I': Q is initialized to the unit matrix, and the
*>                 orthogonal matrix Q is returned;
*>          = 'V': Q must contain an orthogonal matrix Q1 on entry,
*>                 and the product Q1*Q is returned.
*> \endverbatim
*>
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N': do not compute Z;
*>          = 'I': Z is initialized to the unit matrix, and the
*>                 orthogonal matrix Z is returned;
*>          = 'V': Z must contain an orthogonal matrix Z1 on entry,
*>                 and the product Z1*Z is returned.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>          ILO and IHI mark the rows and columns of A which are to be
*>          reduced.  It is assumed that A is already upper triangular
*>          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
*>          normally set by a previous call to DGGBAL; otherwise they
*>          should be set to 1 and N respectively.
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          On entry, the N-by-N general matrix to be reduced.
*>          On exit, the upper triangle and the first subdiagonal of A
*>          are overwritten with the upper Hessenberg matrix H, and the
*>          rest is set to zero.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB, N)
*>          On entry, the N-by-N upper triangular matrix B.
*>          On exit, the upper triangular matrix T = Q**T B Z.  The
*>          elements below the diagonal are set to zero.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
*>          On entry, if COMPQ = 'V', the orthogonal matrix Q1,
*>          typically from the QR factorization of B.
*>          On exit, if COMPQ='I', the orthogonal matrix Q, and if
*>          COMPQ = 'V', the product Q1*Q.
*>          Not referenced if COMPQ='N'.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.
*>          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*>          On entry, if COMPZ = 'V', the orthogonal matrix Z1.
*>          On exit, if COMPZ='I', the orthogonal matrix Z, and if
*>          COMPZ = 'V', the product Z1*Z.
*>          Not referenced if COMPZ='N'.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.
*>          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in]  LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of the array WORK.  LWORK >= 1.
*>          For optimum performance LWORK >= 6*N*NB, where NB is the
*>          optimal blocksize.
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
*>          > 0:  the refinement scheme failed to converge.
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Thijs Steel
*> \author Raf Vandebril
*
*> \ingroup doubleOTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  This routine reduces A to Hessenberg form and maintains B in triangular form
*>  using an iterative variant of Moler and Stewart's original algorithm,
*>  as described by Thijs Steel and Raf Vandebril (ELA 2022)
*>
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGGHD4( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
     $                   LDQ, Z, LDZ, WORK, LWORK, INFO )
      IMPLICIT NONE
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ,
     $                   * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, TEN
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, TEN = 1.0D+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILQ, ILZ, XOVERFLOW
      CHARACTER          COMPQ2, COMPZ2
      INTEGER            ICOMPQ, ICOMPZ, JCOL, JROW, ITAU, IWRK, IERR,
     $                   COUNT, LWORKREQ, ILO2, COUNT2, IMAX
      DOUBLE PRECISION   ERR, ANRM, BNRM, ULP, TOL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLANGE, DLAMCH
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DLANGE, DLAMCH, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARTG, DLASET, DROT, XERBLA, DLAHT0
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Decode COMPQ
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
*     Decode COMPZ
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -11
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGHD4', -INFO )
         RETURN
      END IF
*
*     Find out required workspace
*
      CALL DORMQR( 'L', 'T', N, N, N, B, LDB, WORK, A, LDA, WORK, -1,
     $             IERR )
      LWORKREQ = INT( WORK( 1 ) )+N+N**2
      CALL DGEHRD( N, 1, N, WORK, N, WORK, WORK, -1, IERR )
      LWORKREQ = MAX ( INT( WORK( 1 ) )+N+N**2, LWORKREQ )
      CALL DORMHR( 'R', 'N', N, N, 1, N, WORK, N, WORK, A, LDA, WORK,
     $             -1, IERR )
      LWORKREQ = MAX ( INT( WORK( 1 ) )+N+N**2, LWORKREQ )
      CALL DGEQRF( N, N, B, LDB, WORK, WORK, -1, IERR )
      LWORKREQ = MAX ( INT( WORK( 1 ) )+N, LWORKREQ )
      CALL DLAHT0( COMPQ, COMPZ, N, ILO, IHI, IHI, A, LDA, B, LDB, Q,
     $             LDQ, Z, LDZ, WORK, -1, IERR )
      LWORKREQ = MAX ( INT( WORK( 1 ) ), LWORKREQ )
*
      IF ( LWORK .EQ.-1 ) THEN
         WORK( 1 ) = DBLE( LWORKREQ )
         RETURN
      ELSE IF ( LWORK .LT. LWORKREQ ) THEN
         INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGHD4', -INFO )
         RETURN
      END IF
*
*     Initialize Q and Z if desired.
*
      IF( ICOMPQ.EQ.3 )CALL DLASET( 'FULL', N, N, ZERO, ONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )CALL DLASET( 'FULL', N, N, ZERO, ONE, Z, LDZ )
*
*     Quick return if possible
*
      IF( N.LE.1 )RETURN
*
*     Determine machine parameters
*
      ANRM = DLANGE( 'F', IHI-ILO+1, IHI-ILO+1, A( ILO, ILO ), LDA,
     $   WORK )
      BNRM = DLANGE( 'F', IHI-ILO+1, IHI-ILO+1, B( ILO, ILO ), LDB,
     $   WORK )
      ULP = DLAMCH( 'PRECISION' )
      COUNT = 0
      ILO2 = ILO

*
*     Start of the main loop
*
   40 CONTINUE
*
*     If too many iterations have failed, return with error.
*     Because we default to reducing some columns using DLAHT0,
*     this should be impossible.
*
      IF ( COUNT .GE. 50 ) THEN
         INFO = 1
         GOTO 60
      END IF
      COUNT = COUNT+1
*
*     Reduce an increasing amount of columns using the fallback algorithm
*     DLAHT0. This guarantees convergence.
*     Avoid re-initialization of modified Q and Z.
*
      IF( COUNT .GT. 1) THEN
         IMAX = MIN(IHI, ILO2 + 32*(2**(COUNT-2) - 1))
         COMPQ2 = COMPQ
         COMPZ2 = COMPZ
         IF ( ILQ )
     $      COMPQ2 = 'V'
         IF ( ILZ )
     $      COMPZ2 = 'V'
         CALL DLAHT0( COMPQ2, COMPZ2, N, ILO2, IHI, IMAX, A, LDA,
     $                B, LDB, Q, LDQ, Z, LDZ, WORK, LWORK, IERR )
         ILO2 = IMAX
         IF( ILO2 .GE. IHI ) THEN
            GOTO 60
         END IF
      END IF
*
*     Calculate X = AB^{-1}
*     Perturb B slightly to avoid overflow.
*
      ITAU = ( IHI-ILO2+1 )**2+1
      DO JROW = ILO2, IHI
         WORK( ITAU+JROW-ILO2 ) = B( JROW, JROW )
      END DO

      TOL = 1.0E-15
      COUNT2 = 0
   50 CONTINUE
      CALL DLACPY( 'A', IHI-ILO2+1, IHI-ILO2+1, A( ILO2, ILO2 ), LDA,
     $             WORK, IHI-ILO2+1 )
      CALL DTRSM( 'R', 'U', 'N', 'N', IHI-ILO2+1, IHI-ILO2+1, ONE,
     $            B( ILO2, ILO2 ), LDB, WORK, IHI-ILO2+1 )
      XOVERFLOW = .FALSE.
      DO JROW = 1, (IHI-ILO+1)**2
         IF( DISNAN( WORK(JROW) ) ) THEN
            XOVERFLOW = .TRUE.
            EXIT
         END IF
      END DO
*
*     Calculation of X has overflowed, try again with a perturbed B
*
      IF( XOVERFLOW ) THEN
         IF ( COUNT2 .GE. 5 ) THEN
*           Calculation of X has failed too many times, return with error.
            INFO = 1
            GOTO 60
         END IF
         DO JROW = ILO2, IHI
            IF( ABS( B( JROW, JROW ) ) .LT. TOL*BNRM ) THEN
               B( JROW, JROW ) = TOL*BNRM
            END IF
         END DO
         COUNT2 = COUNT2 + 1
         TOL = TOL * 1.0E+1
         GOTO 50
      END IF
*     Restore original B
      DO JROW = ILO2, IHI
         B( JROW, JROW ) = WORK( ITAU+JROW-ILO2 )
      END DO
*
*     Calculate Q so that Q**T X Q = H
*
      ITAU = ( IHI-ILO2+1 )**2+1
      IWRK = ITAU+IHI-ILO2
      CALL DGEHRD( IHI-ILO2+1, 1, IHI-ILO2+1, WORK, IHI-ILO2+1,
     $             WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )
*
*     Apply Q**T to A and B
*
      CALL DORMHR( 'L', 'T', IHI-ILO2+1, N-ILO2+1, 1, IHI-ILO2+1, WORK,
     $             IHI-ILO2+1, WORK( ITAU ), A( ILO2, ILO2 ), LDA,
     $             WORK( IWRK ), LWORK-IWRK+1, IERR )
      CALL DORMHR( 'L', 'T', IHI-ILO2+1, N-ILO2+1, 1, IHI-ILO2+1, WORK,
     $             IHI-ILO2+1, WORK( ITAU ), B( ILO2, ILO2 ), LDB,
     $             WORK( IWRK ), LWORK-IWRK+1, IERR )
      IF ( ILQ ) CALL DORMHR( 'R', 'N', N, IHI-ILO2+1, 1, IHI-ILO2+1,
     $     WORK, IHI-ILO2+1, WORK( ITAU ), Q( 1, ILO2 ), LDQ,
     $     WORK( N**2+N ), LWORK-IWRK+1, IERR )
*
*     Make B upper triangular
*
      ITAU = 1
      IWRK = ITAU+IHI-ILO2+1

      CALL DGERQF( IHI, IHI-ILO2+1, B( 1, ILO2 ), LDB, WORK( ITAU ),
     $             WORK( IWRK ), LWORK-IWRK+1, IERR )

      CALL DORMRQ( 'R', 'T', IHI, IHI-ILO2+1, IHI-ILO2+1, B( ILO2,
     $             ILO2 ), LDB, WORK( ITAU ), A( 1, ILO2 ), LDA,
     $             WORK( IWRK ), LWORK-IWRK+1, IERR )
      IF ( ILZ ) THEN
         CALL DORMRQ( 'R', 'T', N, IHI-ILO2+1, IHI-ILO2+1, B( ILO2,
     $                ILO2 ), LDB, WORK( ITAU ), Z( 1, ILO2 ), LDZ,
     $                WORK( IWRK ), LWORK-IWRK+1, IERR )
      END IF
      DO JCOL = ILO2, IHI-1
         DO JROW = JCOL+1, IHI
            B( JROW, JCOL ) = ZERO
         END DO
      END DO
*
*     Calculate residual in A
*
      DO WHILE( ILO2 .LT. IHI-1 )
         ERR = ZERO
         DO JROW = ILO2+2, IHI
            ERR = ERR + ABS( A( JROW, ILO2 ) )
         END DO
         IF( ERR .GT. 0.5D2*ANRM*ULP ) GO TO 40
         ILO2 = ILO2 + 1
      END DO

   60 CONTINUE

      DO JCOL = 1, N-1
         DO JROW = JCOL+2, N
            A( JROW, JCOL ) = ZERO
         END DO
      END DO

      RETURN
      END
