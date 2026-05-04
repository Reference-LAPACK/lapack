*> \brief \b DGGEV4
*
* =========== DOCUMENTATION ===========
*
*> \par Purpose:
* =============
*>
*> \details \b Purpose:
*> \verbatim
*>
*> DGGEV4 computes for a pair of N-by-N real nonsymmetric matrices (A,B)
*> the generalized eigenvalues, and optionally, the left and/or right
*> generalized eigenvectors.
*> 
*> This routine performs the following steps:
*> 1) Preprocessing (DGGPRP) 
*> 2) Deflation of infinite eigenvalues (DGGDEF) 
*> 3) Blocked Hessenberg-Triangular Reduction (DGGHD4) 
*> 4) QZ algorithm (DLAQZ0) 
*> 5) Eigenvector computation and normalization.
*> \endverbatim
*
* Arguments:
* ==========
*
*> \param[in] JOBVL
*> \verbatim
*>          JOBVL is CHARACTER*1
*>          = 'N':  do not compute the left generalized eigenvectors;
*>          = 'V':  compute the left generalized eigenvectors.
*> \endverbatim
*>
*> \param[in] JOBVR
*> \verbatim
*>          JOBVR is CHARACTER*1
*>          = 'N':  do not compute the right generalized eigenvectors;
*>          = 'V':  compute the right generalized eigenvectors.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VL, and VR.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          On entry, the matrix A.
*>          On exit, A has been overwritten.
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
*>          On entry, the matrix B.
*>          On exit, B has been overwritten.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] ALPHAR
*> \verbatim
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] ALPHAI
*> \verbatim
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>          (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will be the
*>          generalized eigenvalues. Deflated infinite eigenvalues 
*>          will set BETA to zero.
*> \endverbatim
*>
*> \param[out] VL
*> \verbatim
*>          VL is DOUBLE PRECISION array, dimension (LDVL, N)
*>          If JOBVL = 'V', the left eigenvectors are stored here.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the array VL. LDVL >= 1, and
*>          if JOBVL = 'V', LDVL >= N.
*> \endverbatim
*>
*> \param[out] VR
*> \verbatim
*>          VR is DOUBLE PRECISION array, dimension (LDVR, N)
*>          If JOBVR = 'V', the right eigenvectors are stored here.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR. LDVR >= 1, and
*>          if JOBVR = 'V', LDVR >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the number of deflated 
*>          infinite eigenvalues.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If LWORK = -1, then a workspace query is assumed.
*> \endverbatim
*>
*> \param[out] JPVT
*> \verbatim
*>          JPVT is INTEGER array, dimension (N)
*>          Pivot indices from preprocessing.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
      SUBROUTINE DGGEV4( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR,
     $                   ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK,
     $                   LWORK, IWORK, INFO )
*
* .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
* ..
* .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
* ..
*
* =====================================================================
*
* .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
* ..
* .. Local Scalars ..
      LOGICAL            COMPVL, COMPVR, LQUERY
      CHARACTER          JOBQZ, JOBVEC
*
*     IHI is added to support balancing: DGGBAL sets ILO and IHI to
*     indicate the submatrix that was balanced.  IHI is not needed by
*     the existing pipeline but is required by DGGBAK.
*
      INTEGER            I, IN, J, K, NINFINITE, WORKNEEDED, ILO
      DOUBLE PRECISION   TMP, TOL
* ..
* .. Local Arrays ..
      LOGICAL            LDUMMY( 1 )
      DOUBLE PRECISION   DUMMY( 1 )
* ..
* .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
* ..
* .. External Subroutines ..
*
*     DGGBAL / DGGBAK added for balancing (same as DGGEV3).
*
      EXTERNAL           DGEQRF, DGERQF, DGGBAL, DGGBAK, DGGDEF,
     $                   DGGHD4, DGGPRP, DLAQZ0, DORMQR, DORMRQ,
     $                   DTGEVC, XERBLA
* ..
* .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, SQRT, MIN
* ..
* .. Executable Statements ..
*
      COMPVL = LSAME( JOBVL, 'V' )
      COMPVR = LSAME( JOBVR, 'V' )
      IF( COMPVL .OR. COMPVR ) THEN
         JOBQZ = 'S'
      ELSE
         JOBQZ = 'E'
      END IF
*
* LAPACK style argument checks
      INFO = 0
      IF( .NOT.COMPVL .AND. .NOT.LSAME( JOBVL, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.COMPVR .AND. .NOT.LSAME( JOBVR, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( COMPVL .AND. LDVL.LT.N ) ) THEN
         INFO = -12
      ELSE IF( LDVR.LT.1 .OR. ( COMPVR .AND. LDVR.LT.N ) ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGEV4', -INFO )
         RETURN
      END IF
*
      IF( COMPVL .AND. COMPVR ) THEN
         JOBVEC = 'B'
      ELSE IF( COMPVL ) THEN
         JOBVEC = 'L'
      ELSE IF( COMPVR ) THEN
         JOBVEC = 'R'
      ELSE
         JOBVEC = 'N'
      END IF
*
* Workspace queries
      LQUERY = ( LWORK.EQ.-1 )
      WORKNEEDED = 0
*
      CALL DGGPRP( COMPVL, COMPVR, N, A, LDA, B, LDB, VL, LDVL, VR,
     $             LDVR, DUMMY, -1, IWORK, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DUMMY( 1 ) ) )
      CALL DGGHD4( JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, VL, LDVL,
     $             VR, LDVR, DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DUMMY( 1 ) ) )
      CALL DGEQRF( N, N, A, LDA, DUMMY, DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DUMMY( 1 ) ) )
      CALL DORMQR( 'L', 'T', N, N, N, A, LDA, DUMMY, DUMMY, N,
     $             DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DUMMY( 1 ) ) )
      CALL DGERQF( N, N, A, LDA, DUMMY, DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DUMMY( 1 ) ) )
      CALL DORMRQ( 'L', 'T', N, N, N, A, LDA, DUMMY, DUMMY, N,
     $             DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DUMMY( 1 ) ) )
      CALL DLAQZ0( JOBQZ, JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB,
     $             ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, DUMMY,
     $             -1, 0, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DUMMY( 1 ) ) )
*
*     DGGBAL does not use a workspace query; it always needs 6*N
*     elements.  Ensure the workspace is large enough to accommodate it.
*
      WORKNEEDED = MAX( WORKNEEDED, 6*N )
*
      IF( LQUERY ) THEN
         WORK( 1 ) = DBLE( WORKNEEDED )
         RETURN
      END IF

* ---------------------------------------------------------------
*
* Step 1: Preprocessing
      CALL DGGPRP( COMPVL, COMPVR, N, A, LDA, B, LDB, VL, LDVL, VR,
     $             LDVR, WORK, LWORK, IWORK, INFO )

      TOL = WORK( 2 )
      NINFINITE = 0
      DO K = 1, N
         IF (TOL.LT.ABS(B(K, K))) THEN
            EXIT
         END IF
         NINFINITE = NINFINITE + 1
         DO I = 1, NINFINITE
            B(I, K) = ZERO
         END DO
      END DO

*
* Step 2: Deflation of infinite eigenvalues
      IF( NINFINITE.GT.0 ) THEN
         CALL DGGDEF( COMPVL, COMPVR, N, NINFINITE, A, LDA, B, LDB,
     $                VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
      END IF

      DO I = 1, NINFINITE
         ALPHAR(I) = A(I, I)
         ALPHAI(I) = ZERO
         BETA(I) = ZERO
      END DO

*
* Add small perturbation to B if needed
      DO I = NINFINITE + 1, N
         IF (ABS(B(I, I)).LE.TOL) THEN
            B(I, I) = TOL
         END IF
      END DO
*
* ---------------------------------------------------------------
*
* Step 3: Hessenberg-Triangular Reduction
*
      ILO = MIN( N, NINFINITE + 1 )
      CALL DGGHD4( JOBVL, JOBVR, N, ILO, N, A, LDA, B, LDB,
     $             VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

*
* ---------------------------------------------------------------
*
* Step 4: QZ
*     Use the balanced ILO (not the raw NINFINITE+1) so that DLAQZ0
*     operates on exactly the same submatrix that DGGHD4 reduced to
*     Hessenberg-triangular form.  Rows NINFINITE+1..ILO-1 were
*     identified by DGGBAL as already decoupled and were skipped by
*     DGGHD4; passing NINFINITE+1 here would hand DLAQZ0 a
*     non-Hessenberg leading block.
      IF( NINFINITE.NE.N ) THEN
         CALL DLAQZ0( JOBQZ, JOBVL, JOBVR, N, ILO, N, A, LDA,
     $                B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR,
     $                WORK, LWORK, 0, INFO )
      END IF
*
* Step 5: Eigenvectors
      IF( COMPVL .OR. COMPVR ) THEN
         CALL DTGEVC( JOBVEC, 'B', LDUMMY, N, A, LDA, B, LDB, VL, LDVL,
     $                VR, LDVR, N, IN, WORK, INFO )
      END IF
*
* ---------------------------------------------------------------
*
* Normalise left eigenvectors
      IF( COMPVL ) THEN
         I = 1
         DO WHILE( I.LE.N )
            TMP = ZERO
            IF( ALPHAI( I ).NE.ZERO ) THEN
* Complex case
               DO J = 1, N
                  TMP = MAX(TMP, ABS(VL(J, I)) + ABS(VL(J, I + 1)))
               END DO
               TMP = ONE / TMP
               DO J = 1, N
                  VL( J, I )   = VL( J, I ) * TMP
                  VL( J, I+1 ) = VL( J, I+1 ) * TMP
               END DO
               I = I + 2
            ELSE
* Real case
               DO J = 1, N
                  TMP = MAX(TMP, ABS(VL(J, I)))
               END DO
               TMP = ONE / TMP
               DO J = 1, N
                  VL( J, I ) = VL( J, I ) * TMP
               END DO
               I = I + 1
            END IF
         END DO
      END IF
*
* Normalise right eigenvectors
      IF( COMPVR ) THEN
         I = 1
         DO WHILE( I.LE.N )
            TMP = ZERO
            IF( ALPHAI( I ).NE.ZERO ) THEN
* Complex case
               DO J = 1, N
                  TMP = MAX(TMP, ABS(VR(J, I)) + ABS(VR(J, I + 1)))
               END DO
               TMP = ONE / TMP
               DO J = 1, N
                  VR( J, I )   = VR( J, I ) * TMP
                  VR( J, I+1 ) = VR( J, I+1 ) * TMP
               END DO
               I = I + 2
            ELSE
* Real case
               DO J = 1, N
                  TMP = MAX(TMP, ABS(VR(J, I)))
               END DO
               TMP = ONE / TMP
               DO J = 1, N
                  VR( J, I ) = VR( J, I ) * TMP
               END DO
               I = I + 1
            END IF
         END DO
      END IF
*
      WORK( 1 ) = DBLE( NINFINITE )
      RETURN
      END SUBROUTINE DGGEV4