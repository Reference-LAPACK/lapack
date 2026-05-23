*> \brief \b ZGGEV4
*
* =========== DOCUMENTATION ===========
*
*> \par Purpose:
* =============
*>
*> \details \b Purpose:
*> \verbatim
*>
*> ZGGEV4 computes for a pair of N-by-N complex nonsymmetric matrices 
*> (A,B) the generalized eigenvalues, and optionally, the left and/or 
*> right generalized eigenvectors.
*> 
*> This routine performs the following steps:
*> 1) Preprocessing (ZGGPRP) 
*> 2) Deflation of infinite eigenvalues (ZGGDEF) 
*> 3) Blocked Hessenberg-Triangular Reduction (ZGGHD4) 
*> 4) QZ algorithm (ZLAQZ0) 
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
*>          The order of the matrices A, B, VL, and VR.
*>          N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, N)
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
*>          B is COMPLEX*16 array, dimension (LDB, N)
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
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16 array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is COMPLEX*16 array, dimension (N)
*>          ALPHA(j)/BETA(j), j=1,...,N, will be the
*>          generalized eigenvalues.
*>          Deflated infinite eigenvalues will set BETA to zero.
*> \endverbatim
*>
*> \param[out] VL
*> \verbatim
*>          VL is COMPLEX*16 array, dimension (LDVL, N)
*>          If JOBVL = 'V', the left eigenvectors are stored here.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the array VL.
*>          LDVL >= 1, and if JOBVL = 'V', LDVL >= N.
*> \endverbatim
*>
*> \param[out] VR
*> \verbatim
*>          VR is COMPLEX*16 array, dimension (LDVR, N)
*>          If JOBVR = 'V', the right eigenvectors are stored here.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR.
*>          LDVR >= 1, and if JOBVR = 'V', LDVR >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
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
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (MAX(1, 8*N))
*>          Workspace required for real magnitudes and operations.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N)
*>          Pivot indices from preprocessing.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
      SUBROUTINE ZGGEV4( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA,
     $                   BETA, VL, LDVL, VR, LDVR, WORK,
     $                   LWORK, RWORK, IWORK, INFO )
*
* .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
* ..
* .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )
* ..
*
* =====================================================================
*
* .. Parameters ..
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0D+0, RONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE, X
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), 
     $                     CONE  = ( 1.0D+0, 0.0D+0 ) )
* ..
* .. Local Scalars ..
      LOGICAL            COMPVL, COMPVR, LQUERY
      CHARACTER          JOBQZ, JOBVEC
*
* IHI is added to support balancing: ZGGBAL sets ILO and IHI to
* indicate the submatrix that was balanced. IHI is not needed by
* the existing pipeline but is required by ZGGBAK.
*
      INTEGER            I, IN, J, K, NINFINITE, WORKNEEDED, ILO
      DOUBLE PRECISION   RTMP, TOL
* ..
* .. Local Arrays ..
      LOGICAL            LDUMMY( 1 )
      COMPLEX*16         DUMMY( 1 )
* ..
* .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
* ..
* .. External Subroutines ..
*
* ZGGBAL / ZGGBAK added for balancing (same as ZGGEV3).
*
      EXTERNAL           ZGEQRF, ZGERQF, ZGGBAL, ZGGBAK, ZGGDEF,
     $                   ZGGHD4, ZGGPRP, ZLAQZ0, ZUNMQR, ZUNMRQ,
     $                   ZTGEVC, XERBLA
* ..
* .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, INT, MAX, MIN
* ..
* .. Executable Statements ..
*
      ABS1(X) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )

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
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( COMPVR .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGEV4', -INFO )
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
      CALL ZGGPRP( COMPVL, COMPVR, N, A, LDA, B, LDB, VL, LDVL, VR,
     $             LDVR, DUMMY, -1, RWORK, IWORK, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DBLE( DUMMY( 1 ) ) ) )
      CALL ZGGHD4( JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, VL, LDVL,
     $             VR, LDVR, DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DBLE( DUMMY( 1 ) ) ) )
      CALL ZGEQRF( N, N, A, LDA, DUMMY, DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DBLE( DUMMY( 1 ) ) ) )
      CALL ZUNMQR( 'L', 'C', N, N, N, A, LDA, DUMMY, DUMMY, N,
     $             DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DBLE( DUMMY( 1 ) ) ) )
      CALL ZGERQF( N, N, A, LDA, DUMMY, DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DBLE( DUMMY( 1 ) ) ) )
      CALL ZUNMRQ( 'L', 'C', N, N, N, A, LDA, DUMMY, DUMMY, N,
     $             DUMMY, -1, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DBLE( DUMMY( 1 ) ) ) )
      CALL ZLAQZ0( JOBQZ, JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB,
     $             ALPHA, BETA, VL, LDVL, VR, LDVR, DUMMY,
     $             -1, RWORK, 0, INFO )
      WORKNEEDED = MAX( WORKNEEDED, INT( DBLE( DUMMY( 1 ) ) ) )
*
* ZGGBAL does not use a workspace query; it always needs 6*N
* elements. Ensure the workspace is large enough to accommodate it.
*
      WORKNEEDED = MAX( WORKNEEDED, 6*N )
*
      IF( LQUERY ) THEN
         WORK( 1 ) = DCMPLX( DBLE( WORKNEEDED ), 0.0D+0 )
         RETURN
      END IF

* ---------------------------------------------------------------
*
* Step 1: Preprocessing
      CALL ZGGPRP( COMPVL, COMPVR, N, A, LDA, B, LDB, VL, LDVL, VR,
     $             LDVR, WORK, LWORK, RWORK, IWORK, INFO )

      TOL = RWORK( 2 )
      NINFINITE = 0
      DO K = 1, N
         IF (TOL.LT.ABS(B(K, K))) THEN
            EXIT
         END IF
         NINFINITE = NINFINITE + 1
         DO I = 1, NINFINITE
            B(I, K) = CZERO
         END DO
      END DO

*
* Step 2: Deflation of infinite eigenvalues
      IF( NINFINITE.GT.0 ) THEN
         CALL ZGGDEF( COMPVL, COMPVR, N, NINFINITE, A, LDA, B, LDB,
     $                VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
      END IF

      DO I = 1, NINFINITE
         ALPHA(I) = A(I, I)
         BETA(I) = CZERO
      END DO

*
* Add small perturbation to B if needed
      DO I = NINFINITE + 1, N
         IF (ABS(B(I, I)).LE.TOL) THEN
            B(I, I) = DCMPLX( TOL, RZERO )
         END IF
      END DO
*
* ---------------------------------------------------------------
*
* Step 3: Hessenberg-Triangular Reduction
*
      ILO = MIN( N, NINFINITE + 1 )
      CALL ZGGHD4( JOBVL, JOBVR, N, ILO, N, A, LDA, B, LDB,
     $             VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

*
* ---------------------------------------------------------------
*
* Step 4: QZ
* Use the balanced ILO (not the raw NINFINITE+1) so that ZLAQZ0
* operates on exactly the same submatrix that ZGGHD4 reduced to
* Hessenberg-triangular form. Rows NINFINITE+1..ILO-1 were
* identified by ZGGBAL as already decoupled and were skipped by
* ZGGHD4; passing NINFINITE+1 here would hand ZLAQZ0 a
* non-Hessenberg leading block.
      IF( NINFINITE.NE.N ) THEN
         CALL ZLAQZ0( JOBQZ, JOBVL, JOBVR, N, ILO, N, A, LDA,
     $                B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR,
     $                WORK, LWORK, RWORK, 0, INFO )
      END IF
*
* Step 5: Eigenvectors
      IF( COMPVL .OR. COMPVR ) THEN
         CALL ZTGEVC( JOBVEC, 'B', LDUMMY, N, A, LDA, B, LDB, VL, LDVL,
     $                VR, LDVR, N, IN, WORK, RWORK, INFO )
      END IF
*
* ---------------------------------------------------------------
*
* Normalise left eigenvectors
      IF( COMPVL ) THEN
         DO I = 1, N
            RTMP = RZERO
            DO J = 1, N
               RTMP = MAX( RTMP, ABS1( VL( J, I ) ) )
            END DO
            RTMP = RONE / RTMP
            DO J = 1, N
               VL( J, I ) = VL( J, I ) * DCMPLX( RTMP, RZERO )
            END DO
         END DO
      END IF
*
* Normalise right eigenvectors
      IF( COMPVR ) THEN
         DO I = 1, N
            RTMP = RZERO
            DO J = 1, N
               RTMP = MAX( RTMP, ABS1( VR( J, I ) ) )
            END DO
            RTMP = RONE / RTMP
            DO J = 1, N
               VR( J, I ) = VR( J, I ) * DCMPLX( RTMP, RZERO )
            END DO
         END DO
      END IF
*
      WORK( 1 ) = DCMPLX( DBLE( NINFINITE ), RZERO )
      RETURN
      END SUBROUTINE ZGGEV4