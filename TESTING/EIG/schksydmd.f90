!     This is a test program for checking the implementations of
!     the implementations of the following subroutines
!
!     SSYDMD  for computation of the
!             Dynamic Mode Decomposition (DMD)
!     SSYDMDQ for computation of a
!             QR factorization based compressed DMD
!
!     Developed and supported by:
!     ===========================
!     Developed and coded by Zlatko Drmac, Faculty of Science,
!     University of Zagreb;  drmac@math.hr
!     In cooperation with
!     AIMdyn Inc., Santa Barbara, CA.
!     ========================================================
!............................................................
!............................................................
!
      PROGRAM DMD_TEST
      use iso_fortran_env, only: real32
      IMPLICIT NONE
      integer, parameter :: WP = real32

!............................................................
      REAL(KIND=WP), PARAMETER ::  ONE = 1.0_WP
      REAL(KIND=WP), PARAMETER :: ZERO = 0.0_WP
!............................................................
      REAL(KIND=WP), ALLOCATABLE, DIMENSION(:,:) ::          &
                     A, AC, EIGA, LAMBDA, LAMBDAQ, F, F1, F2,&
                     Z, Z1, S, AU, W, VA, X, X0, Y, Y0, Y1
      REAL(KIND=WP), ALLOCATABLE, DIMENSION(:)   ::          &
                     DA, DL, DR, REIG, REIGA, REIGQ, IEIG,   &
                     IEIGA, IEIGQ,  RES, RES1, RESEX, SINGVX,&
                     SINGVQX, WORK
      INTEGER      , ALLOCATABLE, DIMENSION(:)   ::   IWORK
      REAL(KIND=WP) :: AB(2,2),   WDUMMY(2)
      INTEGER       :: IDUMMY(2), ISEED(4)
      REAL(KIND=WP) :: ANORM, COND, CONDL, CONDR, DMAX, EPS, &
                       TOL, TOL2, SVDIFF, TMP, TMP_AU,       &
                       TMP_FQR, TMP_REZ, TMP_REZQ,  TMP_ZXW, &
                       TMP_EX, XNORM, YNORM
!............................................................
      INTEGER :: K, KQ, LDF, LDS, LDA, LDAU, LDW, LDX, LDY,  &
                 LDZ, LIWORK, LWORK, M, N, L, LLOOP, NRNK
      INTEGER :: i, iJOBREF, iJOBZ, iSCALE, INFO, KDIFF,     &
                 NFAIL, NFAIL_AU, NFAIL_F_QR, NFAIL_REZ,     &
                 NFAIL_REZQ, NFAIL_SVDIFF, NFAIL_TOTAL, NFAILQ_TOTAL, &
                 NFAIL_Z_XV, MODE, MODEL, MODER, WHTEIG, WHTSVD, WHTSYM
      INTEGER    iNRNK, iWHTEIG, iWHTSVD, iWHTSYM, K_TRAJ, LWMINOPT
      CHARACTER(LEN=1) GRADE, JOBREF, JOBZ, PIVTNG, RSIGN,   &
                       SCALE, RESIDS, WANTQ, WANTR

      LOGICAL          TEST_QRDMD
!..... external subroutines (BLAS and LAPACK)
      EXTERNAL SAXPY,  SSYEV, SGEMM, SGEMV, SLACPY, SLASCL
      EXTERNAL SLARNV, SLATMR
!.....external subroutines DMD package, part 1
!     subroutines under test
      EXTERNAL SSYDMD, SSYDMDQ

!..... external functions (BLAS and LAPACK)
      EXTERNAL         SLAMCH, SLANGE, SNRM2
      REAL(KIND=WP) :: SLAMCH, SLANGE, SNRM2
      EXTERNAL ISAMAX
      INTEGER  ISAMAX
      EXTERNAL         LSAME
      LOGICAL          LSAME

      INTRINSIC ABS, INT, MIN, MAX
!............................................................

      ! The test is always in pairs : ( SSYDMD and SSYDMDQ )
      ! because the test includes comparing the results (in pairs).
!.....................................................................................
      TEST_QRDMD = .TRUE. ! This code by default performs tests on SSYDMDQ
                          ! Since the QR factorizations based algorithm is designed for
                          ! single trajectory data, only single trajectory tests will
                          ! be performed with xGEDMDQ.
      WANTQ = 'Q'
      WANTR = 'R'
!.................................................................................

      EPS = SLAMCH( 'P' )  ! machine precision SP

      ! Global counters of failures of some particular tests
      NFAIL      = 0
      NFAIL_REZ  = 0
      NFAIL_REZQ = 0
      NFAIL_Z_XV = 0
      NFAIL_F_QR = 0
      NFAIL_AU   = 0
      KDIFF      = 0
      NFAIL_SVDIFF = 0
      NFAIL_TOTAL  = 0
      NFAILQ_TOTAL = 0


      DO LLOOP = 1, 4

      WRITE(*,*) 'L Loop Index = ', LLOOP

      ! Set the dimensions of the problem ...
      WRITE(*,*) 'M = '
      READ(*,*) M
      WRITE(*,*) M
      ! ... and the number of snapshots.
      WRITE(*,*) 'N = '
      READ(*,*) N
      WRITE(*,*) N

      ! ... Test the dimensions
      IF ( ( MIN(M,N) == 0 ) .OR. ( M < N )  ) THEN
          WRITE(*,*) 'Bad dimensions. Required: M >= N > 0.'
          STOP
      END IF
!.............
      ! The seed inside the LLOOP so that each pass can be reproduced easily.

      ISEED(1) = 4
      ISEED(2) = 3
      ISEED(3) = 2
      ISEED(4) = 1

      LDA = M
      LDF = M
      LDX = MAX(M,N+1)
      LDY = MAX(M,N+1)
      LDW = N
      LDZ = M
      LDAU = MAX(M,N+1)
      LDS = N

      TMP_ZXW  = ZERO
      TMP_AU   = ZERO
      TMP_REZ  = ZERO
      TMP_REZQ = ZERO
      SVDIFF   = ZERO
      TMP_EX   = ZERO

      !
      ! Test the subroutines on real data snapshots. All
      ! computation is done in real arithmetic, even when
      ! Koopman eigenvalues and modes are real.
      !
      ! Allocate memory space
      ALLOCATE( A(LDA,M) )
      ALLOCATE( AC(LDA,M) )
      ALLOCATE( DA(M) )
      ALLOCATE( DL(M) )
      ALLOCATE( F(LDF,N+1) )
      ALLOCATE( F1(LDF,N+1) )
      ALLOCATE( F2(LDF,N+1) )
      ALLOCATE( X(LDX,N) )
      ALLOCATE( X0(LDX,N) )
      ALLOCATE( SINGVX(N) )
      ALLOCATE( SINGVQX(N) )
      ALLOCATE( Y(LDY,N+1) )
      ALLOCATE( Y0(LDY,N+1) )
      ALLOCATE( Y1(LDY,N+1) )
      ALLOCATE( Z(LDZ,N) )
      ALLOCATE( Z1(LDZ,N) )
      ALLOCATE( RES(N)  )
      ALLOCATE( RES1(N) )
      ALLOCATE( RESEX(N) )
      ALLOCATE( REIG(N) )
      ALLOCATE( IEIG(N) )
      ALLOCATE( REIGQ(N) )
      ALLOCATE( IEIGQ(N) )
      ALLOCATE( REIGA(M) )
      ALLOCATE( IEIGA(M) )
      ALLOCATE( VA(LDA,M) )
      ALLOCATE( LAMBDA(N,2) )
      ALLOCATE( LAMBDAQ(N,2) )
      ALLOCATE( EIGA(M,2) )
      ALLOCATE( W(LDW,N) )
      ALLOCATE( AU(LDAU,N) )
      ALLOCATE( S(LDS,LDS) )

      TOL  = M*EPS
      ! This mimics O(M*N)*EPS bound for accumulated roundoff error.
      ! The factor 10 is somewhat arbitrary.
      TOL2 = 10*M*N*EPS

!.............

      DO K_TRAJ = 1, 2
      !  Number of intial conditions in the simulation/trajectories (1 or 2)

      COND = 1.0D8
      DMAX = 1.0D2
      RSIGN = 'F'
      GRADE = 'N'
      MODEL = 6
      CONDL = 1.0D2
      MODER = 6
      CONDR = 1.0D2
      PIVTNG = 'N'

      ! Loop over all parameter MODE values for ZLATMR (+1,..,+6)
      DO MODE = 1, 6

      ALLOCATE( IWORK(2*M) )
      ALLOCATE(DR(N))
      CALL SLATMR( M, M, 'S', ISEED, 'S', DA, MODE, COND, &
                   DMAX, RSIGN, GRADE, DL, MODEL,  CONDL, &
                   DR, MODER, CONDR, PIVTNG, IWORK, M, M, &
                   ZERO, -ONE, 'N', A, LDA, IWORK(M+1), INFO )
      DEALLOCATE(IWORK)
      DEALLOCATE(DR)

           
      LWORK = MAX(1,3*M-1)
      ALLOCATE(WORK(LWORK))
      AC(1:M,1:M)  = A(1:M,1:M)  
      CALL SSYEV( 'V', 'U', M, AC, LDA, REIGA, WORK, LWORK, INFO ) ! LAPACK CALL 
      DEALLOCATE(WORK)
      
      TMP = ABS(REIGA(ISAMAX(M, REIGA, 1))) ! The spectral radius of A
      

      ! Scale A to have the desirable spectral radius.
      CALL SLASCL( 'G', 0, 0, TMP, ONE, M, M, A, M, INFO )
      CALL SLASCL( 'G', 0, 0, TMP, ONE, M, 1, REIGA, M, INFO )

      ! Compute the norm of A
      ANORM = SLANGE( 'F', N, N, A, M, WDUMMY )

      IF ( K_TRAJ == 2 ) THEN
          ! generate data with two inital conditions
      CALL SLARNV(2, ISEED, M, F1(1,1) )
      F1(1:M,1) = 1.0E-10*F1(1:M,1)
      DO i = 1, N/2
         CALL SGEMV( 'N', M, M, ONE, A, M, F1(1,i), 1, ZERO, &
              F1(1,i+1), 1 )
      END DO
      X0(1:M,1:N/2) = F1(1:M,1:N/2)
      Y0(1:M,1:N/2) = F1(1:M,2:N/2+1)

      CALL SLARNV(2, ISEED, M, F1(1,1) )
      DO i = 1, N-N/2
         CALL SGEMV( 'N', M, M, ONE, A, M, F1(1,i), 1, ZERO, &
              F1(1,i+1), 1 )
      END DO
      X0(1:M,N/2+1:N) = F1(1:M,1:N-N/2)
      Y0(1:M,N/2+1:N) = F1(1:M,2:N-N/2+1)
      ELSE
          ! single trajectory
      CALL SLARNV(2, ISEED, M, F(1,1) )
      DO i = 1, N
         CALL SGEMV( 'N', M, M, ONE, A, M, F(1,i), 1, ZERO, &
              F(1,i+1), 1 )
      END DO
      X0(1:M,1:N) = F(1:M,1:N)
      Y0(1:M,1:N) = F(1:M,2:N+1)
      END IF

      XNORM = SLANGE( 'F', M, N, X0, LDX, WDUMMY )
      YNORM = SLANGE( 'F', M, N, Y0, LDX, WDUMMY )
!............................................................

      DO iJOBZ = 1, 4

          SELECT CASE ( iJOBZ )
          CASE(1)
              JOBZ   = 'V' ! Ritz vectors will be computed
              RESIDS = 'R' ! Residuals will be computed
          CASE(2)
              JOBZ   = 'V'
              RESIDS = 'N'
          CASE(3)
              JOBZ   = 'F' ! Ritz vectors in factored form
              RESIDS = 'N'
          CASE(4)
              JOBZ   = 'N'
              RESIDS = 'N'
          END SELECT

      DO iJOBREF = 1, 3

          SELECT CASE ( iJOBREF )
          CASE(1)
              JOBREF = 'R' ! Data for refined Ritz vectors
          CASE(2)
              JOBREF = 'E' ! Exact DMD vectors
          CASE(3)
              JOBREF = 'N'
          END SELECT

      DO iSCALE = 1, 4

          SELECT CASE ( iSCALE )
          CASE(1)
              SCALE = 'S' ! X data normalized
          CASE(2)
              SCALE = 'C' ! X normalized, consist. check
          CASE(3)
              SCALE = 'Y' ! Y data normalized
          CASE(4)
              SCALE = 'N'
          END SELECT

      DO iNRNK = -1, -2, -1
          ! Two truncation strategies. The "-2" case for R&D
          ! purposes only - it uses possibly low accuracy small
          ! singular values, in which case the formulas used in
          ! the DMD are highly sensitive.
          NRNK   = iNRNK

      DO iWHTSVD = 1, 4
          ! Check all four options to compute the POD basis
          ! via the SVD.
          WHTSVD   = iWHTSVD
          
      DO iWHTEIG = 1, 2
          ! Check both symmetric eigensolvers in LAPACK
          WHTEIG = iWHTEIG 
      
      DO iWHTSYM = 1, 2 
          ! Check both symmetrizers of the Rayleigh quotient
          WHTSYM = iWHTSYM    

      DO LWMINOPT = 1, 2
          ! Workspace query for the minimal (1) and for the optimal
          ! (2) workspace lengths determined by workspace query.

       X(1:M,1:N) = X0(1:M,1:N)
       Y(1:M,1:N) = Y0(1:M,1:N)

       ! SSYDMD: Workspace query and workspace allocation
       
       CALL SSYDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD, WHTSYM, WHTEIG, M, &
            N, X, LDX, Y, LDY, NRNK, TOL, K, REIG, Z, &
            LDZ, RES, AU, LDAU, W, LDW, S, LDS, WDUMMY, -1,&
            IDUMMY, -1, INFO )


       LIWORK = IDUMMY(1)
       ALLOCATE( IWORK(LIWORK) )
       LWORK = INT(WDUMMY(LWMINOPT))
       ALLOCATE( WORK(LWORK) )

       ! SSYDMD test: CALL SSYDMD
       
       CALL SSYDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD, WHTSYM, WHTEIG, M, &
            N, X, LDX, Y, LDY, NRNK, TOL, K, REIG, Z, &
            LDZ, RES, AU, LDAU, W, LDW, S, LDS, WORK, LWORK,&
            IWORK, LIWORK, INFO )

       SINGVX(1:N) = WORK(1:N)

       !...... SSYDMD check point
       IF ( LSAME(JOBZ,'V')  ) THEN
          ! Check that Z = X*W, on return from SSYDMD
          ! This checks that the returned aigenvectors in Z are
          ! the product of the SVD'POD basis returned in X
          ! and the eigenvectors of the rayleigh quotient
          ! returned in W
          CALL SGEMM( 'N', 'N', M, K, K, ONE, X, LDX, W, LDW, &
                      ZERO, Z1, LDZ )
          TMP = ZERO
          DO i = 1, K
             CALL SAXPY( M, -ONE, Z(1,i), 1, Z1(1,i), 1)
             TMP = MAX(TMP, SNRM2( M, Z1(1,i), 1 ) )
          END DO
          TMP_ZXW = MAX(TMP_ZXW, TMP )

          IF ( TMP_ZXW > 10*M*EPS ) THEN
              NFAIL_Z_XV = NFAIL_Z_XV + 1
          END IF

       END IF

       !...... SSYDMD check point
       IF ( LSAME(JOBREF,'R') ) THEN
           ! The matrix A*U is returned for computing refined Ritz vectors.
           ! Check that A*U is computed correctly using the formula
           ! A*U = Y * V * inv(SIGMA). This depends on the
           ! accuracy in the computed singular values and vectors of X.
           ! See the paper for an error analysis.
           ! Note that the left singular vectors of the input matrix X
           ! are returned in the array X.
           CALL SGEMM( 'N', 'N', M, K, M, ONE, A, LDA, X, LDX, &
                      ZERO, Z1, LDZ )
           TMP = ZERO
           DO i = 1, K
              CALL SAXPY( M, -ONE, AU(1,i), 1, Z1(1,i), 1)
              TMP = MAX( TMP, SNRM2( M, Z1(1,i),1 ) * &
                       SINGVX(K)/(ANORM*SINGVX(1)) )
           END DO
           TMP_AU = MAX( TMP_AU, TMP )

           IF ( TMP > TOL2 ) THEN
               NFAIL_AU = NFAIL_AU + 1
           END IF

       ELSEIF ( LSAME(JOBREF,'E') ) THEN
       ! The unscaled vectors of the Exact DMD are computed.
       ! This option is included for the sake of completeness,
       ! for users who prefer the Exact DMD vectors. The
       ! returned vectors are in the real form, in the same way
       ! as the Ritz vectors. Here we just save the vectors
       ! and test them separately using a Matlab script.

       CALL SGEMM( 'N', 'N', M, K, M, ONE, A, LDA, AU, LDAU, ZERO, Y1, LDY )
       i=1
       DO WHILE ( i <= K )
        ! have a real eigenvalue with real eigenvector
        CALL SAXPY( M, -REIG(i), AU(1,i), 1, Y1(1,i), 1 )
        RESEX(i) = SNRM2( M, Y1(1,i), 1) / SNRM2(M,AU(1,i),1)
        i = i + 1
       END DO

       END IF

      !...... SSYDMD check point
      IF ( LSAME(RESIDS, 'R') ) THEN
          ! Compare the residuals returned by SSYDMD with the
          ! explicitly computed residuals using the matrix A.
          ! Compute explicitly Y1 = A*Z
          CALL SGEMM( 'N', 'N', M, K, M, ONE, A, LDA, Z, LDZ, ZERO, Y1, LDY )
          ! ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
          ! of the invariant subspaces that correspond to complex conjugate
          ! pairs of eigencalues. (See the description of Z in SSYDMD,)
          i = 1
          DO WHILE ( i <= K )
                ! have a real eigenvalue with real eigenvector
                CALL SAXPY( M, -REIG(i), Z(1,i), 1, Y1(1,i), 1 )
                RES1(i) = SNRM2( M, Y1(1,i), 1)
                i = i + 1
          END DO
          TMP = ZERO
          DO i = 1, K
          TMP = MAX( TMP, ABS(RES(i) - RES1(i)) * &
                    SINGVX(K)/(ANORM*SINGVX(1)) )
          END DO
          TMP_REZ = MAX( TMP_REZ, TMP )

          IF ( TMP > TOL2 ) THEN
              NFAIL_REZ = NFAIL_REZ + 1
          END IF

         IF ( LSAME(JOBREF,'E') ) THEN
            TMP = ZERO
          DO i = 1, K
          TMP = MAX( TMP, ABS(RES1(i) - RESEX(i))/(RES1(i)+RESEX(i)) )
          END DO
          TMP_EX = MAX(TMP_EX,TMP)
         END IF

      END IF

      ! ... store the results for inspection
!!      DO i = 1, K
!!          LAMBDA(i,1) = REIG(i)
!!          LAMBDA(i,2) = IEIG(i)
!!      END DO

      DEALLOCATE(IWORK)
      DEALLOCATE(WORK)

      !======================================================================
      !     Now test the SSYDMDQ, if requested.
      !======================================================================
      IF ( TEST_QRDMD .AND. (K_TRAJ == 1) ) THEN
    
          F1(1:M,1:N+1) = F(1:M,1:N+1)

          ! SSYDMDQ test: Workspace query and workspace allocation
          
          CALL SSYDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, WHTSVD, WHTSYM, WHTEIG, M, N+1, &
                    F1, LDF, X, LDX, Y, LDY, NRNK, TOL, KQ, REIGQ,   &
                    Z, LDZ, RES, AU, LDAU, W, LDW, S, LDS, WDUMMY, &
                    -1, IDUMMY, -1, INFO )
         
          LIWORK = IDUMMY(1)
          ALLOCATE( IWORK(LIWORK) )
          LWORK = INT(WDUMMY(LWMINOPT))
          ALLOCATE(WORK(LWORK))

          ! SSYDMDQ test: CALL SSYDMDQ
          
          CALL SSYDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, WHTSVD, WHTSYM, WHTEIG, M, N+1, &
                    F1, LDF, X, LDX, Y, LDY, NRNK, TOL, KQ, REIGQ,   &
                    Z, LDZ, RES, AU, LDAU, W, LDW, S, LDS, WORK, &
                    LWORK, IWORK, LIWORK, INFO )
         

          SINGVQX(1:KQ) = WORK(MIN(M,N+1)+1: MIN(M,N+1)+KQ)

          !..... SSYDMDQ check point
          IF ( KQ /= K ) THEN
             KDIFF = KDIFF+1
          END IF

          TMP = ZERO
          DO i = 1, MIN(K, KQ)
             TMP = MAX(TMP, ABS(SINGVX(i)-SINGVQX(i)) / &
                                   SINGVX(1) )
          END DO
          SVDIFF = MAX( SVDIFF, TMP )
          IF ( TMP > M*N*EPS ) THEN
             NFAIL_SVDIFF = NFAIL_SVDIFF + 1
          END IF

          !..... SSYDMDQ check point
          IF ( LSAME(WANTQ,'Q') .AND. LSAME(WANTR,'R') ) THEN
             ! Check that the QR factors are computed and returned
             ! as requested. The residual ||F-Q*R||_F / ||F||_F
             ! is compared to M*N*EPS.
             F2 = F
             CALL SGEMM( 'N', 'N', M, N+1, MIN(M,N+1), -ONE, F1, &
                         LDF, Y, LDY, ONE, F2, LDF )
             TMP_FQR = SLANGE( 'F', M, N+1, F2, LDF, WORK ) / &
                   SLANGE( 'F', M, N+1, F,  LDF, WORK )
             IF ( TMP_FQR > TOL2 ) THEN
                 NFAIL_F_QR = NFAIL_F_QR + 1
             END IF
          END IF

          !..... SSYDMDQ checkpoint
          IF ( LSAME(RESIDS, 'R') ) THEN
              ! Compare the residuals returned by SSYDMDQ with the
              ! explicitly computed residuals using the matrix A.
              ! Compute explicitly Y1 = A*Z
              CALL SGEMM( 'N', 'N', M, KQ, M, ONE, A, LDA, Z, LDZ, ZERO, Y1, LDY )
              ! ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
              ! of the invariant subspaces that correspond to complex conjugate
              ! pairs of eigencalues. (See the description of Z in SSYDMDQ)
              i = 1
              DO WHILE ( i <= KQ )
                    ! have a real eigenvalue with real eigenvector
                    CALL SAXPY( M, -REIGQ(i), Z(1,i), 1, Y1(1,i), 1 )
                    ! Y(1:M,i) = Y(1:M,i) - REIG(i)*Z(1:M,i)
                    RES1(i) = SNRM2( M, Y1(1,i), 1)
                    i = i + 1
              END DO
              TMP = ZERO
              DO i = 1, KQ
              TMP = MAX( TMP, ABS(RES(i) - RES1(i)) * &
                  SINGVQX(K)/(ANORM*SINGVQX(1)) )
              END DO
              TMP_REZQ = MAX( TMP_REZQ, TMP )
              IF ( TMP > TOL2 ) THEN
                  NFAIL_REZQ = NFAIL_REZQ + 1
              END IF

          END IF

          DO i = 1, KQ
              LAMBDAQ(i,1) = REIGQ(i)
              LAMBDAQ(i,2) = IEIGQ(i)
          END DO

      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)
      END IF            ! TEST_QRDMD
!======================================================================

      END DO ! LWMINOPT
      !write(*,*) 'LWMINOPT loop completed'
      END DO
      END DO 
      END DO ! WHTSVD LOOP
      !write(*,*) 'WHTSVD loop completed'
      END DO ! NRNK LOOP
      !write(*,*) 'NRNK loop completed'
      END DO ! SCALE LOOP
      !write(*,*) 'SCALE loop completed'
      END DO ! JOBF LOOP
      !write(*,*) 'JOBREF loop completed'
      END DO ! JOBZ LOOP
      !write(*,*) 'JOBZ loop completed'

      END DO ! MODE -6:6
      !write(*,*) 'MODE loop completed'
      END DO ! 1 or 2 trajectories
      !write(*,*) 'trajectories  loop completed'

      DEALLOCATE(A)
      DEALLOCATE(AC)
      DEALLOCATE(DA)
      DEALLOCATE(DL)
      DEALLOCATE(F)
      DEALLOCATE(F1)
      DEALLOCATE(F2)
      DEALLOCATE(X)
      DEALLOCATE(X0)
      DEALLOCATE(SINGVX)
      DEALLOCATE(SINGVQX)
      DEALLOCATE(Y)
      DEALLOCATE(Y0)
      DEALLOCATE(Y1)
      DEALLOCATE(Z)
      DEALLOCATE(Z1)
      DEALLOCATE(RES)
      DEALLOCATE(RES1)
      DEALLOCATE(RESEX)
      DEALLOCATE(REIG)
      DEALLOCATE(IEIG)
      DEALLOCATE(REIGQ)
      DEALLOCATE(IEIGQ)
      DEALLOCATE(REIGA)
      DEALLOCATE(IEIGA)
      DEALLOCATE(VA)
      DEALLOCATE(LAMBDA)
      DEALLOCATE(LAMBDAQ)
      DEALLOCATE(EIGA)
      DEALLOCATE(W)
      DEALLOCATE(AU)
      DEALLOCATE(S)

!............................................................
      !     Generate random M-by-M matrix A. Use DLATMR from
      END DO ! LLOOP


      WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
      WRITE(*,*) ' Test summary for SSYDMD :'
      WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
      WRITE(*,*)
      IF ( NFAIL_Z_XV == 0 ) THEN
          WRITE(*,*) '>>>> Z - U*V test PASSED.'
      ELSE
          WRITE(*,*) 'Z - U*V test FAILED ', NFAIL_Z_XV, ' time(s)'
          WRITE(*,*) 'Max error ||Z-U*V||_F was ', TMP_ZXW
          NFAIL_TOTAL = NFAIL_TOTAL + NFAIL_Z_XV
      END IF
      IF ( NFAIL_AU == 0 ) THEN
          WRITE(*,*) '>>>> A*U test PASSED. '
      ELSE
          WRITE(*,*) 'A*U test FAILED ', NFAIL_AU, ' time(s)'
          WRITE(*,*) 'Max A*U test adjusted error measure was ', TMP_AU
          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAIL_TOTAL = NFAIL_TOTAL + NFAIL_AU
      END IF

      IF ( NFAIL_REZ == 0 ) THEN
          WRITE(*,*) '>>>> Rezidual computation test PASSED.'
      ELSE
          WRITE(*,*) 'Rezidual computation test FAILED ', NFAIL_REZ, 'time(s)'
          WRITE(*,*) 'Max residual computing test adjusted error measure was ', TMP_REZ
          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAIL_TOTAL = NFAIL_TOTAL + NFAIL_REZ
      END IF

      IF ( NFAIL_TOTAL == 0 ) THEN
          WRITE(*,*) '>>>> SSYDMD :: ALL TESTS PASSED.'
      ELSE
          WRITE(*,*) NFAIL_TOTAL, 'FAILURES!'
          WRITE(*,*) '>>>>>>>>>>>>>> SSYDMD :: TESTS FAILED. CHECK THE IMPLEMENTATION.'
      END IF

      IF ( TEST_QRDMD ) THEN
      WRITE(*,*)
      WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
      WRITE(*,*) ' Test summary for SSYDMDQ :'
      WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
      WRITE(*,*)

      IF ( NFAIL_SVDIFF == 0 ) THEN
          WRITE(*,*) '>>>> SSYDMD and SSYDMDQ computed singular &
              &values test PASSED.'
      ELSE
          WRITE(*,*) 'SSYDMD and SSYDMDQ discrepancies in &
              &the singular values unacceptable ', &
              NFAIL_SVDIFF, ' times. Test FAILED.'
          WRITE(*,*) 'The maximal discrepancy in the singular values (relative to the norm) was ', SVDIFF
          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAILQ_TOTAL = NFAILQ_TOTAL + NFAIL_SVDIFF
      END IF

      IF ( NFAIL_F_QR == 0 ) THEN
          WRITE(*,*) '>>>> F - Q*R test PASSED.'
      ELSE
          WRITE(*,*) 'F - Q*R test FAILED ', NFAIL_F_QR, ' time(s)'
          WRITE(*,*) 'The largest relative residual was ', TMP_FQR
          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAILQ_TOTAL = NFAILQ_TOTAL + NFAIL_F_QR
      END IF

      IF ( NFAIL_REZQ == 0 ) THEN
          WRITE(*,*) '>>>> Rezidual computation test PASSED.'
      ELSE
          WRITE(*,*) 'Rezidual computation test FAILED ', NFAIL_REZQ, 'time(s)'
          WRITE(*,*) 'Max residual computing test adjusted error measure was ', TMP_REZQ
          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAILQ_TOTAL = NFAILQ_TOTAL + NFAIL_REZQ
      END IF

      IF ( NFAILQ_TOTAL == 0 ) THEN
          WRITE(*,*) '>>>>>>> SSYDMDQ :: ALL TESTS PASSED.'
      ELSE
         WRITE(*,*) NFAILQ_TOTAL, 'FAILURES!'
         WRITE(*,*) '>>>>>>> SSYDMDQ :: TESTS FAILED. CHECK THE IMPLEMENTATION.'
      END IF

      END IF

      WRITE(*,*)
      WRITE(*,*) 'Test completed.'
      STOP
      END
