! This is a test program for checking the implementations of
! the implementations of the following subroutines
!
! DGEDMD,  SGEDMD, ZGEDMD, CGEDMD  for computation of the
!                        Dynamic Mode Decomposition (DMD)
! DGEDMDQ, SGEDMDQ, ZGEDMDQ, CGEDMDQ for computation of a
!                  QR factorization based compressed DMD
!
! Developed and supported by:
! ===========================
! Developed and coded by Zlatko Drmac, Faculty of Science,
! University of Zagreb;  drmac@math.hr
! In cooperation with
! AIMdyn Inc., Santa Barbara, CA.
!   ========================================================
!   How to run the code (compiler, link info)
!   ========================================================
!   Compile as FORTRAN 90 (or later) and link with BLAS and
!   LAPACK libraries.
!   NOTE: The code is developed and tested on top of the
!   Intel MKL library (versions 2022.0.3 and 2022.2.0),
!   using the Intel Fortran compiler.
!
!   For developers of the C++ implementation
!   ========================================================
!   See the LAPACK++ and Template Numerical Toolkit (TNT)
!
!   Note on a development of the GPU HP implementation
!   ========================================================
!   Work in progress. See CUDA, MAGMA, SLATE.
!   NOTE: The four SVD subroutines used in this code are
!   included as a part of R&D and for the completeness.
!   This was also an opportunity to test those SVD codes.
!   If the scaling option is used all four are essentially
!   equally good. For implementations on HP platforms,
!   one can use whichever SVD is available.
!............................................................
!   NOTE:
!   When using the Intel MKL 2022.0.3 the subroutine xGESVDQ
!   (optionally used in xGEDMD) may cause access violation
!   error for x = S, D, C, Z, but only if called with the
!   work space query. (At least in our Windows 10 MSVS 2019.)
!   The problem can be mitigated by downloading the source
!   code of xGESVDQ from the LAPACK repository and use it
!   localy instead of the one in the MKL. This seems to
!   indicate that the problem is indeed in the MKL.
!   This problem did not appear whith Intel MKL 2022.2.0.
!
!   NOTE:
!   xGESDD seems to have a problem with workspace. In some
!   cases the length of the optimal workspace is returned
!   smaller than the minimal workspace, as specified in the
!   code. As a precaution, all optimal workspaces are
!   set as MAX(minimal, optimal).
!   Latest implementations of complex xGESDD have different
!   length of the real worksapce. We use max value over
!   two versions.
!............................................................
!............................................................
!
      PROGRAM DMD_TEST
!
      use iso_fortran_env, only: real32, real64
      IMPLICIT NONE
      integer, parameter :: WP = real64
      integer, parameter :: SP = real32

!............................................................
      REAL(KIND=WP), PARAMETER ::  ONE = 1.0_WP
      REAL(KIND=WP), PARAMETER :: ZERO = 0.0_WP
      COMPLEX(KIND=WP), PARAMETER ::  ZONE = ( 1.0_WP, 0.0_WP )
      COMPLEX(KIND=WP), PARAMETER :: ZZERO = ( 0.0_WP, 0.0_WP )
!............................................................
      REAL(KIND=WP), ALLOCATABLE, DIMENSION(:,:) ::          &
                     A, AC, EIGA, LAMBDA, LAMBDAQ, F, F1, F2,&
                     Z, Z1, S, AU, W, VA, X, X0, Y, Y0, Y1
      REAL(KIND=WP), ALLOCATABLE, DIMENSION(:)   ::          &
                     DA, DL, DR, REIG, REIGA, REIGQ, IEIG,   &
                     IEIGA, IEIGQ,  RES, RES1, RESEX, SINGVX,&
                     SINGVXsp, SINGVQX, WORK
      INTEGER      , ALLOCATABLE, DIMENSION(:)   ::   IWORK
      REAL(KIND=WP) :: AB(2,2),   WDUMMY(2)
      INTEGER       :: IDUMMY(2), ISEED(4), RJOBDATA(8),     &
                       ZJOBDATA(8)
      REAL(KIND=WP) :: ANORM, COND, CONDL, CONDR, DMAX, EPS, &
                       TOL, TOL2, SVDIFF, TMP, TMP_AU,       &
                       TMP_FQR, TMP_REZ, TMP_REZQ,  TMP_ZXW, &
                       TMP_EX, XNORM, YNORM
!..... for real data, single precision tests
      REAL(KIND=SP), ALLOCATABLE, DIMENSION(:,:) ::          &
                     Asp, AUsp, Fsp, F1sp, LAMBDAsp,         &
                     LAMBDAQsp, Ssp, Wsp, Xsp, Ysp,  Zsp
      REAL(KIND=SP), ALLOCATABLE, DIMENSION(:)   ::          &
                     IEIGsp, IEIGQsp, REIGsp, REIGQsp, RESsp,&
                     SWORK
      REAL(KIND=SP) :: WDUMMYsp(2)
      REAL(KIND=SP) :: EPSsp, TOLsp
      INTEGER :: Ksp, KQsp
!............................................................
      COMPLEX(KIND=WP) ZMAX
      INTEGER LZWORK
      COMPLEX(KIND=WP), ALLOCATABLE, DIMENSION(:,:) ::       &
                        ZA, ZAC, ZAU, ZF, ZF0, ZS, ZVA, ZW,  &
                        ZX, ZY,  ZZ
      COMPLEX(KIND=WP), ALLOCATABLE, DIMENSION(:)   ::       &
                        ZDA, ZDR, ZDL, ZEIGS, ZEIGSA, ZWORK
      COMPLEX(KIND=WP) ZDUMMY(22)
!..... for complex data, single precision test
      COMPLEX(KIND=SP), ALLOCATABLE, DIMENSION(:,:) ::       &
                 ZAUsp, ZFsp, ZSsp, ZWsp, ZXsp, ZYsp, ZZsp
      COMPLEX(KIND=SP), ALLOCATABLE, DIMENSION(:)   ::       &
                        ZEIGSsp, ZWORKsp
      COMPLEX(KIND=SP) ZDUMMYsp(22)
!............................................................
      INTEGER :: K, KQ, LDF, LDS, LDA, LDAU, LDW, LDX, LDY,  &
                 LDZ, LIWORK, LWORK, M, N, NRNK, NRNKsp
      INTEGER :: i, iJOBREF, iJOBZ, iSCALE, INFO, KDIFF,     &
                 NFAIL, NFAIL_AU, NFAIL_F_QR, NFAIL_REZ,     &
                 NFAIL_REZQ, NFAIL_SVDIFF, NFAIL_TOTAL,      &
                 NFAIL_Z_XV, MODE, MODEL, MODER, WHTSVD,     &
                 WHTSVDsp
      INTEGER    iSCALEd, iJOBZd, iNRNK, iWHTSVD,  iTOL,     &
                 iREF, LWMINOPT
      CHARACTER(LEN=1) GRADE, JOBREF, JOBZ, PIVTNG, RSIGN,   &
                       SCALE, RESIDS, WANTQ, WANTR
      LOGICAL          KEYBIN, MANUAL_TEST, TEST_QRDMD,      &
                       TWO_TRAJ, TEST_SP, TESTQ_SP,          &
                       VISUAL_CHECK
      LOGICAL    REAL_DATA
!..... external subroutines (BLAS and LAPACK)
      EXTERNAL DAXPY,  DGEEV, DGEMM, DGEMV, DLACPY, DLASCL
      EXTERNAL DLARNV, DLATMR
      EXTERNAL ZGEEV,  ZGEMV, ZLASCL
      EXTERNAL ZLARNV, ZLATMR
!.....external subroutines DMD package, part 1
!     subroutines under test
      EXTERNAL DGEDMD, DGEDMDQ, SGEDMD, SGEDMDQ
      EXTERNAL ZGEDMD, ZGEDMDQ, CGEDMD, CGEDMDQ
!.....auxiliary, for test only
      EXTERNAL MATRIX_SAVE, MATRIX_SAVE_SP, &
               MATRIX_SAVE_ZRE, MATRIX_SAVE_ZIM, &
               MATRIX_SAVE_CRE, MATRIX_SAVE_CIM

!..... external functions (BLAS and LAPACK)
      EXTERNAL         DLAMCH, DLANGE, DNRM2, SLAMCH
      REAL(KIND=WP) :: DLAMCH, DLANGE, DNRM2
      REAL(KIND=SP) :: SLAMCH
      EXTERNAL IZAMAX
      INTEGER  IZAMAX
      EXTERNAL         LSAME
      LOGICAL          LSAME

      INTRINSIC ABS, INT, MIN, MAX
!............................................................

      EPS   = DLAMCH( 'P' )  ! machine precision DP
!      EPSsp = SLAMCH( 'P' )  ! machine precision SP

!      WRITE(*,fmt='(a)',advance='no')    &
!          'Test REAL subroutines (1), ', &
!          ' or COMPLEX versions (2) >> '
!      READ(*,*) K
!      IF ( K == 1 ) THEN
         REAL_DATA = .TRUE.
!      ELSEIF ( K == 2 ) THEN
!         REAL_DATA = .FALSE.
!      ELSE
!         WRITE(*,*) 'Bad input. Exit.'
!         STOP
!      END IF

!      WRITE(*,fmt='(a)',advance='no') &
!          'Manual individual case studies (1), ', &
!          'or automatic loop over job parameters (2) >> '
!      READ(*,*) K
!      IF ( K == 1 ) THEN
!          KEYBIN       = .TRUE.
!          VISUAL_CHECK = .TRUE.
!          MANUAL_TEST  = .TRUE.
!      ELSEIF ( K == 2 ) THEN
          KEYBIN       = .FALSE.
          VISUAL_CHECK = .FALSE.
          MANUAL_TEST  = .FALSE.
!      ELSE
!          WRITE(*,*) 'Bad input. Exit.'
!          STOP
!      END IF

!      TEST_SP = .TRUE.
!      TEST_SP = TEST_SP .AND. MANUAL_TEST  ! Single precision code is tested
!                                           ! only as manual case studies.
!      TESTQ_SP = .TRUE.
!      TESTQ_SP = TESTQ_SP .AND. MANUAL_TEST

      WRITE(*,fmt='(a)',advance='no') &
          'Number of intial cond./trajectories (1 or 2) >> '
      READ(*,*) K
      IF ( K == 1) THEN
          TWO_TRAJ = .FALSE.
      ELSE IF ( K == 2 ) THEN
          TWO_TRAJ = .TRUE.
      ELSE
          WRITE(*,*) 'Bad input.'
          STOP
      END IF
      IF ( TWO_TRAJ ) THEN
          TEST_QRDMD = .FALSE.
      ELSE                    ! This version (ZGEDMDQ) is for
          TEST_QRDMD = .TRUE. ! one trajectory stream of data.
      END IF

      ! Set the dimensions of the problem ...
      WRITE(*,fmt='(a)',advance='no') &
          'State space dimension  M =  '
      READ(*,*) M
!     ... and the number of snapshots.
      WRITE(*,fmt='(a)',advance='no') &
          'The number of data snapshot pairs  N =  '
      READ(*,*) N
      IF ( ( MIN(M,N) == 0 ) .OR. ( M < N )  ) THEN
          WRITE(*,*) 'Bad dimensions. Required: M >= N > 0.'
          STOP
      END IF

      ISEED(1) = 4
      ISEED(2) = 3
      ISEED(3) = 2
      ISEED(4) = 1

      WRITE(*,fmt='(a)',advance='no') &
          'Parameter MODE for DLATMR (+-1,..,+-6) = '
      READ(*,*) MODE
!..............................................................

      LDA = M
      LDF = M
      LDX = MAX(M,N+1)
      LDY = MAX(M,N+1)
      LDW = N
      LDZ = M
      LDAU = MAX(M,N+1)
      LDS = N

      IF ( REAL_DATA ) THEN
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
      ALLOCATE( Y1(M,N+1) )
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
      ALLOCATE( S(N,N) )
      !IF ( TEST_SP ) THEN
      !    ALLOCATE( Fsp(LDF,N+1) )
      !    ALLOCATE( F1sp(LDF,N+1) )
      !    ALLOCATE( Asp(LDA,M))
      !    ALLOCATE( Xsp(LDX,N) )
      !    ALLOCATE( Ysp(LDY,N+1) )
      !    ALLOCATE( SINGVXsp(N) )
      !    ALLOCATE( REIGsp(N) )
      !    ALLOCATE( IEIGsp(N) )
      !    ALLOCATE( REIGQsp(N) )
      !    ALLOCATE( IEIGQsp(N) )
      !    ALLOCATE( Zsp(LDZ,N) )
      !    ALLOCATE( RESsp(N)  )
      !    ALLOCATE( AUsp(LDAU,N))
      !    ALLOCATE( Ssp(N,N))
      !    ALLOCATE( Wsp(LDW,N) )
      !    ALLOCATE( LAMBDAsp(N,2) )
      !    ALLOCATE( LAMBDAQsp(N,2) )
      !END IF

!............................................................
!     Generate random M-by-M matrix A. Use DLATMR from
!     LAPACK. The test uses a single test matrix and then
!     runs the subroutines under test with different
!     options set by the input parameters.

      COND = 1.0D8
      DMAX = 1.0D2
      RSIGN = 'F'
      GRADE = 'N'
      MODEL = 6
      CONDL = 1.0D2
      MODER = 6
      CONDR = 1.0D2
      PIVTNG = 'N'
      ALLOCATE( IWORK(2*M) )
      ALLOCATE(DR(N))
      CALL DLATMR( M, M, 'S', ISEED, 'N', DA, MODE, COND, &
                   DMAX, RSIGN, GRADE, DL, MODEL,  CONDL, &
                   DR, MODER, CONDR, PIVTNG, IWORK, M, M, &
                   ZERO, -ONE, 'N', A, LDA, IWORK(M+1), INFO )
      DEALLOCATE(IWORK)
      DEALLOCATE(DR)
      ! Store reference eigenvalues and eigenvectors
!      WRITE(*,'(A)', advance='no') &
!      'Compute and store the eigenvalues &
!                 &and eigenvectors of A (DGEEV) . ..'
      LWORK = 4*M+1
      ALLOCATE(WORK(LWORK))
      AC  = A
      Asp = A
      CALL DGEEV( 'N','V', M, AC, M, REIGA, IEIGA, VA, M, &
                  VA, M, WORK, LWORK, WORK, INFO ) ! LAPACK CALL
      DEALLOCATE(WORK)
      TMP = ZERO
      DO i = 1, M
       EIGA(i,1) = REIGA(i)
       EIGA(i,2) = IEIGA(i)
       TMP = MAX( TMP, SQRT(REIGA(i)**2+IEIGA(i)**2))
      END DO
      !TMP = TMP*0.7
      ! Scale A to have the desirable spectral radius.
      CALL DLASCL( 'G', 0, 0, TMP, ONE, M, M, A, M, INFO )
      CALL DLASCL( 'G', 0, 0, TMP, ONE, M, 2, EIGA, M, INFO )

      ! Compute the norm of A
      ANORM = DLANGE( 'F', N, N, A, M, WDUMMY )


      ! Save the eigenvalues and eigenvectors of A for
      ! a post-factum analysis
!      CALL matrix_save( 9,  'eigsA.true',  M, 2, EIGA, M   )
!      CALL matrix_save( 9,  'evecsA.true', M, M, VA,   M   )
!      WRITE(*,*) ' ... done!'
      ! Generate random M-by-1 vector and use it for the
      ! Krylov sequence of length N+1.

!      WRITE(*,*) 'Generate trajectories with random initial states'
      IF ( TWO_TRAJ ) THEN
          ! generate data with two inital conditions
      CALL DLARNV(2, ISEED, M, F1(1,1) )
      F1(1:M,1) = 1.0E-10*F1(1:M,1)
      DO i = 1, N/2
         CALL DGEMV( 'N', M, M, ONE, A, M, F1(1,i), 1, ZERO, &
              F1(1,i+1), 1 )
      END DO
      X0(1:M,1:N/2) = F1(1:M,1:N/2)
      Y0(1:M,1:N/2) = F1(1:M,2:N/2+1)

      CALL DLARNV(2, ISEED, M, F1(1,1) )
      DO i = 1, N-N/2
         CALL DGEMV( 'N', M, M, ONE, A, M, F1(1,i), 1, ZERO, &
              F1(1,i+1), 1 )
      END DO
      X0(1:M,N/2+1:N) = F1(1:M,1:N-N/2)
      Y0(1:M,N/2+1:N) = F1(1:M,2:N-N/2+1)
      ELSE
          ! single trajectory
      CALL DLARNV(2, ISEED, M, F(1,1) )
      F(1:M,1) = 1.0D37*F(1:M,1)
      DO i = 1, N
         CALL DGEMV( 'N', M, M, ONE, A, M, F(1,i), 1, ZERO, &
              F(1,i+1), 1 )
      END DO
      X0(1:M,1:N) = F(1:M,1:N)
      Y0(1:M,1:N) = F(1:M,2:N+1)
      END IF

      XNORM = DLANGE( 'F', M, N, X0, LDX, WDUMMY )
      YNORM = DLANGE( 'F', M, N, Y0, LDX, WDUMMY )
!      WRITE(*,*) 'The norm of X = ', XNORM
!      WRITE(*,*) 'The norm of Y = ', YNORM
!
!      WRITE(*,*) '... random sequence of snapshots generated'
!      CALL matrix_save( 9, 'A.data', M, M, A, M )
!      CALL matrix_save( 9, 'X.data', M, N, X0, LDX )
!      CALL matrix_save( 9, 'Y.data', M, N, Y0, LDY )
!............................................................

      TOL   = M*EPS
!      TOLsp = M*EPSsp

      iTOL = 1

      NFAIL      = 0
      NFAIL_REZ  = 0
      NFAIL_REZQ = 0
      NFAIL_Z_XV = 0
      NFAIL_F_QR = 0
      NFAIL_AU   = 0
      KDIFF      = 0
      NFAIL_SVDIFF = 0
      NFAIL_TOTAL  = 0
      TMP_ZXW  = ZERO
      TMP_AU   = ZERO
      TMP_REZ  = ZERO
      TMP_REZQ = ZERO
      SVDIFF   = ZERO
      TMP_EX   = ZERO
      TOL2     = 10*M*N*EPS ! This mimics O(M*N)*EPS bound &
                            ! for accumulated roundoff error
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
         NRNKsp = iNRNK

      DO iWHTSVD = 1, 4
         ! Check all four options to compute the POD basis
         ! via the SVD.
         WHTSVD   = iWHTSVD
         WHTSVDsp = iWHTSVD

      DO LWMINOPT = 1, 2
         ! Workspace query for the minimal (1) and for the optimal
         ! (2) workspace lengths determined by workspace query.
!...... manual test section
!       IF ( MANUAL_TEST ) THEN
!          ! Override the DO-LOOPS parameters and use the
!          ! loops for repeating manual selection of parameters.
!          RJOBDATA(1:8) = 0 ! for recording job data, used in a
!                           ! Matlab code for testing purposes
!          WRITE(*,*) 'Manual selection of parameters'
!
!          WRITE(*,fmt='(a)',advance='no')   &
!          'Scaling parameter iSCALE or exit &
!          &(1=S, 2=C, 3=Y, 4=N, 5=STOP ) =  '
!          READ(*,*) iSCALEd
!          SELECT CASE ( iSCALEd )
!          CASE(1)
!              SCALE = 'S'
!          CASE(2)
!              SCALE = 'C'
!          CASE(3)
!              SCALE = 'Y'
!          CASE(4)
!              SCALE = 'N'
!          CASE(5)
!              WRITE(*,*) 'Stop manual testing/inspection. '
!              STOP
!          END SELECT
!
!          WRITE(*,fmt='(a)',advance='no') &
!          'Eigenvectors request JOBZ (1=V, 2=F, 3=N ) =  '
!          READ(*,*) iJOBZd
!          SELECT CASE ( iJOBZd )
!          CASE(1)
!              JOBZ   = 'V'
!              RESIDS = 'R'
!              RJOBDATA(6) = 1
!          CASE(2)
!              JOBZ   = 'F'
!              RESIDS = 'N'
!              RJOBDATA(7) = 1
!          CASE(3)
!              JOBZ   = 'N'
!              RESIDS = 'N'
!          END SELECT
!
!          WRITE(*,fmt='(a)',advance='no') &
!          'Eigenvectors refinement JOBF (1=R, 2=E, 3=N ) =  '
!          READ(*,*) iREF
!          SELECT CASE ( iREF )
!          CASE(1)
!              JOBREF = 'R'
!              RJOBDATA(5) = 1
!          CASE(2)
!              JOBREF = 'E'
!              RJOBDATA(4) = 1
!          CASE(3)
!              JOBREF = 'N'
!          END SELECT
!
!          WRITE(*,fmt='(a)',advance='no') &
!         'Select SVD truncation method NRNK ( -1, -2, 0 < K <= N ) = '
!          READ(*,*) NRNK
!          WRITE(*,fmt='(a)',advance='no') &
!          'Select the SVD method (1-4) WHTSVD = '
!          READ(*,*) WHTSVD
!          IF ( TEST_SP ) THEN
!              WRITE(*,fmt='(a)',advance='no') &
!             'SINGLE PREC. Select SVD truncation method NRNK ( -1, -2, 0 < K <= N ) = '
!              READ(*,*) NRNKsp
!              WRITE(*,fmt='(a)',advance='no') &
!             'SINGLE PREC. Select the SVD method (1-4) WHTSVDsp = '
!          READ(*,*) WHTSVDsp
!          END IF
!
!          WRITE(*,fmt='(a)',advance='no') &
!          'Truncation tolerance level (1=M*EPS, 2=sqrt(EPS), 3=EPS) = '
!          READ(*,*) iTOL
!          SELECT CASE (iTOL)
!          CASE(1)
!              TOL = M*EPS
!          CASE(2)
!              TOL = SQRT(EPS)
!          CASE(3)
!              TOL = EPS
!          END SELECT
!
!       END IF
!!...... end manual test section

!...... For each run, take the original data.

       X(1:M,1:N) = X0(1:M,1:N)
       Y(1:M,1:N) = Y0(1:M,1:N)

       ! DGEDMD: Workspace query and workspace allocation
!       WRITE(*,*) 'DGEDMD workspace query ...'
       CALL DGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD, M, &
            N, X, LDX, Y, LDY, NRNK, TOL, K, REIG, IEIG, Z, &
            LDZ, RES, AU, LDAU, W, LDW, S, LDS, WDUMMY, -1, &
            IDUMMY, -1, INFO )
!       IF ( (INFO .EQ. 2) .OR. ( INFO .EQ. 3 ) &
!                          .OR. ( INFO < 0 ) ) THEN
!           WRITE(*,*) 'Call to DGEDMD workspace query failed. &
!           &Check the calling sequence and the code.'
!           WRITE(*,*) 'The error code is ', INFO
!           WRITE(*,*) 'The input parameters were ', &
!           SCALE, JOBZ, RESIDS, JOBREF, WHTSVD, &
!           M, N, LDX, LDY, NRNK, TOL, LDZ, LDAU, LDW, LDS
!           STOP
!       ELSE
!           WRITE(*,*) '... done. Workspace length computed.'
!       END IF
       LIWORK = IDUMMY(1)
       ALLOCATE( IWORK(LIWORK) )
       LWORK = INT(WDUMMY(LWMINOPT))
       ALLOCATE( WORK(LWORK) )

       ! DGEDMD test: CALL DGEDMD
!       WRITE(*,*) 'Run DGEDMD ...'
       CALL DGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD, M, &
            N, X, LDX, Y, LDY, NRNK, TOL, K, REIG, IEIG, Z, &
            LDZ, RES, AU, LDAU, W, LDW, S, LDS, WORK, LWORK,&
            IWORK, LIWORK, INFO )
!       IF ( INFO /= 0 ) THEN
!           WRITE(*,*) 'Call to DGEDMD failed. &
!           &Check the calling sequence and the code.'
!           WRITE(*,*) 'The error code is ', INFO
!           WRITE(*,*) 'The input parameters were ',&
!           SCALE, JOBZ, RESIDS, JOBREF, WHTSVD, &
!           M, N, LDX, LDY, NRNK, TOL
!           STOP
!       END IF
!       WRITE(*,*) 'DGEDMD completed, INFO = ', INFO,'; ', &
!       'The number of computed Ritz pairs K = ', K

       SINGVX(1:N) = WORK(1:N)
!       WRITE(*,*) 'SIGMA(1) = ', SINGVX(1), '; ', &
!                  'SIGMA(K) = ', SINGVX(K), '; ', &
!                  'SIGMA(N) = ', SINGVX(N)
!       IF ( WHTSVD .EQ. 4 ) THEN
!           WRITE(*,*) 'DGEDMD: Scaling factors in the SVD ',&
!           WORK(N+1), WORK(N+2), WORK(N+2)/WORK(N+1)
!       END IF

!...... DGEDMD check point
       IF ( LSAME(JOBZ,'V')  ) THEN
          ! Check that Z = X*W, on return from DGEDMD
          ! This checks that the returned aigenvectors in Z are
          ! the product of the SVD'POD basis returned in X
          ! and the eigenvectors of the rayleigh quotient
          ! returned in W
          CALL DGEMM( 'N', 'N', M, K, K, ONE, X, LDX, W, LDW, &
                      ZERO, Z1, LDZ )
          TMP = ZERO
          DO i = 1, K
             CALL DAXPY( M, -ONE, Z(1,i), 1, Z1(1,i), 1)
             TMP = MAX(TMP, DNRM2( M, Z1(1,i), 1 ) )
          END DO
          TMP_ZXW = MAX(TMP_ZXW, TMP )
!          WRITE(*,*) 'Check #1: Z = X*V :: explicit comp. vs DGEDMD'
!          WRITE(*,*) 'Maximal column-wise error = ', TMP
!          WRITE(*,*) ' JOBZ was ', JOBZ
          IF ( TMP_ZXW <= 10*M*EPS ) THEN
!              WRITE(*,*) ' :) .... OK .........DGEDMD PASSED.'
          ELSE
              NFAIL_Z_XV = NFAIL_Z_XV + 1
!              WRITE(*,*) ':( .................DGEDMD FAILED!', &
!                  'Check the code for implementation errors.'
!              WRITE(*,*) 'The input parameters were ',&
!                 SCALE, JOBZ, RESIDS, JOBREF, WHTSVD, &
!                 M, N, LDX, LDY, NRNK, TOL
          END IF
!          IF ( KEYBIN )  THEN
!              WRITE(*,*) ' Hit ENTER to continue.'
!              READ(*,*)
!          END IF
       END IF
!...... DGEDMD check point
       IF ( LSAME(JOBREF,'R') ) THEN
           ! The matrix A*U is returned for computing refined Ritz vectors.
           ! Check that A*U is computed correctly using the formula
           ! A*U = Y * V * inv(SIGMA). This depends on the
           ! accuracy in the computed singular values and vectors of X.
           ! See the paper for an error analysis.
           ! Note that the left singular vectors of the input matrix X
           ! are returned in the array X.
           CALL DGEMM( 'N', 'N', M, K, M, ONE, A, LDA, X, LDX, &
                      ZERO, Z1, LDZ )
          TMP = ZERO
          DO i = 1, K
             CALL DAXPY( M, -ONE, AU(1,i), 1, Z1(1,i), 1)
            TMP = MAX( TMP, DNRM2( M, Z1(1,i),1 ) * &
                     SINGVX(K)/(ANORM*SINGVX(1)) )
          END DO
          TMP_AU = MAX( TMP_AU, TMP )
!          WRITE(*,*) 'Check #2: A*U :: explicit comp. vs DGEDMD'
!          WRITE(*,*) 'Maximal column-wise error = ', TMP
!          WRITE(*,*) ' JOBZ was ', JOBZ
          IF ( TMP <= TOL2 ) THEN
!              WRITE(*,*) ':) .... OK .........DGEDMD PASSED.', &
!                  ' Hit ENTER to continue.'
          ELSE
              NFAIL_AU = NFAIL_AU + 1
!              WRITE(*,*) ':( .................DGEDMD FAILED!', &
!                  'Check the code for implementation errors.'
!              WRITE(*,*) 'The input parameters were ',&
!                 SCALE, JOBZ, RESIDS, JOBREF, WHTSVD, &
!                 M, N, LDX, LDY, NRNK, TOL
          END IF
!          IF ( KEYBIN )  THEN
!            WRITE(*,*) 'Hit ENTER to continue.'
!            READ(*,*)
!          END IF
       ELSEIF ( LSAME(JOBREF,'E') ) THEN
       ! The unscaled vectors of the Exact DMD are computed.
       ! This option is included for the sake of completeness,
       ! for users who prefer the Exact DMD vectors. The
       ! returned vectors are in the real form, in the same way
       ! as the Ritz vectors. Here we just save the vectors
       ! and test them separately using a Matlab script.

!           IF ( MANUAL_TEST ) THEN
!              CALL matrix_save( 9, 'exvec.exdmd',  M, K, AU, LDAU )
!           ELSE
!               CALL DGEMM( 'N', 'N', M, K, M, ONE, A, LDA, AU, LDAU, ZERO, Y1, M )
!               i=1
!               DO WHILE ( i <= K )
!               IF ( IEIG(i) == ZERO ) THEN
!                ! have a real eigenvalue with real eigenvector
!                CALL DAXPY( M, -REIG(i), AU(1,i), 1, Y1(1,i), 1 )
!                RESEX(i) = DNRM2( M, Y1(1,i), 1) / DNRM2(M,AU(1,i),1)
!                i = i + 1
!               ELSE
!               ! Have a complex conjugate pair
!               ! REIG(i) +- sqrt(-1)*IMEIG(i).
!               ! Since all computation is done in real
!               ! arithmetic, the formula for the residual
!               ! is recast for real representation of the
!               ! complex conjugate eigenpair. See the
!               ! description of RES.
!               AB(1,1) =  REIG(i)
!               AB(2,1) = -IEIG(i)
!               AB(1,2) =  IEIG(i)
!               AB(2,2) =  REIG(i)
!               CALL DGEMM( 'N', 'N', M, 2, 2, -ONE, AU(1,i), &
!                           M, AB, 2, ONE, Y1(1,i), M )
!               RESEX(i)   = DLANGE( 'F', M, 2, Y1(1,i), M, &
!                            WORK )/ DLANGE( 'F', M, 2, AU(1,i), M, &
!                            WORK )
!               RESEX(i+1) = RESEX(i)
!               i = i + 2
!               END IF
!               END DO
!           END IF

       END IF
!...... DGEDMD check point
      IF ( LSAME(RESIDS, 'R') ) THEN
          ! Compare the residuals returned by DGEDMD with the
          ! explicitly computed residuals using the matrix A.
          ! Compute explicitly Y1 = A*Z
          CALL DGEMM( 'N', 'N', M, K, M, ONE, A, LDA, Z, LDZ, ZERO, Y1, M )
          ! ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
          ! of the invariant subspaces that correspond to complex conjugate
          ! pairs of eigencalues. (See the description of Z in DGEDMD,)
          i = 1
          DO WHILE ( i <= K )
            IF ( IEIG(i) == ZERO ) THEN
                ! have a real eigenvalue with real eigenvector
                CALL DAXPY( M, -REIG(i), Z(1,i), 1, Y1(1,i), 1 )
                RES1(i) = DNRM2( M, Y1(1,i), 1)
                i = i + 1
            ELSE
               ! Have a complex conjugate pair
               ! REIG(i) +- sqrt(-1)*IMEIG(i).
               ! Since all computation is done in real
               ! arithmetic, the formula for the residual
               ! is recast for real representation of the
               ! complex conjugate eigenpair. See the
               ! description of RES.
               AB(1,1) =  REIG(i)
               AB(2,1) = -IEIG(i)
               AB(1,2) =  IEIG(i)
               AB(2,2) =  REIG(i)
               CALL DGEMM( 'N', 'N', M, 2, 2, -ONE, Z(1,i), &
                           M, AB, 2, ONE, Y1(1,i), M )
               RES1(i)   = DLANGE( 'F', M, 2, Y1(1,i), M, &
                                  WORK )
               RES1(i+1) = RES1(i)
               i = i + 2
            END IF
          END DO
          TMP = ZERO
          DO i = 1, K
          TMP = MAX( TMP, ABS(RES(i) - RES1(i)) * &
                    SINGVX(K)/(ANORM*SINGVX(1)) )
          END DO
          TMP_REZ = MAX( TMP_REZ, TMP )
!          WRITE(*,*) &
!          'Check #3: Residuals :: explicit comp. vs DGEDMD'
!          WRITE(*,*) &
!          'Maximal rel. differences of the residuals = ', TMP
          IF ( TMP <= TOL2 ) THEN
!              WRITE(*,*) ':) .... OK ..........DGEDMD PASSED.'
          ELSE
              NFAIL_REZ = NFAIL_REZ + 1
!              WRITE(*,*) ':( ..................DGEDMD FAILED!', &
!                  'Check the code for implementation errors.'
!              WRITE(*,*) 'The input parameters were ',&
!                 SCALE, JOBZ, RESIDS, JOBREF, WHTSVD, &
!                 M, N, LDX, LDY, NRNK, TOL
          END IF
!         IF (KEYBIN ) THEN
!             WRITE(*,*) 'Hit ENTER to continue.'
!             READ(*,*)
!         END IF

         IF ( LSAME(JOBREF,'E') ) THEN
            TMP = ZERO
          DO i = 1, K
          TMP = MAX( TMP, ABS(RES1(i) - RESEX(i))/(RES1(i)+RESEX(i)) )
          END DO
          TMP_EX = MAX(TMP_EX,TMP)
         END IF

      END IF
!   ... store the results for inspection
      DO i = 1, K
          LAMBDA(i,1) = REIG(i)
          LAMBDA(i,2) = IEIG(i)
      END DO
!      IF ( MANUAL_TEST ) THEN
!      CALL matrix_save( 9, 'eigsdmd.dmd',  K, 2, LAMBDA, N )
!      CALL matrix_save( 9, 'evecsdmd.dmd', M, K, Z,      M )
!      IF ( LSAME(RESIDS, 'R') ) THEN
!         CALL matrix_save( 9, 'resdmd.dmd',   K, 1, RES, N )
!      END IF
!      CALL matrix_save( 9,  'svdx.dmd', N, 1,SINGVX, N  )
!      END IF
      DEALLOCATE(IWORK)
      DEALLOCATE(WORK)
!======================================================================
!     Now test the SGEDMD, if requested.
!======================================================================
!      IF ( TEST_SP ) THEN
!          RJOBDATA(1) = 1
!          ! The single precision code is tested
!          Xsp(1:M,1:N) = X0(1:M,1:N)
!          Ysp(1:M,1:N) = Y0(1:M,1:N)
!          ! SGEDMD: Workspace query and workspace allocation
!          WRITE(*,*) 'SGEDMD workspace query ...'
!          CALL SGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVDsp, M, &
!            N, Xsp, LDX, Ysp, LDY, NRNKsp, TOLsp, Ksp, REIGsp, IEIGsp, Zsp, &
!            LDZ, RESsp, AUsp, LDAU, Wsp, LDW, Ssp, LDS, WDUMMYsp, -1,&
!            IDUMMY, -1, INFO )
!          IF ( (INFO .EQ. 2) .OR. ( INFO .EQ. 3 ) &
!                          .OR. ( INFO < 0 ) ) THEN
!              WRITE(*,*) 'Call to SGEDMD workspace inquiry failed. &
!              &Check the calling sequence and the code.'
!              WRITE(*,*) 'The error code is ', INFO
!              WRITE(*,*) 'The input parameters were ', &
!              SCALE, JOBZ, RESIDS, JOBREF, WHTSVDsp, &
!              M, N, LDX, LDY, NRNKsp, TOLsp, LDZ, LDAU, LDW, LDS
!              STOP
!           ELSE
!               WRITE(*,*) '... done. Workspace length computed.'
!           END IF
!
!           LIWORK = IDUMMY(1)
!           LWORK  = INT( WDUMMYsp(LWMINOPT) )
!           ALLOCATE( IWORK(LIWORK) )
!           ALLOCATE( SWORK(LWORK)  )
!
!          CALL SGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVDsp, M, &
!            N, Xsp, LDX, Ysp, LDY, NRNKsp, TOLsp, Ksp, REIGsp, IEIGsp, Zsp, &
!            LDZ, RESsp, AUsp, LDAU, Wsp, LDW, Ssp, LDS, SWORK, LWORK,&
!            IWORK, LIWORK, INFO )
!          WRITE(*,*) 'SGEDMD completed, INFO, K = ', INFO, Ksp
!          LAMBDAsp(1:Ksp,1)=REIGsp(1:Ksp)
!          LAMBDAsp(1:Ksp,2)=IEIGsp(1:Ksp)
!          SINGVXsp = SWORK(1:N)
!          IF ( WHTSVD .EQ. 4 ) THEN
!             WRITE(*,*) 'SGEDMD: Scaling factors in the SVD ',&
!             SWORK(N+1), SWORK(N+2), SWORK(N+2)/SWORK(N+1)
!          END IF
!
!          IF ( MANUAL_TEST ) THEN
!          CALL matrix_save_sp( 9, 'eigsdmdsp.dmd',  Ksp, 2, LAMBDAsp, N )
!          IF ( LSAME(RESIDS, 'R') ) THEN
!              CALL matrix_save_sp( 9, 'evecsdmdsp.dmd',   M, Ksp, Zsp,    M )
!              CALL matrix_save_sp( 9, 'resdmdsp.dmd',  Ksp, 1, RESsp, N )
!          END IF
!              CALL matrix_save( 9,  'svdxsp.dmd', N, 1,SINGVXsp, N  )
!          END IF
!          DEALLOCATE(SWORK)
!          DEALLOCATE(IWORK)
!      END IF
!======================================================================
!     Now test the DGEDMDQ, if requested.
!======================================================================
      IF ( TEST_QRDMD ) THEN
          RJOBDATA(2) = 1
          WANTQ = 'Q'   ! Always return Q and R from the
          WANTR = 'R'   ! initial QR factorization
          F1 = F

          ! DGEDMDQ test: Workspace query and workspace allocation
!          WRITE(*,*) 'DGEDMDQ workspace query ...'
          CALL DGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, &
               JOBREF, WHTSVD, M, N+1, F1, LDF, X, LDX, Y, &
               LDY, NRNK, TOL, KQ, REIGQ, IEIGQ, Z, LDZ,   &
               RES, AU, LDAU, W, LDW, S, LDS, WDUMMY,      &
               -1, IDUMMY, -1, INFO )
          LIWORK = IDUMMY(1)
          ALLOCATE( IWORK(LIWORK) )
          LWORK = INT(WDUMMY(LWMINOPT))
          ALLOCATE(WORK(LWORK))
          ! DGEDMDQ test: CALL DGEDMDQ
          CALL DGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, &
               JOBREF, WHTSVD, M, N+1, F1, LDF, X, LDX, Y, &
               LDY, NRNK, TOL, KQ, REIGQ, IEIGQ, Z, LDZ,   &
               RES, AU, LDAU, W, LDW, S, LDS,              &
               WORK, LWORK, IWORK, LIWORK, INFO )
!          IF ( INFO /= 0 ) THEN
!             WRITE(*,*) 'Call to DGEDMDQ failed. &
!             &Check the calling sequence and the code.'
!             WRITE(*,*) 'The error code is ', INFO
!             WRITE(*,*) 'The input parameters were ',&
!             SCALE, JOBZ, RESIDS, WANTQ, WANTR, WHTSVD, &
!             M, N, LDX, LDY, NRNK, TOL
!             STOP
!          END IF
!          WRITE(*,*)'DGEDMDQ completed, INFO, K = ', INFO, KQ
          SINGVQX(1:KQ) = WORK(MIN(M,N+1)+1: MIN(M,N+1)+KQ)
!..... DGEDMDQ check point
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
!..... DGEDMDQ check point
      IF ( LSAME(WANTQ,'Q') .AND. LSAME(WANTR,'R') ) THEN
         ! Check that the QR factors are computed and returned
         ! as requested. The residual ||F-Q*R||_F / ||F||_F
         ! is compared to M*N*EPS.
         F2 = F
         CALL DGEMM( 'N', 'N', M, N+1, MIN(M,N+1), -ONE, F1, &
                     LDF, Y, LDY, ONE, F2, LDF )
         TMP_FQR = DLANGE( 'F', M, N+1, F2, LDF, WORK ) / &
               DLANGE( 'F', M, N+1, F,  LDF, WORK )
!         WRITE(*,*) ' F = Q*R test :: error = ', TMP_FQR
         IF ( TMP_FQR > TOL2 ) THEN
             NFAIL_F_QR = NFAIL_F_QR + 1
!         ELSE
!             WRITE(*,*) '........ PASSED.'
         END IF
      END IF
!..... DGEDMDQ checkpoint
      IF ( LSAME(RESIDS, 'R') ) THEN
          ! Compare the residuals returned by DGEDMDQ with the
          ! explicitly computed residuals using the matrix A.
          ! Compute explicitly Y1 = A*Z
          CALL DGEMM( 'N', 'N', M, KQ, M, ONE, A, M, Z, M, ZERO, Y1, M )
          ! ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
          ! of the invariant subspaces that correspond to complex conjugate
          ! pairs of eigencalues. (See the description of Z in DGEDMDQ)
          i = 1
          DO WHILE ( i <= KQ )
            IF ( IEIGQ(i) == ZERO ) THEN
                ! have a real eigenvalue with real eigenvector
                CALL DAXPY( M, -REIGQ(i), Z(1,i), 1, Y1(1,i), 1 )
                ! Y(1:M,i) = Y(1:M,i) - REIG(i)*Z(1:M,i)
                RES1(i) = DNRM2( M, Y1(1,i), 1)
                i = i + 1
            ELSE
               ! Have a complex conjugate pair
               ! REIG(i) +- sqrt(-1)*IMEIG(i).
               ! Since all computation is done in real
               ! arithmetic, the formula for the residual
               ! is recast for real representation of the
               ! complex conjugate eigenpair. See the
               ! description of RES.
               AB(1,1) =  REIGQ(i)
               AB(2,1) = -IEIGQ(i)
               AB(1,2) =  IEIGQ(i)
               AB(2,2) =  REIGQ(i)
               CALL DGEMM( 'N', 'N', M, 2, 2, -ONE, Z(1,i), &
                           M, AB, 2, ONE, Y1(1,i), M )             ! BLAS CALL
               ! Y(1:M,i:i+1) = Y(1:M,i:i+1) - Z(1:M,i:i+1) * AB   ! INTRINSIC
               RES1(i)   = DLANGE( 'F', M, 2, Y1(1,i), M, &
                                  WORK )                           ! LAPACK CALL
               RES1(i+1) = RES1(i)
               i = i + 2
            END IF
          END DO
          TMP = ZERO
          DO i = 1, KQ
          TMP = MAX( TMP, ABS(RES(i) - RES1(i)) * &
              SINGVQX(K)/(ANORM*SINGVQX(1)) )
          END DO
          TMP_REZQ = MAX( TMP_REZQ, TMP )
!          WRITE(*,*) &
!          'Check #2: Residuals :: explicit comp. vs DGEDMDQ'
!          WRITE(*,*) &
!          'Maximal rel. differences of the residuals = ', TMP
          IF ( TMP <= TOL2 ) THEN
!              WRITE(*,*) '.... OK ........ DGEDMDQ PASSED.'
          ELSE
              NFAIL_REZQ = NFAIL_REZQ + 1
!              WRITE(*,*) '................ DGEDMDQ FAILED!', &
!                  'Check the code for implementation errors.'
          END IF
!         IF ( KEYBIN ) THEN
!             WRITE(*,*) ' Hit ENTER to continue.'
!             READ(*,*)
!         END IF

      END IF

      DO i = 1, KQ
          LAMBDAQ(i,1) = REIGQ(i)
          LAMBDAQ(i,2) = IEIGQ(i)
      END DO

!      IF ( MANUAL_TEST ) THEN
!      CALL matrix_save( 9,  'eigsqrdmd.qrdmd', KQ, 2,LAMBDAQ, N  )
!      IF ( LSAME(RESIDS, 'R') ) THEN
!      CALL matrix_save( 9, 'evecsqrdmd.qrdmd', M, KQ,    Z, M )
!      CALL matrix_save( 9,'resqrdmd.qrdmd', KQ, 1,   RES,    KQ )
!      WRITE(*,*) 'Output saved to files,'
!      END IF
!      IF ( LSAME(WANTQ,'Q') ) THEN
!          CALL matrix_save( 9, 'QQ.qrdmd', M, N+1, F1, LDF )
!      END IF
!      IF ( LSAME(JOBREF,'E') ) THEN
!          CALL matrix_save( 9, 'EXZ.qrdmd', (N+1), KQ, AU, LDAU )
!      ELSE IF ( LSAME(JOBREF,'R') ) THEN
!          CALL matrix_save( 9, 'AU.qrdmd', (N+1), KQ, AU, LDAU )
!          CALL matrix_save( 9, 'XU.qrdmd', (N+1), KQ,  X, LDX )
!      END IF
!      END IF

      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)
      END IF            ! TEST_QRDMD
!======================================================================

!======================================================================
!     Now test the SGEDMDQ, if requested.
!======================================================================
!      IF ( TESTQ_SP ) THEN
!          RJOBDATA(3) = 1
!          ! The single precision code is tested
!          F1sp(1:M,1:N+1) = F(1:M,1:N+1)
!          CALL SGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, WHTSVD, M, N+1, &
!                    F1sp, LDF, Xsp, LDX, Ysp, LDY, NRNKsp, TOLsp, KQsp, REIGQsp,   &
!                    IEIGQsp, Zsp, LDZ, RESsp, AUsp, LDAU, Wsp, LDW, Ssp, LDS, WDUMMYsp, &
!                    -1, IDUMMY, -1, INFO )
!          LIWORK = IDUMMY(1)
!          LWORK  = INT( WDUMMYsp(LWMINOPT) )
!          ALLOCATE( IWORK(LIWORK) )
!          ALLOCATE( SWORK(LWORK)  )
!
!          CALL SGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, WHTSVD, M, N+1, &
!                    F1sp, LDF, Xsp, LDX, Ysp, LDY, NRNKsp, TOLsp, KQsp, REIGQsp,   &
!                    IEIGQsp, Zsp, LDZ, RESsp, AUsp, LDAU, Wsp, LDW, Ssp, LDS, SWORK, &
!                    LWORK, IWORK, LIWORK, INFO )
!          WRITE(*,*) 'SGEDMDQ completed, INFO, K = ', INFO, KQsp
!
!          DO i = 1, KQsp
!             LAMBDAQsp(i,1) = REIGQsp(i)
!             LAMBDAQsp(i,2) = IEIGQsp(i)
!          END DO
!
!          IF ( MANUAL_TEST ) THEN
!          CALL matrix_save_sp( 9,  'eigsqrdmdsp.qrdmd', KQsp, 2,LAMBDAQsp, N  )
!          IF ( LSAME(RESIDS, 'R') ) THEN
!               CALL matrix_save_sp( 9, 'evecsqrdmdsp.qrdmd', M, KQsp,    Zsp, M )
!               CALL matrix_save_sp( 9,'resqrdmdsp.qrdmd', KQsp, 1,   RESsp,    KQsp )
!               WRITE(*,*) 'Output saved to files,'
!          END IF
!          IF ( LSAME(WANTQ,'Q') ) THEN
!             CALL matrix_save_sp( 9, 'QQsp.qrdmd', M, N+1, F1sp, LDF )
!          END IF
!          IF ( LSAME(JOBREF,'E') ) THEN
!              CALL matrix_save_sp( 9, 'EXZsp.qrdmd', (N+1), KQsp, AUsp, LDAU )
!          ELSE IF ( LSAME(JOBREF,'R') ) THEN
!              CALL matrix_save_sp( 9, 'AUsp.qrdmd', (N+1), KQsp, AUsp, LDAU )
!              CALL matrix_save_sp( 9, 'XUsp.qrdmd', (N+1), KQsp,  Xsp, LDX )
!          END IF
!          END IF
!          DEALLOCATE(IWORK)
!          DEALLOCATE(SWORK)
!      END IF ! TEST_SP
!================================================================================

!    WRITE(*,*) 'The scaling option was SCALE = ', SCALE
!    WRITE(*,*) 'The eigenvectors/modes option JOBZ was = ', JOBZ
!    WRITE(*,*) 'The residuals option RESIDS was = ', RESIDS
!    WRITE(*,*) 'The refinement data option was JOBREF = ', JOBREF
!    WRITE(*,*) 'The numerical rank option NRNK = ', NRNK
!    WRITE(*,*) 'The tolerance level was - ', iTOL, TOL
!    WRITE(*,*) 'The SVD option was WHTSVD = ', WHTSVD
!    WRITE(*,*) 'The number of modes computed by DGEDMD = ', K
!    IF ( TEST_QRDMD ) &
!    WRITE(*,*) 'The number of modes computed by DGEDMDQ = ', KQ
!    WRITE(*,*)
!    IF ( MANUAL_TEST ) THEN
!         OPEN( UNIT   = 9,         FILE   = 'RJOBDATA.data', &
!               STATUS = 'replace', ACTION = 'write' )
!          DO i = 1, 8
!             WRITE(9, FMT="(I4)", advance="NO") RJOBDATA(i)
!          END DO
!          CLOSE(9)
!          WRITE(*,*) 'Output saved to files.'
!      END IF
!    IF ( VISUAL_CHECK ) THEN
!    WRITE(*,*) 'Pause for visual inspection (Matlab script)'
!    WRITE(*,*) 'Hit Enter to continue.'
!    READ(*,*)
!    END IF
      END DO ! LWMINOPT
      END DO ! WHTSVD LOOP
      END DO ! NRNK LOOP
      END DO ! SCALE LOOP
      END DO ! JOBF LOOP
      END DO ! JOBZ LOOP


!      WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
!      WRITE(*,*) ' Test summary for DGEDMD :'
!      WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
!      WRITE(*,*)
      IF ( NFAIL_Z_XV == 0 ) THEN
!          WRITE(*,*) 'Z - U*V test PASSED.'
      ELSE
!          WRITE(*,*) 'Z - U*V test FAILED ', NFAIL_Z_XV, ' time(s)'
!          WRITE(*,*) 'Max error ||Z-U*V||_F was ', TMP_ZXW
          NFAIL_TOTAL = NFAIL_TOTAL + 1
      END IF
!      WRITE(*,*)
      IF ( NFAIL_AU == 0 ) THEN
!          WRITE(*,*) 'A*U test PASSED. '
      ELSE
!          WRITE(*,*) 'A*U test FAILED ', NFAIL_AU, ' time(s)'
!          WRITE(*,*) 'Max A*U test adjusted error measure was ', TMP_AU
!          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAIL_TOTAL = NFAIL_TOTAL + 1
      END IF
!      WRITE(*,*)
      IF ( NFAIL_REZ == 0 ) THEN
!          WRITE(*,*) 'Rezidual computation test PASSED.'
      ELSE
!          WRITE(*,*) 'Rezidual computation test FAILED ', NFAIL_REZ, 'time(s)'
!          WRITE(*,*) 'Max residual computing test adjusted error measure was ', TMP_REZ
!          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAIL_TOTAL = NFAIL_TOTAL + 1
      END IF
!      WRITE(*,*) 'When Exact DMD vectors were coputed and residuals of the DMD &
!          &vectors were available, then the largest relative difference of the &
!          & residuals over all vectors in all cases was ', TMP_EX, '. Use the &
!          &MANUAL_TEST option to visualize the residuals.   '

      IF ( TEST_QRDMD ) THEN
!      WRITE(*,*)
!      WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
!      WRITE(*,*) ' Test summary for DGEDMDQ :'
!      WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
!      WRITE(*,*)
      IF ( NFAIL_F_QR == 0 ) THEN
!          WRITE(*,*) 'F - Q*R test PASSED.'
      ELSE
!          WRITE(*,*) 'F - Q*R test FAILED ', NFAIL_F_QR, ' time(s)'
!          WRITE(*,*) 'The largest relative residual was ', TMP_FQR
!          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAIL_TOTAL = NFAIL_TOTAL + 1
      END IF
!      WRITE(*,*)
      IF ( KDIFF == 0 ) THEN
!          WRITE(*,*) 'DGESVD and DGESVDQ agreed in the number &
!              &of modes in all test cases. PASSED.'
      ELSE
!          WRITE(*,*) 'DGESVD and DGESVDQ computed different &
!              &number of modes ', KDIFF, ' times.'
      END IF
!      WRITE(*,*)
      IF ( NFAIL_SVDIFF == 0 ) THEN
!          WRITE(*,*) 'DGESVD and DGESVDQ computed singular &
!              &values test:'
!          WRITE(*,*) 'The maximal discrepancy in the singular &
!              &values was ', SVDIFF, ' Test PASSED.'
      ELSE
!          WRITE(*,*) 'DGESVD and DGESVDQ discrepancies in &
!              &the singular values unacceptable ', &
!              NFAIL_SVDIFF, ' times'
!          WRITE(*,*) 'The maximal discrepancy in the singular values was ', SVDIFF
!          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAIL_TOTAL = NFAIL_TOTAL + 1
      END IF
!      WRITE(*,*)
      IF ( NFAIL_REZQ == 0 ) THEN
!          WRITE(*,*) 'Rezidual computation test PASSED.'
      ELSE
!          WRITE(*,*) 'Rezidual computation test FAILED ', NFAIL_REZQ, 'time(s)'
!          WRITE(*,*) 'Max residual computing test adjusted error measure was ', TMP_REZQ
!          WRITE(*,*) 'It should be up to O(M*N) times EPS, EPS = ', EPS
          NFAIL_TOTAL = NFAIL_TOTAL + 1
      END IF
!      WRITE(*,*)
      END IF
      IF ( NFAIL_TOTAL == 0 ) THEN
          WRITE(*,*) '>>>>>>> ALL TESTS PASSED.'
      ELSE
          WRITE(*,*) NFAIL_TOTAL, &
          '>>>>>>> TESTS FAILED. CHECK THE IMPLEMENTATION.'
      END IF


      ELSE

!============================================================================
!============================================================================
! Test complex data subroutines
!============================================================================
!============================================================================

!
!
!      COND = 1.0D8
!      ZMAX = (1.0D2,1.0D2)
!      RSIGN = 'F'
!      GRADE = 'N'
!      MODEL = 6
!      CONDL = 1.0D2
!      MODER = 6
!      CONDR = 1.0D2
!      PIVTNG = 'N'
!      ALLOCATE( IWORK(2*M) )
!      ALLOCATE( ZDA(M) )
!      ALLOCATE( ZA(LDA,M) )
!      ALLOCATE( ZAC(LDA,M) )
!      ALLOCATE(ZDL(M))
!      ALLOCATE(ZDR(M))
!      CALL ZLATMR( M, M, 'S', ISEED, 'N', ZDA, MODE, COND, &
!                   ZMAX, RSIGN, GRADE, ZDL, MODEL,  CONDL, &
!                   ZDR, MODER, CONDR, PIVTNG, IWORK, M, M, &
!                   ZERO, -ONE, 'N', ZA, LDA, IWORK(M+1), INFO )
!      DEALLOCATE(IWORK)
!
!
!      ! Store reference eigenvalues and eigenvectors
!      WRITE(*,'(A)', advance='no') &
!      'Compute and store the eigenvalues &
!                 &and eigenvectors of A (ZGEEV) . ..'
!
!      LZWORK = MAX(1,2*M)
!      ALLOCATE( ZEIGSA(M) )
!      ALLOCATE(ZVA(LDA,M))
!      ALLOCATE( ZWORK(LZWORK) )
!      ALLOCATE( WORK(2*M) )
!      ZAC = ZA
!      CALL ZGEEV( 'N','V', M, ZAC, LDA, ZEIGSA, ZVA, M, &
!                  ZVA, M, ZWORK, LZWORK, WORK, INFO ) ! LAPACK CALL
!      WRITE(*,*) '... completed.'
!      DEALLOCATE(WORK)
!      DEALLOCATE(ZWORK)
!
!      TMP = ABS(ZEIGSA(IZAMAX(M, ZEIGSA, 1))) ! The spectral radius of ZA
!
!      CALL ZLASCL( 'G',0, 0, TMP, ONE, M, M, &
!                   ZA, LDA, INFO )
!      CALL ZLASCL( 'G',0, 0, TMP, ONE, M, 1, &
!                   ZEIGSA, M, INFO )
!
!      IF ( MANUAL_TEST ) THEN
!      CALL matrix_save_Zre( 9, 'ZAre.data', M, M, ZA, LDA )
!      CALL matrix_save_Zim( 9, 'ZAim.data', M, M, ZA, LDA )
!      CALL matrix_save_Zre( 9, 'ZEIGSAre.data', M, 1, ZEIGSA, M )
!      CALL matrix_save_Zim( 9, 'ZEIGSAim.data', M, 1, ZEIGSA, M )
!      END IF
!
!      ALLOCATE( ZF(LDF,N+1) )
!      ALLOCATE( ZF0(LDF,N+1) )
!      CALL ZLARNV(2, ISEED, M, ZF(1,1) )
!      DO i = 1, N
!         CALL ZGEMV( 'N', M, M, ZONE, ZA, M, ZF(1,i), 1, ZZERO, &
!              ZF(1,i+1), 1 )
!      END DO
!      ZF0(1:M,1:N+1) = ZF(1:M,1:N+1)
!!........................................................................
!
!      ALLOCATE( ZX(LDX,N) )
!      ALLOCATE( ZY(LDY,N+1) )
!      ALLOCATE( ZAU(LDAU,N) )
!      ALLOCATE( ZW(LDW,N) )
!      ALLOCATE( ZS(LDS,N) )
!      ALLOCATE( ZZ(LDZ,N) )
!      ALLOCATE( RES(N) )
!      ALLOCATE( ZEIGS(N) )
!
!      ZX(1:M,1:N) = ZF0(1:M,1:N)
!      ZY(1:M,1:N) = ZF0(1:M,2:N+1)
!
!      IF ( TEST_SP ) THEN
!          ALLOCATE( ZXsp(LDX,N) )
!          ALLOCATE( ZYsp(LDY,N+1) )
!          ALLOCATE( ZFsp(LDF,N+1) )
!          ALLOCATE( ZEIGSsp(N))
!          ALLOCATE( RESsp(N) )
!          ALLOCATE( ZWsp(LDW,N) )
!          ALLOCATE( ZAUsp(LDAU,N) )
!          ALLOCATE( ZZsp(LDZ,N) )
!          ALLOCATE( ZSsp(LDS,N) )
!          ZXsp(1:M,1:N) = ZX(1:M,1:N)
!          ZYsp(1:M,1:N) = ZY(1:M,1:N)
!          ZFsp(1:M,1:N+1) = ZF(1:M,1:N+1)
!      END IF
!
!
!      IF ( MANUAL_TEST ) THEN
!      CALL matrix_save_Zre( 9, 'ZXre.data', M, N, ZX, LDX  )
!      CALL matrix_save_Zim( 9, 'ZXim.data', M, N, ZX, LDX )
!      CALL matrix_save_Zre( 9, 'ZYre.data', M, N, ZY, LDY  )
!      CALL matrix_save_Zim( 9, 'ZYim.data', M, N, ZY, LDY )
!      END IF
!
!      TOL   = M*EPS
!      TOLsp = M*EPSsp
!      ! here nested doo loops
!
!            DO iJOBZ = 1, 4
!
!          SELECT CASE ( iJOBZ )
!          CASE(1)
!              JOBZ   = 'V'
!              RESIDS = 'R'
!          CASE(2)
!              JOBZ   = 'V'
!              RESIDS = 'N'
!          CASE(3)
!              JOBZ   = 'F'
!              RESIDS = 'N'
!          CASE(4)
!              JOBZ   = 'N'
!              RESIDS = 'N'
!          END SELECT
!
!      DO iJOBREF = 1, 3
!
!          SELECT CASE ( iJOBREF )
!          CASE(1)
!              JOBREF = 'R'
!          CASE(2)
!              JOBREF = 'E'
!          CASE(3)
!              JOBREF = 'N'
!          END SELECT
!
!      DO iSCALE = 1, 4
!
!          SELECT CASE ( iSCALE )
!          CASE(1)
!              SCALE = 'S'
!          CASE(2)
!              SCALE = 'C'
!          CASE(3)
!              SCALE = 'Y'
!          CASE(4)
!              SCALE = 'N'
!          END SELECT
!
!      DO iNRNK = -1, -2, -1
!         NRNK   = iNRNK
!         NRNKsp = iNRNK
!
!      DO iWHTSVD = 1, 4
!         ! Check all four options to compute the POD basis
!         ! via the SVD.
!         WHTSVD   = iWHTSVD
!         WHTSVDsp = iWHTSVD
!      DO LWMINOPT = 1, 2
!         ! Workspace query for the minimal (1) and for the optimal
!         ! (2) workspace lengths determined by workspace query.
!
!       IF ( MANUAL_TEST ) THEN
!          ! Override the DO-LOOPS parameters and use the
!          ! loops for repeating manual selection of parameters.
!          ZJOBDATA(1:8) = 0 ! for recording job data, used in a
!                           ! Matlab code for testing purposes
!
!          WRITE(*,*) 'Manual selection of parameters'
!          WRITE(*,fmt='(a)',advance='no') &
!          'Scaling parameter iSCALE (1=S, 2=C, 3=Y, 4=N ) =  '
!          READ(*,*) iSCALEd
!          SELECT CASE ( iSCALEd )
!          CASE(1)
!              SCALE = 'S'
!          CASE(2)
!              SCALE = 'C'
!          CASE(3)
!              SCALE = 'Y'
!          CASE(4)
!              SCALE = 'N'
!          END SELECT
!
!          WRITE(*,fmt='(a)',advance='no') &
!          'Eigenvectors request JOBZ (1=V, 2=F, 3=N ) =  '
!          READ(*,*) iJOBZd
!          SELECT CASE ( iJOBZd )
!          CASE(1)
!              JOBZ   = 'V'
!              RESIDS = 'R'
!              ZJOBDATA(1) = 1
!              ZJOBDATA(2) = 1
!          CASE(2)
!              JOBZ   = 'F'
!              RESIDS = 'N'
!              ZJOBDATA(4) = 1
!          CASE(3)
!              JOBZ   = 'N'
!              RESIDS = 'N'
!          END SELECT
!
!          WRITE(*,fmt='(a)',advance='no') &
!          'Eigenvectors refinement JOBF (1=RefineData, 2=ExactDMD, 3=N ) =  '
!          READ(*,*) iREF
!          SELECT CASE ( iREF )
!          CASE(1)
!              JOBREF = 'R'
!              ZJOBDATA(5) = 1
!          CASE(2)
!              JOBREF = 'E'
!              ZJOBDATA(3) = 1
!          CASE(3)
!              JOBREF = 'N'
!          END SELECT
!
!          WRITE(*,fmt='(a)',advance='no') &
!         'Select SVD truncation method NRNK ( -1, -2, 0 < K <= N ) = '
!          READ(*,*) NRNK
!          WRITE(*,fmt='(a)',advance='no') &
!          'Select the SVD method (1-4) WHTSVD = '
!          READ(*,*) WHTSVD
!          IF ( TEST_SP ) THEN
!              WRITE(*,fmt='(a)',advance='no') &
!             'SINGLE PREC. Select SVD truncation method NRNK ( -1, -2, 0 < K <= N ) = '
!              READ(*,*) NRNKsp
!              WRITE(*,fmt='(a)',advance='no') &
!             'SINGLE PREC. Select the SVD method (1-4) WHTSVDsp = '
!          READ(*,*) WHTSVDsp
!          END IF
!
!          WRITE(*,fmt='(a)',advance='no') &
!          'Truncation tolerance level (1=M*EPS, 2=sqrt(EPS), 3=EPS) = '
!          READ(*,*) iTOL
!          SELECT CASE (iTOL)
!          CASE(1)
!              TOL = M*EPS
!          CASE(2)
!              TOL = SQRT(EPS)
!          CASE(3)
!              TOL = EPS
!          END SELECT
!
!       END IF
!!...... end manual setup of test options
!
!
!
!      ZX(1:M,1:N) = ZF0(1:M,1:N)
!      ZY(1:M,1:N) = ZF0(1:M,2:N+1)
!      ZF(1:M,1:N+1) = ZF0(1:M,1:N+1)
!      WRITE(*,*) 'Test data generated.'
!      WRITE(*,*) 'Call ZGEDMD for workspace length ...'
!
!      CALL ZGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD,  &
!                   M,  N, ZX, LDX, ZY, LDY, NRNK, TOL,  &
!                   K, ZEIGS, ZZ, LDZ,  RES,  &
!                   ZAU, LDAU, ZW,  LDW,   ZS, LDS,        &
!                   ZDUMMY, -1, WDUMMY, LWORK, IDUMMY, LIWORK, INFO )
!      WRITE(*,*) 'Done! INFO = ', INFO
!
!     LZWORK = INT(ZDUMMY(LWMINOPT))
!     ALLOCATE(ZWORK(LZWORK))
!     LIWORK = IDUMMY(1)
!     ALLOCATE(IWORK(LIWORK))
!     LWORK = INT(WDUMMY(1))
!     ALLOCATE(WORK(LWORK))
!     WRITE(*,*) 'Run ZGEDMD ....'
!     CALL ZGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD,   &
!                   M,  N, ZX, LDX, ZY, LDY, NRNK, TOL,   &
!                   K, ZEIGS, ZZ, LDZ,  RES, ZAU, LDAU,   &
!                   ZW,  LDW,   ZS, LDS,ZWORK, LZWORK,    &
!                   WORK, LWORK, IWORK, LIWORK, INFO )
!      WRITE(*,*) 'Done! INFO = ', INFO
!
!
!      IF ( MANUAL_TEST ) THEN
!
!      CALL matrix_save_Zre( 9, 'ZEIGSre.data', K, 1, ZEIGS, K )
!      CALL matrix_save_Zim( 9, 'ZEIGSim.data', K, 1, ZEIGS, K )
!      IF ( LSAME(JOBZ,'V') ) THEN
!      CALL matrix_save_Zre( 9, 'ZEVECre.data', M, K, ZZ, LDZ )
!      CALL matrix_save_Zim( 9, 'ZEVECim.data', M, K, ZZ, LDZ )
!      END IF
!      IF ( LSAME(RESIDS,'R') ) &
!          CALL matrix_save( 9, 'ZRES.data', K, 1, RES, K )
!      IF ( LSAME(JOBREF,'E') ) THEN
!           ! The unscaled vectors of the Exact DMD are computed.
!           ! This option is included for the sake of completeness,
!           ! for users who prefer the Exact DMD vectors. The
!           ! returned vectors are saved in the real form, and can be
!           ! tested separately using a Matlab script.
!           CALL matrix_save_Zre( 9, 'ZEXVECre.exdmd',  M, K, ZAU, LDAU )
!           CALL matrix_save_Zim( 9, 'ZEXVECim.exdmd',  M, K, ZAU, LDAU )
!      END IF
!      IF ( LSAME(JOBREF,'R') ) THEN
!           ! The returned vectors are saved in the real form, and
!           ! can be tested separately using a Matlab script.
!           CALL matrix_save_Zre( 9, 'ZREFre.exdmd',  M, K, ZAU, LDAU )
!           CALL matrix_save_Zim( 9, 'ZREFim.exdmd',  M, K, ZAU, LDAU )
!      END IF
!           ! The eigenvectors are returned in factored form.
!      IF ( LSAME(JOBZ,'F') .OR. LSAME(JOBZ,'V')) THEN
!           CALL matrix_save_Zre( 9, 'ZWre.dmd',  K, K, ZW, LDW )
!           CALL matrix_save_Zim( 9, 'ZWim.dmd',  K, K, ZW, LDW )
!           CALL matrix_save_Zre( 9, 'ZUkre.dmd',  M, K, ZX, LDX )
!           CALL matrix_save_Zim( 9, 'ZUkim.dmd',  M, K, ZX, LDX )
!      END IF
!
!      END IF ! manual test ... .saving the results for case study/analysis
!
!      DEALLOCATE( ZWORK )
!      DEALLOCATE( WORK  )
!      DEALLOCATE( IWORK )
!
!      IF ( TEST_QRDMD ) THEN
!         ZJOBDATA(6) = 1
!         WANTQ = 'Q'
!         WANTR = 'R'
!        WRITE(*,*) 'Call ZGEDMDQ for workspace length ...'
!      CALL ZGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, &
!                    WHTSVD, M, N+1, ZF, LDF,  ZX, LDX,  ZY, LDY,  &
!                    NRNK,  TOL, K, ZEIGS, ZZ, LDZ, RES,  ZAU,  &
!                    LDAU, ZW, LDW, ZS, LDS, ZDUMMY, -1,   &
!                    WDUMMY,  -1, IDUMMY, -1, INFO )
!      WRITE(*,*) 'Done! INFO = ', INFO
!
!      LZWORK = INT(ZDUMMY(LWMINOPT))
!      ALLOCATE( ZWORK(LZWORK) )
!      LIWORK = IDUMMY(1)
!      ALLOCATE(IWORK(LIWORK))
!      LWORK = INT(WDUMMY(1))
!      ALLOCATE(WORK(LWORK))
!
!      WRITE(*,*) 'Run ZGEDMDQ ....'
!      CALL ZGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, &
!                    WHTSVD, M, N+1, ZF, LDF,  ZX, LDX,  ZY, LDY,  &
!                    NRNK,  TOL, KQ, ZEIGS, ZZ, LDZ, RES,  ZAU,  &
!                    LDAU, ZW, LDW, ZS, LDS, ZWORK, LZWORK,   &
!                    WORK,  LWORK, IWORK, LIWORK, INFO )
!      WRITE(*,*) 'Done! INFO = ', INFO
!
!    IF ( MANUAL_TEST ) THEN
!       CALL matrix_save_Zre( 9, 'QZEIGSre.data', K, 1, ZEIGS, K )
!       CALL matrix_save_Zim( 9, 'QZEIGSim.data', K, 1, ZEIGS, K )
!
!      IF ( LSAME(JOBZ,'V') ) THEN
!      CALL matrix_save_Zre( 9, 'QZEVECre.data', M, K, (ZZ), LDZ )
!      CALL matrix_save_Zim( 9, 'QZEVECim.data', M, K, (ZZ), LDZ )
!      END IF
!
!       IF ( LSAME(RESIDS, 'R') ) THEN
!       CALL matrix_save( 9,'QZRES.zgedmdq', KQ, 1,   RES,    KQ )
!       END IF
!
!       IF ( LSAME(WANTQ,'Q') ) THEN
!         CALL matrix_save_Zre( 9, 'QZQQre.data', M, MIN(M,N+1), ZF, LDF )
!         CALL matrix_save_Zim( 9, 'QZQQim.data', M, MIN(M,N+1), ZF, LDF )
!       END IF
!
!       IF ( LSAME(JOBREF,'E') ) THEN
!           ! The unscaled vectors of the Exact DMD are computed.
!           ! This option is included for the sake of completeness,
!           ! for users who prefer the Exact DMD vectors. The
!           ! returned vectors are saved in the real form, and can be
!           ! tested separately using a Matlab script.
!           CALL matrix_save_Zre( 9, 'QZEXVECre.exdmd',  MIN(M,N+1), KQ, ZAU, LDAU )
!           CALL matrix_save_Zim( 9, 'QZEXVECim.exdmd',  MIN(M,N+1), KQ, ZAU, LDAU )
!       ELSE IF ( LSAME(JOBREF,'R') ) THEN
!           CALL matrix_save_Zre( 9, 'QZAUre.data',  N+1, KQ, ZAU, LDAU )
!           CALL matrix_save_Zim( 9, 'QZAUim.data',  N+1, KQ, ZAU, LDAU )
!           CALL matrix_save_Zre( 9, 'QZXUre.qrdmd', (N+1), KQ,  ZX, LDX )
!           CALL matrix_save_Zim( 9, 'QZXUim.qrdmd', (N+1), KQ,  ZX, LDX )
!       END IF
!
!       IF ( LSAME(JOBZ,'F') .OR. LSAME(JOBZ,'V')) THEN
!           CALL matrix_save_Zre( 9, 'QZWre.dmd',  K, K, ZW, LDW )
!           CALL matrix_save_Zim( 9, 'QZWim.dmd',  K, K, ZW, LDW )
!           CALL matrix_save_Zre( 9, 'QZUkre.dmd',  M, K, ZX, LDX )
!           CALL matrix_save_Zim( 9, 'QZUkim.dmd',  M, K, ZX, LDX )
!      END IF
!
!
!    END IF  ! MANUAL_TEST
!    END IF ! TEST_QRDMD
!
!    !DEALLOCATE(IWORK)
!    !DEALLOCATE(ZWORK)
!    !DEALLOCATE(WORK)
!    !DEALLOCATE(ZF)
!
!      IF ( TEST_SP ) THEN
!       ZJOBDATA(7) = 1
!         TOLsp = M*EPSsp
!         ZXsp(1:M,1:N) = ZF0(1:M,1:N)
!         ZYsp(1:M,1:N) = ZF0(1:M,2:N+1)
!          WRITE(*,*) 'Call CGEDMD for workspace length ...'
!
!      CALL CGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD,  &
!                   M,  N, ZXsp, LDX, ZYsp, LDY, NRNK, TOLsp,  &
!                   K, ZEIGSsp, ZZsp, LDZ,  RESsp,  &
!                   ZAUsp, LDAU, ZWsp,  LDW,   ZSsp, LDS,        &
!                   ZDUMMYsp, -1, WDUMMYsp, LWORK, IDUMMY, LIWORK, INFO )
!      WRITE(*,*) 'Done! INFO = ', INFO
!
!     LZWORK = INT(ZDUMMYsp(LWMINOPT))
!     ALLOCATE(ZWORKsp(LZWORK))
!     LIWORK = IDUMMY(1)
!     ALLOCATE(IWORK(LIWORK))
!     LWORK = INT(WDUMMYsp(1))
!     ALLOCATE(SWORK(LWORK))
!
! !     ALLOCATE( ZAUsp(LDAU,N) )
! !     ALLOCATE( ZWsp(LDW,N) )
! !     ALLOCATE( ZSsp(LDS,N) )
! !     ALLOCATE( ZZsp(LDZ,N) )
! !     ALLOCATE( RESsp(N) )
! !     ALLOCATE( ZEIGSsp(N) )
!
!
!     WRITE(*,*) 'Run CGEDMD ....'
!
!     CALL CGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD,  &
!                   M,  N, ZXsp, LDX, ZYsp, LDY, NRNK, TOLsp,  &
!                   K, ZEIGSsp, ZZsp, LDZ,  RESsp,  &
!                   ZAUsp, LDAU, ZWsp,  LDW,   ZSsp, LDS,        &
!                   ZWORKsp, LZWORK, SWORK, LWORK, IWORK, LIWORK, INFO )
!     WRITE(*,*) 'Done! INFO = ', INFO
!
!    IF ( MANUAL_TEST ) THEN
!
!      CALL matrix_save_Cre( 9, 'CEIGSre.data', K, 1, (ZEIGSsp), K )
!      CALL matrix_save_Cim( 9, 'CEIGSim.data', K, 1, (ZEIGSsp), K )
!      IF ( LSAME(JOBZ,'V') ) THEN
!      CALL matrix_save_Cre( 9, 'CEVECre.data', M, K, (ZZsp), LDZ )
!      CALL matrix_save_Cim( 9, 'CEVECim.data', M, K, (ZZsp), LDZ )
!      END IF
!      IF ( LSAME(RESIDS,'R') ) &
!          CALL matrix_save_sp( 9, 'CRES.data', K, 1, RESsp, K )
!
!      IF ( LSAME(JOBREF,'E') ) THEN
!           ! The unscaled vectors of the Exact DMD are computed.
!           ! This option is included for the sake of completeness,
!           ! for users who prefer the Exact DMD vectors. The
!           ! returned vectors are saved in the real form, and can be
!           ! tested separately using a Matlab script.
!           CALL matrix_save_Cre( 9, 'CEXVECre.exdmd',  M, K, ZAUsp, LDAU )
!           CALL matrix_save_Cim( 9, 'CEXVECim.exdmd',  M, K, ZAUsp, LDAU )
!      END IF
!      IF ( LSAME(JOBREF,'R') ) THEN
!           ! The returned vectors are saved in the real form, and
!           ! can be tested separately using a Matlab script.
!           CALL matrix_save_Cre( 9, 'CREFre.exdmd',  M, K, ZAUsp, LDAU )
!           CALL matrix_save_Cim( 9, 'CREFim.exdmd',  M, K, ZAUsp, LDAU )
!      END IF
!           ! The eigenvectors are returned in factored form.
!      IF ( LSAME(JOBZ,'F') .OR. LSAME(JOBZ,'V')) THEN
!           CALL matrix_save_Cre( 9, 'CWre.dmd',  K, K, ZWsp, LDW )
!           CALL matrix_save_Cim( 9, 'CWim.dmd',  K, K, ZWsp, LDW )
!           CALL matrix_save_Cre( 9, 'CUkre.dmd',  M, K, ZXsp, LDX )
!           CALL matrix_save_Cim( 9, 'CUkim.dmd',  M, K, ZXsp, LDX )
!      END IF
!
!    END IF ! manual test ... .saving the results for case study/analysis
!
!
!    DEALLOCATE(ZWORKsp)
!    DEALLOCATE(SWORK)
!    DEALLOCATE(IWORK)
!
!    IF ( TEST_QRDMD ) THEN
!      ZJOBDATA(8) = 1
!      ZFsp(1:M,1:N+1) = ZF0(1:M,1:N+1)
!      WRITE(*,*) 'Call CGEDMDQ for workspace length ...'
!      CALL CGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, &
!                    WHTSVD, M, N+1, ZFsp, LDF,  ZXsp, LDX,  ZYsp, LDY,  &
!                    NRNK,  TOLsp, K, ZEIGSsp, ZZsp, LDZ, RESsp,  ZAUsp,  &
!                    LDAU, ZWsp, LDW, ZSsp, LDS, ZDUMMYsp, -1,   &
!                    WDUMMYsp,  -1, IDUMMY, -1, INFO )
!      WRITE(*,*) 'Done! INFO = ', INFO
!
!     LZWORK = INT(ZDUMMYsp(LWMINOPT))
!     ALLOCATE(ZWORKsp(LZWORK))
!     LIWORK = IDUMMY(1)
!     ALLOCATE(IWORK(LIWORK))
!     LWORK = INT(WDUMMYsp(1))
!     ALLOCATE(SWORK(LWORK))
!
!     CALL CGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, &
!                    WHTSVD, M, N+1, ZFsp, LDF,  ZXsp, LDX,  ZYsp, LDY,  &
!                    NRNK,  TOLsp, KQ, ZEIGSsp, ZZsp, LDZ, RESsp,  ZAUsp,  &
!                    LDAU, ZWsp, LDW, ZSsp, LDS, ZWORKsp, LZWORK,   &
!                    SWORK,  LWORK, IWORK, LIWORK, INFO )
!      WRITE(*,*) 'Done! INFO = ', INFO
!
!      IF ( MANUAL_TEST ) THEN
!      CALL matrix_save_Cre( 9, 'CQEIGSre.data', KQ, 1, (ZEIGSsp), KQ )
!      CALL matrix_save_Cim( 9, 'CQEIGSim.data', KQ, 1, (ZEIGSsp), KQ )
!      IF ( LSAME(JOBZ,'V') ) THEN
!      CALL matrix_save_Cre( 9, 'CQEVECre.data', M, KQ, (ZZsp), LDZ )
!      CALL matrix_save_Cim( 9, 'CQEVECim.data', M, KQ, (ZZsp), LDZ )
!      END IF
!      IF ( LSAME(RESIDS,'R') ) &
!          CALL matrix_save_sp( 9, 'CQRES.data', KQ, 1, RESsp, KQ )
!
!      IF ( LSAME(RESIDS, 'R') ) THEN
!       CALL matrix_save_sp( 9,'QCRES.zgedmdq', KQ, 1,   RESsp,    KQ )
!       END IF
!
!       IF ( LSAME(WANTQ,'Q') ) THEN
!         CALL matrix_save_Cre( 9, 'QCQQre.data', M, MIN(M,N+1), ZFsp, LDF )
!         CALL matrix_save_Cim( 9, 'QCQQim.data', M, MIN(M,N+1), ZFsp, LDF )
!       END IF
!
!       IF ( LSAME(JOBREF,'E') ) THEN
!           ! The unscaled vectors of the Exact DMD are computed.
!           ! This option is included for the sake of completeness,
!           ! for users who prefer the Exact DMD vectors. The
!           ! returned vectors are saved in the real form, and can be
!           ! tested separately using a Matlab script.
!           CALL matrix_save_Cre( 9, 'QCEXVECre.exdmd',  MIN(M,N+1), KQ, ZAUsp, LDAU )
!           CALL matrix_save_Cim( 9, 'QCEXVECim.exdmd',  MIN(M,N+1), KQ, ZAUsp, LDAU )
!       ELSE IF ( LSAME(JOBREF,'R') ) THEN
!           CALL matrix_save_Cre( 9, 'QCAUre.data',  N+1, KQ, ZAUsp, LDAU )
!           CALL matrix_save_Cim( 9, 'QCAUim.data',  N+1, KQ, ZAUsp, LDAU )
!           CALL matrix_save_Cre( 9, 'QCXUre.qrdmd', (N+1), KQ,  ZXsp, LDX )
!           CALL matrix_save_Cim( 9, 'QCXUim.qrdmd', (N+1), KQ,  ZXsp, LDX )
!       END IF
!
!       IF ( LSAME(JOBZ,'F') .OR. LSAME(JOBZ,'V')) THEN
!           CALL matrix_save_Cre( 9, 'QCWre.dmd',  K, K, ZWsp, LDW )
!           CALL matrix_save_Cim( 9, 'QCWim.dmd',  K, K, ZWsp, LDW )
!           CALL matrix_save_Cre( 9, 'QCUkre.dmd',  M, K, ZXsp, LDX )
!           CALL matrix_save_Cim( 9, 'QCUkim.dmd',  M, K, ZXsp, LDX )
!       END IF
!      END IF
!      DEALLOCATE(IWORK)
!      DEALLOCATE(SWORK)
!      DEALLOCATE(ZWORKsp)
!    END IF ! TEST_QRDMD SINGLE COMPLEX
!
!
!      END IF ! TEST_SP
!     IF ( MANUAL_TEST ) THEN
!         OPEN( UNIT   = 9,         FILE   = 'ZJOBDATA.data', &
!               STATUS = 'replace', ACTION = 'write' )
!          DO i = 1, 8
!             WRITE(9, FMT="(I4)", advance="NO") ZJOBDATA(i)
!          END DO
!          CLOSE(9)
!          WRITE(*,*) 'Output saved to files.'
!     END IF
!      END DO
!      END DO
!      END DO
!      END DO
!      END DO
!      END DO

END IF ! REAL OR COMPLEX DATA



    STOP

    END

