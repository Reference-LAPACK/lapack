! This is a test program for checking the implementations of
! the implementations of the following subroutines
!
! ZGEDMD,  for computation of the
!          Dynamic Mode Decomposition (DMD)
! ZGEDMDQ, for computation of a
!          QR factorization based compressed DMD
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
      use iso_fortran_env, only: real32
      IMPLICIT NONE
      integer, parameter :: WP = real32

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
      REAL(KIND=WP), ALLOCATABLE, DIMENSION(:,:) ::          &
                     Asp, AUsp, Fsp, F1sp, LAMBDAsp,         &
                     LAMBDAQsp, Ssp, Wsp, Xsp, Ysp,  Zsp
      REAL(KIND=WP), ALLOCATABLE, DIMENSION(:)   ::          &
                     IEIGsp, IEIGQsp, REIGsp, REIGQsp, RESsp,&
                     SWORK
      REAL(KIND=WP) :: WDUMMYsp(2)
      REAL(KIND=WP) :: EPSsp, TOLsp
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
      COMPLEX(KIND=WP), ALLOCATABLE, DIMENSION(:,:) ::       &
                 ZAUsp, ZFsp, ZSsp, ZWsp, ZXsp, ZYsp, ZZsp
      COMPLEX(KIND=WP), ALLOCATABLE, DIMENSION(:)   ::       &
                        ZEIGSsp, ZWORKsp
      COMPLEX(KIND=WP) ZDUMMYsp(22)
!............................................................
      INTEGER :: K, KQ, LDF, LDS, LDA, LDAU, LDW, LDX, LDY,  &
                 LDZ, LIWORK, LWORK, M, N, L, LLOOP, NRNK, NRNKsp
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

!..... external functions (BLAS and LAPACK)
      EXTERNAL         DLAMCH, DLANGE, DNRM2, SLAMCH
      REAL(KIND=WP) :: DLAMCH, DLANGE, DNRM2
      REAL(KIND=WP) :: SLAMCH
      EXTERNAL IZAMAX
      INTEGER  IZAMAX
      EXTERNAL         LSAME
      LOGICAL          LSAME

      INTRINSIC ABS, INT, MIN, MAX
!............................................................

      EPS   = DLAMCH( 'P' )  ! machine precision DP
      EPSsp = SLAMCH( 'P' )  ! machine precision SP

         REAL_DATA = .FALSE.
          KEYBIN       = .FALSE.
          VISUAL_CHECK = .FALSE.
          MANUAL_TEST  = .FALSE.

      TEST_SP = .TRUE.
      TESTQ_SP = .TRUE.

      DO K = 1, 2

      IF ( K == 1) THEN
          TWO_TRAJ = .FALSE.
      ELSE IF ( K == 2 ) THEN
          TWO_TRAJ = .TRUE.
      END IF

      IF ( TWO_TRAJ ) THEN
          TEST_QRDMD = .FALSE.
      ELSE                    ! This version (ZGEDMDQ) is for
          TEST_QRDMD = .TRUE. ! one trajectory stream of data.
      END IF

      ! Set the dimensions of the problem ...
      !READ( NIN, FMT = * )( IOLDSD( I ), I = 1, 4 )
      READ(*,*) L
      WRITE(*,*) 'L = ', L
      DO LLOOP = 1, L
      WRITE(*,*) 'L Loop Index = ', LLOOP
      READ(*,*) M
      WRITE(*,*) 'M = ', M
!     ... and the number of snapshots.
      READ(*,*) N
      WRITE(*,*) 'N = ', N
      IF ( ( MIN(M,N) == 0 ) .OR. ( M < N )  ) THEN
          WRITE(*,*) 'Bad dimensions. Required: M >= N > 0.'
          STOP
      END IF

      ISEED(1) = 4
      ISEED(2) = 3
      ISEED(3) = 2
      ISEED(4) = 1

      ! Loop over all parameter MODE values for DLATMR (+-1,..,+-6)
      DO MODE = -6, 6, 1
!..............................................................

      LDA = M
      LDF = M
      LDX = MAX(M,N+1)
      LDY = MAX(M,N+1)
      LDW = N
      LDZ = M
      LDAU = MAX(M,N+1)
      LDS = N

!============================================================================
!============================================================================
! Test complex data subroutines
!============================================================================
!============================================================================

      COND = 1.0D8
      ZMAX = (1.0D2,1.0D2)
      RSIGN = 'F'
      GRADE = 'N'
      MODEL = 6
      CONDL = 1.0D2
      MODER = 6
      CONDR = 1.0D2
      PIVTNG = 'N'
      ALLOCATE( IWORK(2*M) )
      ALLOCATE( ZDA(M) )
      ALLOCATE( ZA(LDA,M) )
      ALLOCATE( ZAC(LDA,M) )
      ALLOCATE(ZDL(M))
      ALLOCATE(ZDR(M))
      CALL CLATMR( M, M, 'S', ISEED, 'N', ZDA, MODE, COND, &
                   ZMAX, RSIGN, GRADE, ZDL, MODEL,  CONDL, &
                   ZDR, MODER, CONDR, PIVTNG, IWORK, M, M, &
                   ZERO, -ONE, 'N', ZA, LDA, IWORK(M+1), INFO )
      DEALLOCATE(IWORK)


      !'Compute and store the eigenvalues &
      !           &and eigenvectors of A (ZGEEV) . ..'

      LZWORK = MAX(1,2*M)
      ALLOCATE( ZEIGSA(M) )
      ALLOCATE(ZVA(LDA,M))
      ALLOCATE( ZWORK(LZWORK) )
      ALLOCATE( WORK(2*M) )
      ZAC = ZA
      CALL CGEEV( 'N','V', M, ZAC, LDA, ZEIGSA, ZVA, M, &
                  ZVA, M, ZWORK, LZWORK, WORK, INFO ) ! LAPACK CALL
      DEALLOCATE(WORK)
      DEALLOCATE(ZWORK)

      TMP = ABS(ZEIGSA(IZAMAX(M, ZEIGSA, 1))) ! The spectral radius of ZA

      CALL CLASCL( 'G',0, 0, TMP, ONE, M, M, &
                   ZA, LDA, INFO )
      CALL CLASCL( 'G',0, 0, TMP, ONE, M, 1, &
                   ZEIGSA, M, INFO )

      ALLOCATE( ZF(LDF,N+1) )
      ALLOCATE( ZF0(LDF,N+1) )
      CALL CLARNV(2, ISEED, M, ZF(1,1) )
      DO i = 1, N
         CALL CGEMV( 'N', M, M, ZONE, ZA, M, ZF(1,i), 1, ZZERO, &
              ZF(1,i+1), 1 )
      END DO
      ZF0(1:M,1:N+1) = ZF(1:M,1:N+1)
!........................................................................

      ALLOCATE( ZX(LDX,N) )
      ALLOCATE( ZY(LDY,N+1) )
      ALLOCATE( ZAU(LDAU,N) )
      ALLOCATE( ZW(LDW,N) )
      ALLOCATE( ZS(LDS,N) )
      ALLOCATE( ZZ(LDZ,N) )
      ALLOCATE( RES(N) )
      ALLOCATE( ZEIGS(N) )

      ZX(1:M,1:N) = ZF0(1:M,1:N)
      ZY(1:M,1:N) = ZF0(1:M,2:N+1)

      IF ( TEST_SP ) THEN
          ALLOCATE( ZXsp(LDX,N) )
          ALLOCATE( ZYsp(LDY,N+1) )
          ALLOCATE( ZFsp(LDF,N+1) )
          ALLOCATE( ZEIGSsp(N))
          ALLOCATE( RESsp(N) )
          ALLOCATE( ZWsp(LDW,N) )
          ALLOCATE( ZAUsp(LDAU,N) )
          ALLOCATE( ZZsp(LDZ,N) )
          ALLOCATE( ZSsp(LDS,N) )
          ZXsp(1:M,1:N) = ZX(1:M,1:N)
          ZYsp(1:M,1:N) = ZY(1:M,1:N)
          ZFsp(1:M,1:N+1) = ZF(1:M,1:N+1)
      END IF

      TOL   = M*EPS
      TOLsp = M*EPSsp
      ! here nested doo loops

            DO iJOBZ = 1, 4

          SELECT CASE ( iJOBZ )
          CASE(1)
              JOBZ   = 'V'
              RESIDS = 'R'
          CASE(2)
              JOBZ   = 'V'
              RESIDS = 'N'
          CASE(3)
              JOBZ   = 'F'
              RESIDS = 'N'
          CASE(4)
              JOBZ   = 'N'
              RESIDS = 'N'
          END SELECT

      DO iJOBREF = 1, 3

          SELECT CASE ( iJOBREF )
          CASE(1)
              JOBREF = 'R'
          CASE(2)
              JOBREF = 'E'
          CASE(3)
              JOBREF = 'N'
          END SELECT

      DO iSCALE = 1, 4

          SELECT CASE ( iSCALE )
          CASE(1)
              SCALE = 'S'
          CASE(2)
              SCALE = 'C'
          CASE(3)
              SCALE = 'Y'
          CASE(4)
              SCALE = 'N'
          END SELECT

      DO iNRNK = -1, -2, -1
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

      ZX(1:M,1:N) = ZF0(1:M,1:N)
      ZY(1:M,1:N) = ZF0(1:M,2:N+1)
      ZF(1:M,1:N+1) = ZF0(1:M,1:N+1)

      CALL CGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD,  &
                   M,  N, ZX, LDX, ZY, LDY, NRNK, TOL,  &
                   K, ZEIGS, ZZ, LDZ,  RES,  &
                   ZAU, LDAU, ZW,  LDW,   ZS, LDS,        &
                   ZDUMMY, -1, WDUMMY, LWORK, IDUMMY, LIWORK, INFO )

     LZWORK = INT(ZDUMMY(LWMINOPT))
     ALLOCATE(ZWORK(LZWORK))
     LIWORK = IDUMMY(1)
     ALLOCATE(IWORK(LIWORK))
     LWORK = INT(WDUMMY(1))
     ALLOCATE(WORK(LWORK))

     CALL CGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD,   &
                   M,  N, ZX, LDX, ZY, LDY, NRNK, TOL,   &
                   K, ZEIGS, ZZ, LDZ,  RES, ZAU, LDAU,   &
                   ZW,  LDW,   ZS, LDS,ZWORK, LZWORK,    &
                   WORK, LWORK, IWORK, LIWORK, INFO )

      DEALLOCATE( ZWORK )
      DEALLOCATE( WORK  )
      DEALLOCATE( IWORK )

      IF ( TEST_QRDMD ) THEN
         ZJOBDATA(6) = 1
         WANTQ = 'Q'
         WANTR = 'R'

         CALL CGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, &
                       WHTSVD, M, N+1, ZF, LDF,  ZX, LDX,  ZY, LDY,  &
                       NRNK,  TOL, K, ZEIGS, ZZ, LDZ, RES,  ZAU,  &
                       LDAU, ZW, LDW, ZS, LDS, ZDUMMY, LZWORK,   &
                       WDUMMY,  -1, IDUMMY, -1, INFO )

         LZWORK = INT(ZDUMMY(LWMINOPT))
         ALLOCATE( ZWORK(LZWORK) )
         LIWORK = IDUMMY(1)
         ALLOCATE(IWORK(LIWORK))
         LWORK = INT(WDUMMY(1))
         ALLOCATE(WORK(LWORK))

         CALL CGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, &
                       WHTSVD, M, N+1, ZF, LDF,  ZX, LDX,  ZY, LDY,  &
                       NRNK,  TOL, KQ, ZEIGS, ZZ, LDZ, RES,  ZAU,  &
                       LDAU, ZW, LDW, ZS, LDS, ZWORK, LZWORK,   &
                       WORK,  LWORK, IWORK, LIWORK, INFO )

         DEALLOCATE(ZWORK)
         DEALLOCATE(IWORK)
         DEALLOCATE(WORK)

      END IF ! TEST_QRDMD

      IF ( TEST_SP ) THEN
         ZJOBDATA(7) = 1
           TOLsp = M*EPSsp
           ZXsp(1:M,1:N) = ZF0(1:M,1:N)
           ZYsp(1:M,1:N) = ZF0(1:M,2:N+1)

        CALL CGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD,   &
                     M,  N, ZXsp, LDX, ZYsp, LDY, NRNK,     &
                     TOLsp, K, ZEIGSsp, ZZsp, LDZ,  RESsp,  &
                     ZAUsp, LDAU, ZWsp,  LDW,   ZSsp, LDS,  &
                     ZDUMMYsp, -1, WDUMMYsp, LWORK, IDUMMY, &
                     LIWORK, INFO)

       LZWORK = INT(ZDUMMYsp(LWMINOPT))
       ALLOCATE(ZWORKsp(LZWORK))
       LIWORK = IDUMMY(1)
       ALLOCATE(IWORK(LIWORK))
       LWORK = INT(WDUMMYsp(1))
       ALLOCATE(SWORK(LWORK))

       CALL CGEDMD( SCALE, JOBZ, RESIDS, JOBREF, WHTSVD,  &
                     M,  N, ZXsp, LDX, ZYsp, LDY, NRNK, TOLsp,  &
                     K, ZEIGSsp, ZZsp, LDZ,  RESsp,  &
                     ZAUsp, LDAU, ZWsp,  LDW,   ZSsp, LDS,        &
                     ZWORKsp, LZWORK, SWORK, LWORK,      &
                     IWORK, LIWORK, INFO )

       DEALLOCATE(ZWORKsp)
       DEALLOCATE(SWORK)
       DEALLOCATE(IWORK)

       IF ( TEST_QRDMD ) THEN
       ZJOBDATA(8) = 1
       ZFsp(1:M,1:N+1) = ZF0(1:M,1:N+1)

       CALL CGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, &
                     WHTSVD, M, N+1, ZFsp, LDF,  ZXsp, LDX,  ZYsp,  &
                     LDY, NRNK,  TOLsp, K, ZEIGSsp, ZZsp, LDZ,  &
                     RESsp, ZAUsp, LDAU, ZWsp, LDW, ZSsp, LDS,  &
                     ZDUMMYsp, -1, WDUMMYsp,  -1, IDUMMY, -1, INFO )

      LZWORK = INT(ZDUMMYsp(LWMINOPT))
      ALLOCATE(ZWORKsp(LZWORK))
      LIWORK = IDUMMY(1)
      ALLOCATE(IWORK(LIWORK))
      LWORK = INT(WDUMMYsp(1))
      ALLOCATE(SWORK(LWORK))

      CALL CGEDMDQ( SCALE, JOBZ, RESIDS, WANTQ, WANTR, JOBREF, &
                    WHTSVD, M, N+1, ZFsp, LDF,  ZXsp, LDX,  ZYsp, &
                    LDY, NRNK,  TOLsp, KQ, ZEIGSsp, ZZsp, LDZ, RESsp, &
                    ZAUsp, LDAU, ZWsp, LDW, ZSsp, LDS, ZWORKsp,   &
                    LZWORK, SWORK,  LWORK, IWORK, LIWORK, INFO )

      DEALLOCATE(IWORK)
      DEALLOCATE(SWORK)
      DEALLOCATE(ZWORKsp)
      END IF ! TEST_QRDMD SINGLE COMPLEX


      END IF ! TEST_SP
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

      DEALLOCATE(ZDA)
      DEALLOCATE(ZA)
      DEALLOCATE(ZAC)
      DEALLOCATE(ZDL)
      DEALLOCATE(ZDR)
      DEALLOCATE(ZEIGSA)
      DEALLOCATE(ZVA)
      DEALLOCATE(ZF)
      DEALLOCATE(ZF0)
      DEALLOCATE(ZX)
      DEALLOCATE(ZY)
      DEALLOCATE(ZAU)
      DEALLOCATE(ZW)
      DEALLOCATE(ZS)
      DEALLOCATE(ZZ)
      DEALLOCATE(RES)
      DEALLOCATE(ZEIGS)

      IF ( TEST_SP ) THEN
          DEALLOCATE(ZXsp)
          DEALLOCATE(ZYsp)
          DEALLOCATE(ZFsp)
          DEALLOCATE(ZEIGSsp)
          DEALLOCATE(RESsp)
          DEALLOCATE(ZWsp)
          DEALLOCATE(ZAUsp)
          DEALLOCATE(ZZsp)
          DEALLOCATE(ZSsp)
      END IF

      END DO
      END DO
      END DO

     STOP
     END
