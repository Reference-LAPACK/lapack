      SUBROUTINE ZLA_HERFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA,
     $                                AF, LDAF, IPIV, COLEQU, C, B, LDB,
     $                                Y, LDY, BERR_OUT, N_NORMS, ERRS_N,
     $                                ERRS_C, RES, AYB, DY, Y_TAIL,
     $                                RCOND, ITHRESH, RTHRESH, DZ_UB,
     $                                IGNORE_CWISE, INFO )
*
*     -- LAPACK routine (version 3.2)                                 --
*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and --
*     -- Jason Riedy of Univ. of California Berkeley.                 --
*     -- November 2008                                                --
*
*     -- LAPACK is a software package provided by Univ. of Tennessee, --
*     -- Univ. of California Berkeley and NAG Ltd.                    --
*
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE,
     $                   N_NORMS, ITHRESH
      CHARACTER          UPLO
      LOGICAL            COLEQU, IGNORE_CWISE
      DOUBLE PRECISION   RTHRESH, DZ_UB
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )
      DOUBLE PRECISION   C( * ), AYB( * ), RCOND, BERR_OUT( * ),
     $                   ERRS_N( NRHS, * ), ERRS_C( NRHS, * )
*     ..
*
*  Purpose
*  =======
* 
*  ZLA_HERFSX_EXTENDED computes ... .
*
*  Arguments
*  =========
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            UPLO2, CNT, I, J, X_STATE, Z_STATE,
     $                   Y_PREC_STATE
      DOUBLE PRECISION   YK, DYK, YMIN, NORMY, NORMX, NORMDX, DXRAT,
     $                   DZRAT, PREVNORMDX, PREV_DZ_Z, DXRATMAX,
     $                   DZRATMAX, DX_X, DZ_Z, FINAL_DX_X, FINAL_DZ_Z,
     $                   EPS, HUGEVAL, INCR_THRESH
      LOGICAL            INCR_PREC
      COMPLEX*16         ZDUM
*     ..
*     .. Parameters ..
      INTEGER            UNSTABLE_STATE, WORKING_STATE, CONV_STATE,
     $                   NOPROG_STATE, BASE_RESIDUAL, EXTRA_RESIDUAL,
     $                   EXTRA_Y
      PARAMETER          ( UNSTABLE_STATE = 0, WORKING_STATE = 1,
     $                   CONV_STATE = 2, NOPROG_STATE = 3 )
      PARAMETER          ( BASE_RESIDUAL = 0, EXTRA_RESIDUAL = 1,
     $                   EXTRA_Y = 2 )
      INTEGER            FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I
      INTEGER            RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I
      INTEGER            CMP_ERR_I, PIV_GROWTH_I
      PARAMETER          ( FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2,
     $                   BERR_I = 3 )
      PARAMETER          ( RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 )
      PARAMETER          ( CMP_RCOND_I = 7, CMP_ERR_I = 8,
     $                   PIV_GROWTH_I = 9 )
      INTEGER            LA_LINRX_ITREF_I, LA_LINRX_ITHRESH_I,
     $                   LA_LINRX_CWISE_I
      PARAMETER          ( LA_LINRX_ITREF_I = 1,
     $                   LA_LINRX_ITHRESH_I = 2 )
      PARAMETER          ( LA_LINRX_CWISE_I = 3 )
      INTEGER            LA_LINRX_TRUST_I, LA_LINRX_ERR_I,
     $                   LA_LINRX_RCOND_I
      PARAMETER          ( LA_LINRX_TRUST_I = 1, LA_LINRX_ERR_I = 2 )
      PARAMETER          ( LA_LINRX_RCOND_I = 3 )
      INTEGER            LA_LINRX_MAX_N_ERRS
      PARAMETER          ( LA_LINRX_MAX_N_ERRS = 3 )      
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           ILAUPLO
      INTEGER            ILAUPLO
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZCOPY, ZHETRS, ZHEMV, BLAS_ZHEMV_X,
     $                   BLAS_ZHEMV2_X, ZLA_HEAMV, ZLA_WWADDW,
     $                   ZLA_LIN_BERR
      DOUBLE PRECISION   DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, DIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      IF (INFO.NE.0) RETURN
      EPS = DLAMCH( 'Epsilon' )
      HUGEVAL = DLAMCH( 'Overflow' )
*     Force HUGEVAL to Inf
      HUGEVAL = HUGEVAL * HUGEVAL
*     Using HUGEVAL may lead to spurious underflows.
      INCR_THRESH = DBLE( N ) * EPS

      IF ( LSAME ( UPLO, 'L' ) ) THEN
         UPLO2 = ILAUPLO( 'L' )
      ELSE
         UPLO2 = ILAUPLO( 'U' )
      ENDIF

      DO J = 1, NRHS
         Y_PREC_STATE = EXTRA_RESIDUAL
         IF ( Y_PREC_STATE .EQ. EXTRA_Y ) THEN
            DO I = 1, N
               Y_TAIL( I ) = 0.0D+0
            END DO
         END IF

         DXRAT = 0.0D+0
         DXRATMAX = 0.0D+0
         DZRAT = 0.0D+0
         DZRATMAX = 0.0D+0
         FINAL_DX_X = HUGEVAL
         FINAL_DZ_Z = HUGEVAL
         PREVNORMDX = HUGEVAL
         PREV_DZ_Z = HUGEVAL
         DZ_Z = HUGEVAL
         DX_X = HUGEVAL

         X_STATE = WORKING_STATE
         Z_STATE = UNSTABLE_STATE
         INCR_PREC = .FALSE.

         DO CNT = 1, ITHRESH
*
*         Compute residual RES = B_s - op(A_s) * Y,
*             op(A) = A, A**T, or A**H depending on TRANS (and type).
*
            CALL ZCOPY( N, B( 1, J ), 1, RES, 1 )
            IF ( Y_PREC_STATE .EQ. BASE_RESIDUAL ) THEN
               CALL ZHEMV( UPLO, N, DCMPLX(-1.0D+0), A, LDA, Y( 1, J ), 
     $              1, DCMPLX(1.0D+0), RES, 1 )
            ELSE IF ( Y_PREC_STATE .EQ. EXTRA_RESIDUAL ) THEN
               CALL BLAS_ZHEMV_X( UPLO2, N, DCMPLX(-1.0D+0), A, LDA,
     $              Y( 1, J ), 1, DCMPLX(1.0D+0), RES, 1, PREC_TYPE)
            ELSE
               CALL BLAS_ZHEMV2_X(UPLO2, N, DCMPLX(-1.0D+0), A, LDA,
     $              Y(1, J), Y_TAIL, 1, DCMPLX(1.0D+0), RES, 1, 
     $     PREC_TYPE)
            END IF

!         XXX: RES is no longer needed.
            CALL ZCOPY( N, RES, 1, DY, 1 )
            CALL ZHETRS( UPLO, N, NRHS, AF, LDAF, IPIV, DY, N, INFO )
*
*         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT.
*
            NORMX = 0.0D+0
            NORMY = 0.0D+0
            NORMDX = 0.0D+0
            DZ_Z = 0.0D+0
            YMIN = HUGEVAL

            DO I = 1, N
               YK = CABS1( Y( I, J ) )
               DYK = CABS1( DY( I ) )

               IF (YK .NE. 0.0D+0) THEN
                  DZ_Z = MAX( DZ_Z, DYK / YK )
               ELSE IF ( DYK .NE. 0.0D+0 ) THEN
                  DZ_Z = HUGEVAL
               END IF

               YMIN = MIN( YMIN, YK )

               NORMY = MAX( NORMY, YK )

               IF ( COLEQU ) THEN
                  NORMX = MAX( NORMX, YK * C( I ) )
                  NORMDX = MAX( NORMDX, DYK * C( I ) )
               ELSE
                  NORMX = NORMY
                  NORMDX = MAX( NORMDX, DYK )
               END IF
            END DO

            IF ( NORMX .NE. 0.0D+0 ) THEN
               DX_X = NORMDX / NORMX
            ELSE IF ( NORMDX .EQ. 0.0D+0 ) THEN
               DX_X = 0.0D+0
            ELSE
               DX_X = HUGEVAL
            END IF

            DXRAT = NORMDX / PREVNORMDX
            DZRAT = DZ_Z / PREV_DZ_Z
*
*         Check termination criteria.
*
            IF ( YMIN*RCOND .LT. INCR_THRESH*NORMY
     $           .AND. Y_PREC_STATE .LT. EXTRA_Y )
     $           INCR_PREC = .TRUE.

            IF ( X_STATE .EQ. NOPROG_STATE .AND. DXRAT .LE. RTHRESH )
     $           X_STATE = WORKING_STATE
            IF ( X_STATE .EQ. WORKING_STATE ) THEN
               IF ( DX_X .LE. EPS ) THEN
                  X_STATE = CONV_STATE
               ELSE IF ( DXRAT .GT. RTHRESH ) THEN
                  IF ( Y_PREC_STATE .NE. EXTRA_Y ) THEN
                     INCR_PREC = .TRUE.
                  ELSE
                     X_STATE = NOPROG_STATE
                  END IF
               ELSE
                  IF (DXRAT .GT. DXRATMAX) DXRATMAX = DXRAT
               END IF
               IF ( X_STATE .GT. WORKING_STATE ) FINAL_DX_X = DX_X
            END IF

            IF ( Z_STATE .EQ. UNSTABLE_STATE .AND. DZ_Z .LE. DZ_UB )
     $           Z_STATE = WORKING_STATE
            IF ( Z_STATE .EQ. NOPROG_STATE .AND. DZRAT .LE. RTHRESH )
     $           Z_STATE = WORKING_STATE
            IF ( Z_STATE .EQ. WORKING_STATE ) THEN
               IF ( DZ_Z .LE. EPS ) THEN
                  Z_STATE = CONV_STATE
               ELSE IF ( DZ_Z .GT. DZ_UB ) THEN
                  Z_STATE = UNSTABLE_STATE
                  DZRATMAX = 0.0D+0
                  FINAL_DZ_Z = HUGEVAL
               ELSE IF ( DZRAT .GT. RTHRESH ) THEN
                  IF ( Y_PREC_STATE .NE. EXTRA_Y ) THEN
                     INCR_PREC = .TRUE.
                  ELSE
                     Z_STATE = NOPROG_STATE
                  END IF
               ELSE
                  IF ( DZRAT .GT. DZRATMAX ) DZRATMAX = DZRAT
               END IF
               IF ( Z_STATE .GT. WORKING_STATE ) FINAL_DZ_Z = DZ_Z
            END IF

            IF ( X_STATE.NE.WORKING_STATE.AND.
     $           ( IGNORE_CWISE.OR.Z_STATE.NE.WORKING_STATE ) )
     $           GOTO 666

            IF ( INCR_PREC ) THEN
               INCR_PREC = .FALSE.
               Y_PREC_STATE = Y_PREC_STATE + 1
               DO I = 1, N
                  Y_TAIL( I ) = 0.0D+0
               END DO
            END IF

            PREVNORMDX = NORMDX
            PREV_DZ_Z = DZ_Z
*
*           Update soluton.
*
            IF ( Y_PREC_STATE .LT. EXTRA_Y ) THEN
               CALL ZAXPY( N, DCMPLX(1.0D+0), DY, 1, Y(1,J), 1 )
            ELSE
               CALL ZLA_WWADDW( N, Y(1,J), Y_TAIL, DY )
            END IF

         END DO
*        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT.
 666     CONTINUE
*
*     Set final_* when cnt hits ithresh.
*
         IF ( X_STATE .EQ. WORKING_STATE ) FINAL_DX_X = DX_X
         IF ( Z_STATE .EQ. WORKING_STATE ) FINAL_DZ_Z = DZ_Z
*
*     Compute error bounds.
*
         IF ( N_NORMS .GE. 1 ) THEN
            ERRS_N( J, LA_LINRX_ERR_I ) = FINAL_DX_X / (1 - DXRATMAX)
         END IF
         IF (N_NORMS .GE. 2) THEN
            ERRS_C( J, LA_LINRX_ERR_I ) = FINAL_DZ_Z / (1 - DZRATMAX)
         END IF
*
*     Compute componentwise relative backward error from formula
*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
*     where abs(Z) is the componentwise absolute value of the matrix
*     or vector Z.
*
*         Compute residual RES = B_s - op(A_s) * Y,
*             op(A) = A, A**T, or A**H depending on TRANS (and type).
*
         CALL ZCOPY( N, B( 1, J ), 1, RES, 1 )
         CALL ZHEMV( UPLO, N, DCMPLX(-1.0D+0), A, LDA, Y(1,J), 1,
     $        DCMPLX(1.0D+0), RES, 1 )

         DO I = 1, N
            AYB( I ) = CABS1( B( I, J ) )
         END DO
*
*     Compute abs(op(A_s))*abs(Y) + abs(B_s).
*
         CALL ZLA_HEAMV( UPLO2, N, 1.0D+0,
     $        A, LDA, Y(1, J), 1, 1.0D+0, AYB, 1 )

         CALL ZLA_LIN_BERR( N, N, 1, RES, AYB, BERR_OUT( J ) )
*
*     End of loop for each RHS.
*
      END DO
*
      RETURN
      END
