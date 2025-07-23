*> \brief \b DKGT01
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DKGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D,
*                          WORK, RESULT )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            ITYPE, LDA, LDB, LDZ, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), D( * ), RESULT( * ),
*      $                   WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DKGT01 checks a decomposition of the form
*>
*>    A Z   =  B Z D or
*>    A B Z =  Z D or
*>    B A Z =  Z D
*>
*> where A is a skew-symmetric matrix, B is
*> skew-symmetric positive definite, Z is orthogonal, and D is diagonal.
*>
*> One of the following test ratios is computed:
*>
*> ITYPE = 1:  RESULT(1) = | A Z - B Z D | / ( |A| |Z| n ulp )
*>
*> ITYPE = 2:  RESULT(1) = | A B Z - Z D | / ( |A| |Z| n ulp )
*>
*> ITYPE = 3:  RESULT(1) = | B A Z - Z D | / ( |A| |Z| n ulp )
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ITYPE
*> \verbatim
*>          ITYPE is INTEGER
*>          The form of the skew-symmetric generalized eigenproblem.
*>          = 1:  A*z = (lambda)*B*z
*>          = 2:  A*B*z = (lambda)*z
*>          = 3:  B*A*z = (lambda)*z
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          skew-symmetric matrices A and B is stored.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of eigenvalues found.  0 <= M <= N.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          The original skew-symmetric matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB, N)
*>          The original symmetric positive definite matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, M)
*>          The computed eigenvectors of the generalized eigenproblem.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (M)
*>          The computed eigenvalues of the generalized eigenproblem.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N*N)
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is DOUBLE PRECISION array, dimension (1)
*>          The test ratio as described above.
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
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE DKGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D,
     $                   WORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            ITYPE, LDA, LDB, LDZ, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), D( * ), RESULT( * ),
     $                   WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, ULP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANKY
      EXTERNAL           DLAMCH, DLANGE, DLANKY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DSCAL, DAXPY, DSYMM, DKYMM
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      ULP = DLAMCH( 'Epsilon' )
*
*     Compute product of 1-norms of A and Z.
*
      ANORM = DLANKY( '1', UPLO, N, A, LDA, WORK )*
     $        DLANGE( '1', N, M, Z, LDZ, WORK )
      IF( ANORM.EQ.ZERO )
     $   ANORM = ONE
*
      IF( ITYPE.EQ.1 ) THEN
*
*        Norm of AZ - BZD
*
         CALL DKYMM( 'Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO,
     $               WORK, N )
         DO 10 I = 1, M-1
            CALL DCOPY( N, Z( 1, I+1 ), 1, WORK(N**2+(I-1)*N+1), 1 )
            CALL DSCAL( N, D( I ), WORK(N**2+(I-1)*N+1), 1 )
   10    CONTINUE
         DO 20 I = 2, M-1
            CALL DAXPY( N, -D( I-1 ), Z( 1, I-1 ), 1,
     $                  WORK(N**2+(I-1)*N+1), 1 )
   20    CONTINUE
         CALL DCOPY( N, Z( 1, M-1 ), 1, WORK(N**2+(M-1)*N+1), 1 )
         CALL DSCAL( N, -D( M-1 ), WORK(N**2+(M-1)*N+1), 1 )
         CALL DSYMM( 'Left', UPLO, N, M, ONE, B, LDB, WORK(N**2+1),
     $               N, -ONE, WORK, N )
*
         RESULT( 1 ) = ( DLANGE( '1', N, M, WORK, N, WORK ) / ANORM ) /
     $                 ( N*ULP )
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Norm of ABZ - ZD
*
         CALL DSYMM( 'Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, ZERO,
     $               WORK, N )
         DO 30 I = 1, M-1
            CALL DCOPY( N, Z( 1, I+1 ), 1, WORK(N**2+(I-1)*N+1), 1 )
            CALL DSCAL( N, D( I ), WORK(N**2+(I-1)*N+1), 1 )
   30    CONTINUE
         DO 40 I = 2, M-1
            CALL DAXPY( N, -D( I-1 ), Z( 1, I-1 ), 1,
     $                  WORK(N**2+(I-1)*N+1), 1 )
   40    CONTINUE
         CALL DCOPY( N, Z( 1, M-1 ), 1, WORK(N**2+(M-1)*N+1), 1 )
         CALL DSCAL( N, -D( M-1 ), WORK(N**2+(M-1)*N+1), 1 )
         CALL DKYMM( 'Left', UPLO, N, M, ONE, A, LDA, WORK, N, -ONE,
     $               WORK(N**2+1), N )
*
         RESULT( 1 ) = ( DLANGE( '1', N, M, WORK(N**2+1), N, WORK )
     $                 / ANORM ) / ( N*ULP )
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Norm of BAZ - ZD
*
         CALL DKYMM( 'Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO,
     $               WORK, N )
         DO 50 I = 1, M-1
            CALL DCOPY( N, Z( 1, I+1 ), 1, WORK(N**2+(I-1)*N+1), 1 )
            CALL DSCAL( N, D( I ), WORK(N**2+(I-1)*N+1), 1 )
   50    CONTINUE
         DO 60 I = 2, M-1
            CALL DAXPY( N, -D( I-1 ), Z( 1, I-1 ), 1,
     $                  WORK(N**2+(I-1)*N+1), 1 )
   60    CONTINUE
         CALL DCOPY( N, Z( 1, M-1 ), 1, WORK(N**2+(M-1)*N+1), 1 )
         CALL DSCAL( N, -D( M-1 ), WORK(N**2+(M-1)*N+1), 1 )
         CALL DSYMM( 'Left', UPLO, N, M, ONE, B, LDB, WORK, N, -ONE,
     $               WORK(N**2+1), N )
*
         RESULT( 1 ) = ( DLANGE( '1', N, M, WORK(N**2+1), N, WORK )
     $                 / ANORM ) / ( N*ULP )
      END IF
*
      RETURN
*
*     End of DKGT01
*
      END
