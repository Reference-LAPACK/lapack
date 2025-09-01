*> \brief \b DKTEQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DKTEQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dkteqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dkteqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dkteqr.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DKTEQR( COMPZ, N, E, Z, LDZ, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          COMPZ
*       INTEGER            INFO, LDZ, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   E( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DKTEQR computes all eigenvalues and, optionally, eigenvectors of a
*> skew-symmetric tridiagonal matrix using the implicit double shift
*> QL or QR method.
*> The eigenvectors of a full skew-symmetric matrix can be found if
*> DKYTRD has been used to reduce this matrix to tridiagonal form.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N':  Compute eigenvalues only.
*>          = 'V':  Compute eigenvalues and eigenvectors of the original
*>                  skew-symmetric matrix.  On entry, Z must contain the
*>                  orthogonal matrix used to reduce the original matrix
*>                  to tridiagonal form.
*>          = 'I':  Compute eigenvalues and eigenvectors of the
*>                  tridiagonal matrix.  Z is initialized to the identity
*>                  matrix.
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
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          On entry, the (n-1) lower subdiagonal elements of the
*>          tridiagonal matrix.
*>          On exit, the (n-1) lower subdiagonal elements of the
*>          block diagonal matrix. If INFO = 0, the matrix consists
*>          of 2-by-2 skew-symmetric blocks, and zeros.
*>          The values in E, which represent blocks, are always
*>          positive, and sorted in descending order.
*>          The eigenvalues of each blocks can be evaluated directly.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*>          On entry, if  COMPZ = 'V', then Z contains the orthogonal
*>          matrix used in the reduction to tridiagonal form.
*>          On exit, if INFO = 0, then if  COMPZ = 'V', Z is the
*>          orthogonal matrix transforming the original skew-symmetric
*>          matrix to the block diagonal matrix, and if COMPZ = 'I',
*>          Z is the orthogonal matrix transforming the skew-symmetric
*>          tridiagonal matrix to the block diagonal matrix.
*>          The eigenvectors of corresponding matrix can be evaluated
*>          directly.
*>          If COMPZ = 'N', then Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1, and if
*>          eigenvectors are desired, then  LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (max(1,2*N-4))
*>          If COMPZ = 'N', then WORK is not referenced.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  the algorithm has failed to find all the eigenvalues in
*>                a total of 30*N iterations; if INFO = i, then i
*>                elements of E have not converged to zero; on exit
*>                E contain the elements of a skew-symmetric tridiagonal
*>                matrix which is orthogonally similar to the original
*>                matrix.
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
*> \ingroup steqr
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Derived from subroutine dsteqr.
*>
*> \endverbatim
*
*> \par Contributors:
*  ==================
*>
*>    Shuo Zheng, China, Jul 2025 \n
*>
*  =====================================================================
      SUBROUTINE DKTEQR( COMPZ, N, E, Z, LDZ, WORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM3, LSV, M, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, EPS, EPS2, P, R, VA, VB, E3,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLAPY2, DLANKT
      EXTERNAL           LSAME, DLAMCH, DLAPY2, DLANKT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET,
     $                   DLASRT, DSWAP, DSCAL, DLASR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DKTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0)
     $   RETURN
*
      IF( N.EQ.1) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
      IF( N.EQ.2) THEN
         IF( ICOMPZ.EQ.2 ) THEN
            Z( 1, 1 ) = ONE
            Z( 1, 2 ) = ZERO
            Z( 2, 1 ) = ZERO
            Z( 2, 2 ) = ONE
         END IF
         IF( E(1).LT.ZERO ) THEN
            E(1) = -E(1)
            CALL DSWAP( N, Z( 1, 1 ), 1, Z( 1, 2 ), 1 )
         END IF
         RETURN
      END IF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( ABS( E( M+
     $          1 ) ) )*EPS .AND. M.EQ.L1 ) THEN
               E( M ) = ZERO
               GO TO 30
            ELSEIF( TST.LE.( ABS( E( M-
     $          1 ) ) )*EPS .AND. M.EQ.NM1 ) THEN
               E( M ) = ZERO
               GO TO 30
            ELSEIF( TST.LE.( SQRT( ABS( E( M-1 ) ) )*
     $          SQRT( ABS( E( M+1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANKT( 'M', LEND-L+1, E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
*     Choose between QL and QR iteration
*
	   IF( L.NE.LEND ) THEN
      	IF( ABS( E( LEND-1 ) ).LT.ABS( E( L ) ) ) THEN
         	LEND = LSV
         	L = LENDSV
      	END IF
      END IF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND .AND. L.NE.LEND-1 ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
     		   IF( TST.LE.( EPS2*ABS( E( M+1 ) ) )*ABS( E( M+1 ) )+
     $             SAFMIN .AND. M.EQ.L) THEN
                  GO TO 60
     		   ELSEIF( TST.LE.( EPS2*ABS( E( M-1 ) ) )*ABS( E( M-1 ) )+
     $             SAFMIN .AND. M.EQ.LENDM1 ) THEN
                  GO TO 60
     		   ELSEIF( TST.LE.( EPS2*ABS( E( M-1 ) ) )*ABS( E( M+1 ) )+
     $             SAFMIN ) THEN
                  GO TO 60
     		   END IF
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
*
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, get its eigensystem directly
*
         IF( M.EQ.L+1 ) THEN
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
*
*        Exit if all iteratives have been done
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        If remaining matrix is 3-by-3, get its eigensystem directly
*
         IF( M.EQ.L+2 ) THEN
            IF ( MOD( JTOT, 10 ).EQ.0 ) THEN
               B = E(L)*E(L) * (ONE - MIN(ABS(E(L+1)/E(L)), ONE))
            ELSE
               B = E(L)*E(L)
            END IF
            P = -E(M-1)*E(M-1) + B
            R = E(M-1)*E(M-2)
            S = SIGN(DLAPY2( P, R ), P)
*
            IF(S.EQ.ZERO) THEN
               VA = -ONE
               VB = ZERO
            ELSE
               VA = -P/S
               VB = -R/S
            END IF
*
*           Update E.
*
            TEMP = E(M-1)
            E(M-1) = VA*E(M-1) - VB*E(M-2)
            E(M-2) = -VB*TEMP - VA*E(M-2)
*
*           If eigenvectors are desired, then update Z initially.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( M-2 ) = -VA
               WORK( N+M-4 ) = VB
               CALL DLASR( 'R', 'G', 'B', N, 3, WORK( M-2 ),
     $                     WORK( N+M-4 ), Z(1, M-2), LDZ )
            END IF
*
            I = L + 1
*
*           Update E.
*
            E(I) = -E(I)
            E(I-1) = -E(I-1)
*
*           If eigenvectors are desired, then update Z.
*
            IF( ICOMPZ.GT.0 ) THEN
               CALL DSCAL(N, -ONE, Z(1, I), 1)
            END IF
*
            GO TO 40
         END IF
*
*        Form shift and set initial values.
*
         IF ( MOD( JTOT, 10 ).EQ.0 ) THEN
            B = E(L)*E(L) * (ONE - MIN(ABS(E(L+1)/E(L)), ONE))
         ELSE
            B = E(L)*E(L)
         END IF
         P = -E(M-1)*E(M-1) + B
         R = E(M-1)*E(M-2)
         S = SIGN(DLAPY2( P, R ), P)
*
         IF(S.EQ.ZERO) THEN
            VA = -ONE
            VB = ZERO
         ELSE
            VA = -P/S
            VB = -R/S
         END IF
*
*        Update E.
*
         TEMP = E(M-1)
         E(M-1) = VA*E(M-1) - VB*E(M-2)
         E(M-2) = -VB*TEMP - VA*E(M-2)
         E3 = E(M-3)
         E(M-3) = -VA*E(M-3)
*
*        If eigenvectors are desired, then update Z initially.
*
         IF( ICOMPZ.GT.0 ) THEN
            WORK( M-2 ) = -VA
            WORK( N+M-4 ) = VB
            CALL DLASR( 'R', 'G', 'B', N, 3, WORK( M-2 ),
     $                  WORK( N+M-4 ), Z(1, M-2), LDZ )
         END IF    
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L+3, -1
*
*           Set iterative values.
*
            P = E(I)
            R = VB*E3
            S = SIGN(DLAPY2( P, R ), P)
            E(I) = -S
*
            IF(S.EQ.ZERO) THEN
               VA = -ONE
               VB = ZERO
            ELSE
               VA = -P/S
               VB = -R/S
            END IF
*
*           Update E.
*
            TEMP = E(I-1)
            E(I-1) = VA*E(I-1) - VB*E(I-2)
            E(I-2) = -VB*TEMP - VA*E(I-2)
            E3 = E(I-3)
            E(I-3) = -VA*E(I-3)
*
*           If eigenvectors are desired, then update Z.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I-2 ) = -VA
               WORK( N+I-4 ) = VB
            END IF    
*
  70     CONTINUE
*
         IF( ICOMPZ.GT.0 ) THEN
            CALL DLASR( 'R', 'G', 'B', N, M-L-1, WORK( L+1 ),
     $                  WORK( N+L-1 ), Z( 1, L+1 ), LDZ )
         END IF
*
         I = L + 2
*
*        Set iterative values.
*
         P = E(I)
         R = VB*E3
         S = SIGN(DLAPY2( P, R ), P)
         E(I) = -S
*
         IF(S.EQ.ZERO) THEN
            VA = -ONE
            VB = ZERO
         ELSE
            VA = -P/S
            VB = -R/S
         END IF
*
*        Update E.
*
         TEMP = E(I-1)
         E(I-1) = VA*E(I-1) - VB*E(I-2)
         E(I-2) = -VB*TEMP - VA*E(I-2)
*
*        If eigenvectors are desired, then update Z.
*
         IF( ICOMPZ.GT.0 ) THEN
            WORK( I-2 ) = -VA
            WORK( N+I-4 ) = VB
            CALL DLASR( 'R', 'G', 'B', N, 3, WORK( I-2 ),
     $                  WORK( N+I-4 ), Z(1, I-2), LDZ )
         END IF
*
         I = L + 1
*
*        Update E.
*
         E(I) = -E(I)
         E(I-1) = -E(I-1)
*
*        If eigenvectors are desired, then update Z.
*
         IF( ICOMPZ.GT.0 ) THEN
            CALL DSCAL(N, -ONE, Z(1, I), 1)
         END IF
*
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND .AND. L.NE.LEND+1 ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
     		   IF( TST.LE.( EPS2*ABS( E( M-2 ) ) )*ABS( E( M-2 ) )+
     $             SAFMIN .AND. M.EQ.L) THEN
                  GO TO 110
     		   ELSEIF( TST.LE.( EPS2*ABS( E( M ) ) )*ABS( E( M ) )+
     $             SAFMIN .AND. M.EQ.LENDP1 ) THEN
                  GO TO 110
     		   ELSEIF( TST.LE.( EPS2*ABS( E( M-2 ) ) )*ABS( E( M ) )+
     $             SAFMIN ) THEN
                  GO TO 110
     		   END IF
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
*
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, get its eigensystem directly
*
         IF( M.EQ.L-1 ) THEN
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
*
*        Exit if all iteratives have been done
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        If remaining matrix is 3-by-3, get its eigensystem directly
*
         IF( M.EQ.L-2 ) THEN
            IF ( MOD( JTOT, 10 ).EQ.0 ) THEN
               B = E(L-1)*E(L-1) * (ONE - MIN(ABS(E(L-2)/E(L-1)), ONE))
            ELSE
               B = E(L-1)*E(L-1)
            END IF
            P = -E(M)*E(M) + B
            R = E(M)*E(M+1)
            S = SIGN(DLAPY2( P, R ), P)
*
            IF(S.EQ.ZERO) THEN
               VA = -ONE
               VB = ZERO
            ELSE
               VA = -P/S
               VB = -R/S
            END IF
*
*           Update E.
*
            TEMP = E(M)
            E(M) = VA*E(M) - VB*E(M+1)
            E(M+1) = -VB*TEMP - VA*E(M+1)
*
*           If eigenvectors are desired, then update Z initially.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( M ) = VA
               WORK( N+M-2 ) = VB
               CALL DLASR( 'R', 'G', 'F', N, 3, WORK( M ),
     $                     WORK( N+M-2 ), Z( 1, M ), LDZ )
            END IF
*
            I = L - 1
*
*           Update E.
*
            E(I-1) = -E(I-1)
            E(I) = -E(I)
*
*           If eigenvectors are desired, then update Z.
*
            IF( ICOMPZ.GT.0 ) THEN
               CALL DSCAL(N, -ONE, Z(1, I), 1)
            END IF
*
            GO TO 90
         END IF
*
*        Form shift and set initial values.
*
         IF ( MOD( JTOT, 10 ).EQ.0 ) THEN
            B = E(L-1)*E(L-1) * (ONE - MIN(ABS(E(L-2)/E(L-1)), ONE))
         ELSE
            B = E(L-1)*E(L-1)
         END IF
         P = -E(M)*E(M) + B
         R = E(M)*E(M+1)
         S = SIGN(DLAPY2( P, R ), P)
*
         IF(S.EQ.ZERO) THEN
            VA = -ONE
            VB = ZERO
         ELSE
            VA = -P/S
            VB = -R/S
         END IF
*
*        Update E.
*
         TEMP = E(M)
         E(M) = VA*E(M) - VB*E(M+1)
         E(M+1) = -VB*TEMP - VA*E(M+1)
         E3 = E(M+2)
         E(M+2) = -VA*E(M+2)
*
*        If eigenvectors are desired, then update Z initially.
*
         IF( ICOMPZ.GT.0 ) THEN
            WORK( M ) = VA
            WORK( N+M-2 ) = VB
            CALL DLASR( 'R', 'G', 'F', N, 3, WORK( M ),
     $                  WORK( N+M-2 ), Z( 1, M ), LDZ )
         END IF    
*
*        Inner loop
*
         LM3 = L - 3
         DO 120 I = M + 1, LM3
*
*           Set iterative values.
*
            P = E(I-1)
            R = VB*E3
            S = SIGN(DLAPY2( P, R ), P)
            E(I-1) = -S
*
            IF(S.EQ.ZERO) THEN
               VA = -ONE
               VB = ZERO
            ELSE
               VA = -P/S
               VB = -R/S
            END IF
*
*           Update E.
*
            TEMP = E(I)
            E(I) = VA*E(I) - VB*E(I+1)
            E(I+1) = -VB*TEMP - VA*E(I+1)
            E3 = E(I+2)
            E(I+2) = -VA*E(I+2)
*
*           If eigenvectors are desired, then update Z.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = VA
               WORK( N+I-2 ) = VB
            END IF    
*
  120    CONTINUE
*
         IF( ICOMPZ.GT.0 ) THEN
            CALL DLASR( 'R', 'G', 'F', N, L-M-1, WORK( M+1 ),
     $                  WORK( N+M-1 ), Z( 1, M+1 ), LDZ )
         END IF
*
         I = L - 2
*
*        Set iterative values.
*
         P = E(I-1)
         R = VB*E3
         S = SIGN(DLAPY2( P, R ), P)
         E(I-1) = -S
*
         IF(S.EQ.ZERO) THEN
            VA = -ONE
            VB = ZERO
         ELSE
            VA = -P/S
            VB = -R/S
         END IF
*
*        Update E.
*
         TEMP = E(I)
         E(I) = VA*E(I) - VB*E(I+1)
         E(I+1) = -VB*TEMP - VA*E(I+1)
*
*        If eigenvectors are desired, then update Z.
*
         IF( ICOMPZ.GT.0 ) THEN
            WORK( I ) = VA
            WORK( N+I-2 ) = VB
            CALL DLASR( 'R', 'G', 'F', N, 3, WORK( I ),
     $                  WORK( N+I-2 ), Z( 1, I ), LDZ )
         END IF
*
         I = L - 1
*
*        Update E.
*
         E(I-1) = -E(I-1)
         E(I) = -E(I)
*
*        If eigenvectors are desired, then update Z.
*
         IF( ICOMPZ.GT.0 ) THEN
            CALL DSCAL(N, -ONE, Z(1, I), 1)
         END IF
*
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
*
      END IF
*
*     Undo scaling if necessary
*
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1,
     $                E(LSV), N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1,
     $                E(LSV), N, INFO )
      END IF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      GO TO 190
*
*     Order blocks.
*     Use Selection Sort to minimize swaps of eigenvectors
*
  160 CONTINUE
      II = 1
      DO WHILE(II.LT.(N-1))
         IF(E(II).EQ.ZERO) THEN
            DO K = II+1,N-1,2
               IF(E(K).EQ.ZERO) THEN
                  DO I = II, K-2
                     E(I) = E(I+1)
                     IF( ICOMPZ.GT.0 ) THEN
                        CALL DSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ), 1 )
                     END IF
                  END DO
                  E(K-1) = ZERO
                  II = K+1
                  EXIT
               ELSEIF(MOD(N,2).EQ.1 .AND. K.EQ.(N-1)) THEN
                  DO I = II, K-1
                     E(I) = E(I+1)
                     IF( ICOMPZ.GT.0 ) THEN
                        CALL DSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ), 1 )
                     END IF
                  END DO
                  IF( ICOMPZ.GT.0 ) THEN
                     CALL DSWAP( N, Z( 1, K ), 1, Z( 1, K+1 ), 1 )
                  END IF
                  E(K) = ZERO
                  II = K+1
                  EXIT
               ELSEIF(MOD(N,2).EQ.0 .AND. K.EQ.(N-2)) THEN
                  DO I = II, K-1
                     E(I) = E(I+1)
                     IF( ICOMPZ.GT.0 ) THEN
                        CALL DSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ), 1 )
                     END IF
                  END DO
                  IF( ICOMPZ.GT.0 ) THEN
                     CALL DSWAP( N, Z( 1, K ), 1, Z( 1, K+1 ), 1 )
                  END IF
                  E(K) = ZERO
                  II = K+1
                  EXIT
               END IF
            END DO
            IF (II.LT.(N-1)) THEN
               CYCLE
            END IF
         END IF
         II = II+2
      END DO
*
      DO 180 II = 1, N-1, 2
         I = II
         P = ABS(E(II))
         DO 170 K = II+2, N-1, 2
            IF(ABS(E(K)).GT.P) THEN
               I = K
               P = ABS(E(K))
            END IF
  170    CONTINUE
         IF(I.NE.II) THEN
            CALL DSWAP( 1, E( I ), 1, E( II ), 1 )
            IF( ICOMPZ.GT.0 ) THEN
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, II ), 1 )
               CALL DSWAP( N, Z( 1, I+1 ), 1, Z( 1, II+1 ), 1 )
            END IF
         END IF
         IF(E(II).LT.ZERO) THEN
            E(II) = -E(II)
            IF( ICOMPZ.GT.0 ) THEN
               CALL DSWAP( N, Z( 1, II ), 1, Z( 1, II+1 ), 1 )
            END IF
         END IF
  180 CONTINUE
*
  190 CONTINUE
      RETURN
*
*     End of DKTEQR
*
      END
