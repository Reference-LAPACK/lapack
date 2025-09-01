*> \brief \b CGRQTS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T,
*                          BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDB, LWORK, M, P, N
*       ..
*       .. Array Arguments ..
*       REAL               RESULT( 4 ), RWORK( * )
*       COMPLEX            A( LDA, * ), AF( LDA, * ), R( LDA, * ),
*      $                   Q( LDA, * ), B( LDB, * ), BF( LDB, * ),
*      $                   T( LDB, * ),  Z( LDB, * ), BWK( LDB, * ),
*      $                   TAUA( * ), TAUB( * ), WORK( LWORK )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CGRQTS tests CGGRQF, which computes the GRQ factorization of an
*> M-by-N matrix A and a P-by-N matrix B: A = R*Q and B = Z*T*Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          The M-by-N matrix A.
*> \endverbatim
*>
*> \param[out] AF
*> \verbatim
*>          AF is COMPLEX array, dimension (LDA,N)
*>          Details of the GRQ factorization of A and B, as returned
*>          by CGGRQF, see CGGRQF for further details.
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is COMPLEX array, dimension (LDA,N)
*>          The N-by-N unitary matrix Q.
*> \endverbatim
*>
*> \param[out] R
*> \verbatim
*>          R is COMPLEX array, dimension (LDA,MAX(M,N))
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the arrays A, AF, R and Q.
*>          LDA >= max(M,N).
*> \endverbatim
*>
*> \param[out] TAUA
*> \verbatim
*>          TAUA is COMPLEX array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors, as returned
*>          by SGGQRC.
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,N)
*>          On entry, the P-by-N matrix A.
*> \endverbatim
*>
*> \param[out] BF
*> \verbatim
*>          BF is COMPLEX array, dimension (LDB,N)
*>          Details of the GQR factorization of A and B, as returned
*>          by CGGRQF, see CGGRQF for further details.
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is REAL array, dimension (LDB,P)
*>          The P-by-P unitary matrix Z.
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX array, dimension (LDB,max(P,N))
*> \endverbatim
*>
*> \param[out] BWK
*> \verbatim
*>          BWK is COMPLEX array, dimension (LDB,N)
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the arrays B, BF, Z and T.
*>          LDB >= max(P,N).
*> \endverbatim
*>
*> \param[out] TAUB
*> \verbatim
*>          TAUB is COMPLEX array, dimension (min(P,N))
*>          The scalar factors of the elementary reflectors, as returned
*>          by SGGRQF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK, LWORK >= max(M,P,N)**2.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (M)
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is REAL array, dimension (4)
*>          The test ratios:
*>            RESULT(1) = norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP)
*>            RESULT(2) = norm( T*Q - Z'*B ) / (MAX(P,N)*norm(B)*ULP)
*>            RESULT(3) = norm( I - Q'*Q ) / ( N*ULP )
*>            RESULT(4) = norm( I - Z'*Z ) / ( P*ULP )
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
*> \ingroup complex_eig
*
*  =====================================================================
      SUBROUTINE CGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T,
     $                   BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, P, N
*     ..
*     .. Array Arguments ..
      REAL               RESULT( 4 ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), R( LDA, * ),
     $                   Q( LDA, * ), B( LDB, * ), BF( LDB, * ),
     $                   T( LDB, * ),  Z( LDB, * ), BWK( LDB, * ),
     $                   TAUA( * ), TAUB( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            CROGUE
      PARAMETER          ( CROGUE = ( -1.0E+10, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
      REAL               ANORM, BNORM, ULP, UNFL, RESID
*     ..
*     .. External Functions ..
      REAL               SLAMCH, CLANGE, CLANHE
      EXTERNAL           SLAMCH, CLANGE, CLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CGGRQF, CLACPY, CLASET, CUNGQR,
     $                   CUNGRQ, CHERK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      ULP = SLAMCH( 'Precision' )
      UNFL = SLAMCH( 'Safe minimum' )
*
*     Copy the matrix A to the array AF.
*
      CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL CLACPY( 'Full', P, N, B, LDB, BF, LDB )
*
      ANORM = MAX( CLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( CLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL CGGRQF( M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK,
     $             LWORK, INFO )
*
*     Generate the N-by-N matrix Q
*
      CALL CLASET( 'Full', N, N, CROGUE, CROGUE, Q, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 .AND. M.LT.N )
     $      CALL CLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )
         IF( M.GT.1 )
     $      CALL CLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA,
     $                   Q( N-M+2, N-M+1 ), LDA )
      ELSE
         IF( N.GT.1 )
     $      CALL CLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA,
     $                   Q( 2, 1 ), LDA )
      END IF
      CALL CUNGRQ( N, N, MIN( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO )
*
*     Generate the P-by-P matrix Z
*
      CALL CLASET( 'Full', P, P, CROGUE, CROGUE, Z, LDB )
      IF( P.GT.1 )
     $   CALL CLACPY( 'Lower', P-1, N, BF( 2,1 ), LDB, Z( 2,1 ), LDB )
      CALL CUNGQR( P, P, MIN( P,N ), Z, LDB, TAUB, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL CLASET( 'Full', M, N, CZERO, CZERO, R, LDA )
      IF( M.LE.N )THEN
         CALL CLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ),
     $                LDA )
      ELSE
         CALL CLACPY( 'Full', M-N, N, AF, LDA, R, LDA )
         CALL CLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ),
     $                LDA )
      END IF
*
*     Copy T
*
      CALL CLASET( 'Full', P, N, CZERO, CZERO, T, LDB )
      CALL CLACPY( 'Upper', P, N, BF, LDB, T, LDB )
*
*     Compute R - A*Q'
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', M, N, N, -CONE,
     $            A, LDA, Q, LDA, CONE, R, LDA )
*
*     Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .
*
      RESID = CLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL(MAX(1,M,N) ) ) / ANORM ) / ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute T*Q - Z'*B
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', P, N, P, CONE,
     $           Z, LDB, B, LDB, CZERO, BWK, LDB )
      CALL CGEMM( 'No transpose', 'No transpose', P, N, N, CONE, T, LDB,
     $            Q, LDA, -CONE, BWK, LDB )
*
*     Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
*
      RESID = CLANGE( '1', P, N, BWK, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / REAL( MAX( 1,P,M ) ) )/BNORM ) / ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL CLASET( 'Full', N, N, CZERO, CONE, R, LDA )
      CALL CHERK( 'Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R,
     $            LDA )
*
*     Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = CLANHE( '1', 'Upper', N, R, LDA, RWORK )
      RESULT( 3 ) = ( RESID / REAL( MAX( 1,N ) ) ) / ULP
*
*     Compute I - Z'*Z
*
      CALL CLASET( 'Full', P, P, CZERO, CONE, T, LDB )
      CALL CHERK( 'Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB,
     $            ONE, T, LDB )
*
*     Compute norm( I - Z'*Z ) / ( P*ULP ) .
*
      RESID = CLANHE( '1', 'Upper', P, T, LDB, RWORK )
      RESULT( 4 ) = ( RESID / REAL( MAX( 1,P ) ) ) / ULP
*
      RETURN
*
*     End of CGRQTS
*
      END
