*> \brief \b SKTEDC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SKTEDC + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sktedc.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sktedc.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sktedc.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKTEDC( COMPZ, N, E, Z, LDZ, WORK, LWORK, IWORK,
*                          LIWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          COMPZ
*       INTEGER            INFO, LDZ, LIWORK, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       REAL               E( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKTEDC computes all eigenvalues and, optionally, eigenvectors of a
*> skew-symmetric tridiagonal matrix using the divide and conquer method.
*> The eigenvectors of a full real skew-symmetric matrix can also be
*> found if SKYTRD has been used to reduce this matrix to tridiagonal form.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N':  Compute eigenvalues only.
*>          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*>          = 'V':  Compute eigenvectors of original dense skew-symmetric
*>                  matrix also.  On entry, Z contains the orthogonal
*>                  matrix used to reduce the original matrix to
*>                  tridiagonal form.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The dimension of the skew-symmetric tridiagonal matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is REAL array, dimension (N-1)
*>          On entry, the subdiagonal elements of the tridiagonal matrix.
*>          On exit, if INFO = 0, the block form of eigenvalues in absolute
*>          descending order.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is REAL array, dimension (LDZ,N)
*>          On entry, if COMPZ = 'V', then Z contains the orthogonal
*>          matrix used in the reduction to tridiagonal form.
*>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*>          orthonormal eigenvectors of the original skew-symmetric matrix,
*>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*>          of the skew-symmetric tridiagonal matrix.
*>          If  COMPZ = 'N', then Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1.
*>          If eigenvectors are desired, then LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
*>          If COMPZ = 'V' and N > 1 then LWORK must be at least
*>                         ( 2*N**2 + 2*N + floor(N/2) - 4 ).
*>          If COMPZ = 'I' and N > 1 then LWORK must be at least
*>                         ( N**2 + 2*N + floor(N/2) - 4 ).
*>          Note that for COMPZ = 'I' or 'V', then if N/2 is less than or
*>          equal to the minimum divide size, usually 25, then LWORK need
*>          only be max(1, 2*N - 4).
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*> \endverbatim
*>
*> \param[in] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          The dimension of the array IWORK.
*>          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
*>          If COMPZ = 'V' or 'I', and N > 1 then LIWORK must be at least ( 4*N ).
*>          Note that for COMPZ = 'I' or 'V', then if N/2 is less than or
*>          equal to the minimum divide size, usually 25, then LIWORK
*>          need only be 1.
*>
*>          If LIWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal size of the IWORK array,
*>          returns this value as the first entry of the IWORK array, and
*>          no error message related to LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  if INFO = 1, a singular value did not converge in bidiagonal SVD.
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
*> \ingroup ktedc
*
*> \par Contributors:
*  ==================
*>
*> Jeff Rutter, Computer Science Division, University of California
*> at Berkeley, USA \n
*>  Modified by Francoise Tisseur, University of Tennessee
*>
*  =====================================================================
      SUBROUTINE SKTEDC( COMPZ, N, E, Z, LDZ, WORK, LWORK, IWORK,
     $                   LIWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            FINISH, I, ICOMPZ, II, K, LIWMIN,
     $                   LWMIN, M, SMLSIZ, START, STOREZ, STRTRW
      REAL               EPS, ORGNRM, P, TINY
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANKT, SROUNDUP_LWORK
      EXTERNAL           ILAENV, LSAME, SLAMCH, SLANKT,
     $                   SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SLACPY, SLAEDK0, SLASCL, SLASET,
     $                   SKTEQR, SSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MOD, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
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
      ELSE IF( ( LDZ.LT.1 ) .OR.
     $         ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -5
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         SMLSIZ = ILAENV( 9, 'SKTEDC', ' ', 0, 0, 0, 0 )
         IF( N.LE.1 .OR. ICOMPZ.EQ.0 ) THEN
            LIWMIN = 1
            LWMIN = 1
         ELSE IF( N/2.LE.SMLSIZ ) THEN
            LIWMIN = 1
            LWMIN = MAX(1, 2*N - 4)
         ELSE
            IF( ICOMPZ.EQ.1 ) THEN
               LWMIN = 2*N**2 + 2*N + N/2 - 4
            ELSE IF( ICOMPZ.EQ.2 ) THEN
               LWMIN = N**2 + 2*N + N/2 - 4
            END IF
            LIWMIN = 4*N
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT. LQUERY ) THEN
            INFO = -7
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT. LQUERY ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SKTEDC', -INFO )
         RETURN
      ELSE IF (LQUERY) THEN
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
            IF( ICOMPZ.NE.0 ) THEN
               CALL SSWAP( N, Z( 1, 1 ), 1, Z( 1, 2 ), 1 )
            END IF
         END IF
         RETURN
      END IF
*
*     If COMPZ = 'N', use SKTEQR to compute the eigenvalues.
*
      IF( ICOMPZ.EQ.0 ) THEN
         CALL SKTEQR( 'N', N, E, Z, LDZ, WORK, INFO )
         GO TO 70
      END IF
*
*     If N is smaller than the minimum divide size (SMLSIZ+1), then
*     solve the problem with another solver.
*
      IF( N/2.LE.SMLSIZ ) THEN
*
         CALL SKTEQR( COMPZ, N, E, Z, LDZ, WORK, INFO )
*
      ELSE
*
*        If COMPZ = 'V', the Z matrix must be stored elsewhere for later
*        use.
*
         IF( ICOMPZ.EQ.1 ) THEN
            STOREZ = 1 + N*N
         ELSE
            STOREZ = 1
         END IF
*
         IF( ICOMPZ.EQ.2 ) THEN
            CALL SLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
         END IF
*
*        Scale.
*
         ORGNRM = SLANKT( 'M', N, E )
         IF( ORGNRM.EQ.ZERO )
     $      GO TO 70
*
         EPS = SLAMCH( 'Epsilon' )
*
         START = 1
*
*        while ( START < N )
*
   10    CONTINUE
         IF( START.LT.N ) THEN
*
*           Let FINISH be the position of the next subdiagonal entry
*           such that E( FINISH ) <= TINY or FINISH = N if no such
*           subdiagonal exists.  The matrix identified by the elements
*           between START and FINISH constitutes an independent
*           sub-problem.
*
            FINISH = START
   20       CONTINUE
            IF( FINISH.LT.N ) THEN
               IF( FINISH.EQ.1 ) THEN
                  TINY = EPS*ABS( E( 2 ) )
               ELSEIF( FINISH.EQ.N-1 ) THEN
                  TINY = EPS*ABS( E( N-2 ) )
               ELSE
                  TINY = EPS*SQRT( ABS( E( FINISH-1 ) ) )*
     $                   SQRT( ABS( E( FINISH+1 ) ) )
               END IF
               IF( ABS( E( FINISH ) ).GT.TINY ) THEN
                  FINISH = FINISH + 1
                  GO TO 20
               END IF
            END IF
*
*           (Sub) Problem determined.  Compute its size and solve it.
*
            M = FINISH - START + 1
            IF( M.EQ.1 ) THEN
               START = FINISH + 1
               GO TO 10
            END IF
            IF( M/2.GT.SMLSIZ ) THEN
*
*              Scale.
*
               ORGNRM = SLANKT( 'M', M, E( START ) )
               CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1,
     $                      E( START ), M-1, INFO )
*
               IF( ICOMPZ.EQ.1 ) THEN
                  STRTRW = 1
               ELSE
                  STRTRW = START
               END IF
               CALL SLAEDK0( ICOMPZ, N, M, E( START ),
     $                       Z( STRTRW, START ), LDZ, WORK( 1 ), N,
     $                       WORK( STOREZ ), IWORK, INFO )
               IF( INFO.NE.0 ) THEN
                  GO TO 70
               END IF
*
*              Scale back.
*
               CALL SLASCL( 'G', 0, 0, ONE, ORGNRM, M-1, 1,
     $                      E( START ),  M-1, INFO )
*
            ELSE
               IF( ICOMPZ.EQ.1 ) THEN
*
*                 Since QR won't update a Z matrix which is larger than
*                 the length of E plus 1, we must solve the sub-problem
*                 in a workspace and then multiply back into Z.
*
                  CALL SKTEQR( 'I', M, E( START ), WORK,
     $                         M, WORK( M*M+1 ), INFO )
                  CALL SLACPY( 'A', N, M, Z( 1, START ), LDZ,
     $                         WORK( STOREZ ), N )
                  CALL SGEMM( 'N', 'N', N, M, M, ONE,
     $                        WORK( STOREZ ), N, WORK, M, ZERO,
     $                        Z( 1, START ), LDZ )
               ELSE IF( ICOMPZ.EQ.2 ) THEN
                  CALL SKTEQR( 'I', M, E( START ),
     $                         Z( START, START ), LDZ, WORK, INFO )
               ELSE
                  CALL SKTEQR( 'N', M, E( START ),
     $                         Z( START, START ), LDZ, WORK, INFO )
               END IF
               IF( INFO.NE.0 ) THEN
                  GO TO 70
               END IF
            END IF
*
            START = FINISH + 1
            GO TO 10
         END IF
*
*        endwhile
*
*
*        Order blocks.
*        Use Selection Sort to minimize swaps of eigenvectors
*
         II = 1
         DO WHILE(II.LT.(N-1))
            IF(E(II).EQ.ZERO) THEN
               DO K = II+1,N-1,2
                  IF(E(K).EQ.ZERO) THEN
                     DO I = II, K-2
                        E(I) = E(I+1)
                        IF( ICOMPZ.GT.0 ) THEN
                           CALL SSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ),
     $                                 1 )
                        END IF
                     END DO
                     E(K-1) = ZERO
                     II = K+1
                     EXIT
                  ELSEIF(MOD(N,2).EQ.1 .AND. K.EQ.(N-1)) THEN
                     DO I = II, K-1
                        E(I) = E(I+1)
                        IF( ICOMPZ.GT.0 ) THEN
                           CALL SSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ),
     $                                 1 )
                        END IF
                     END DO
                     IF( ICOMPZ.GT.0 ) THEN
                        CALL SSWAP( N, Z( 1, K ), 1, Z( 1, K+1 ), 1 )
                     END IF
                     E(K) = ZERO
                     II = K+1
                     EXIT
                  ELSEIF(MOD(N,2).EQ.0 .AND. K.EQ.(N-2)) THEN
                     DO I = II, K-1
                        E(I) = E(I+1)
                        IF( ICOMPZ.GT.0 ) THEN
                           CALL SSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ),
     $                                 1 )
                        END IF
                     END DO
                     IF( ICOMPZ.GT.0 ) THEN
                        CALL SSWAP( N, Z( 1, K ), 1, Z( 1, K+1 ), 1 )
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
         DO 60 II = 1, N-1, 2
            I = II
            P = ABS(E(II))
            DO 50 K = II+2, N-1, 2
               IF(ABS(E(K)).GT.P) THEN
                  I = K
                  P = ABS(E(K))
               END IF
  50        CONTINUE
            IF(I.NE.II) THEN
               CALL SSWAP( 1, E( I ), 1, E( II ), 1 )
               IF( ICOMPZ.GT.0 ) THEN
                  CALL SSWAP( N, Z( 1, I ), 1, Z( 1, II ), 1 )
                  CALL SSWAP( N, Z( 1, I+1 ), 1, Z( 1, II+1 ), 1 )
               END IF
            END IF
            IF(E(II).LT.ZERO) THEN
               E(II) = -E(II)
               IF( ICOMPZ.GT.0 ) THEN
                  CALL SSWAP( N, Z( 1, II ), 1, Z( 1, II+1 ), 1 )
               END IF
            END IF
  60     CONTINUE
      END IF
*
   70 CONTINUE
      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of SKTEDC
*
      END
