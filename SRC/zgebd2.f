*> \brief \b ZGEBD2 reduces a general matrix to bidiagonal form using an unblocked algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZGEBD2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgebd2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgebd2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgebd2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   D( * ), E( * )
*       COMPLEX*16         A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEBD2 reduces a complex general m by n matrix A to upper or lower
*> real bidiagonal form B by a unitary transformation: Q**H * A * P = B.
*>
*> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows in the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns in the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the m by n general matrix to be reduced.
*>          On exit,
*>          if m >= n, the diagonal and the first superdiagonal are
*>            overwritten with the upper bidiagonal matrix B; the
*>            elements below the diagonal, with the array TAUQ, represent
*>            the unitary matrix Q as a product of elementary
*>            reflectors, and the elements above the first superdiagonal,
*>            with the array TAUP, represent the unitary matrix P as
*>            a product of elementary reflectors;
*>          if m < n, the diagonal and the first subdiagonal are
*>            overwritten with the lower bidiagonal matrix B; the
*>            elements below the first subdiagonal, with the array TAUQ,
*>            represent the unitary matrix Q as a product of
*>            elementary reflectors, and the elements above the diagonal,
*>            with the array TAUP, represent the unitary matrix P as
*>            a product of elementary reflectors.
*>          See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (min(M,N))
*>          The diagonal elements of the bidiagonal matrix B:
*>          D(i) = A(i,i).
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (min(M,N)-1)
*>          The off-diagonal elements of the bidiagonal matrix B:
*>          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*>          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
*> \endverbatim
*>
*> \param[out] TAUQ
*> \verbatim
*>          TAUQ is COMPLEX*16 array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors which
*>          represent the unitary matrix Q. See Further Details.
*> \endverbatim
*>
*> \param[out] TAUP
*> \verbatim
*>          TAUP is COMPLEX*16 array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors which
*>          represent the unitary matrix P. See Further Details.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (max(M,N))
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value.
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
*> \ingroup gebd2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrices Q and P are represented as products of elementary
*>  reflectors:
*>
*>  If m >= n,
*>
*>     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
*>
*>  Each H(i) and G(i) has the form:
*>
*>     H(i) = I - tauq * v * v**H  and G(i) = I - taup * u * u**H
*>
*>  where tauq and taup are complex scalars, and v and u are complex
*>  vectors; v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in
*>  A(i+1:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in
*>  A(i,i+2:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
*>
*>  If m < n,
*>
*>     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
*>
*>  Each H(i) and G(i) has the form:
*>
*>     H(i) = I - tauq * v * v**H  and G(i) = I - taup * u * u**H
*>
*>  where tauq and taup are complex scalars, v and u are complex vectors;
*>  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
*>  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
*>  tauq is stored in TAUQ(i) and taup in TAUP(i).
*>
*>  The contents of A on exit are illustrated by the following examples:
*>
*>  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*>
*>    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
*>    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
*>    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
*>    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
*>    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
*>    (  v1  v2  v3  v4  v5 )
*>
*>  where d and e denote diagonal and off-diagonal elements of B, vi
*>  denotes an element of the vector defining H(i), and ui an element of
*>  the vector defining G(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLARF1F, ZLARFG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'ZGEBD2', -INFO )
         RETURN
      END IF
*
      IF( M.GE.N ) THEN
*
*        Reduce to upper bidiagonal form
*
         DO 10 I = 1, N
*
*           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
            ALPHA = A( I, I )
            CALL ZLARFG( M-I+1, ALPHA, A( MIN( I+1, M ), I ), 1,
     $                   TAUQ( I ) )
            D( I ) = DBLE( ALPHA )
*
*           Apply H(i)**H to A(i:m,i+1:n) from the left
*
            IF( I.LT.N )
     $         CALL ZLARF1F( 'Left', M-I+1, N-I, A( I, I ), 1,
     $                     DCONJG( TAUQ( I ) ), A( I, I+1 ), LDA, WORK )
            A( I, I ) = D( I )
*
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector G(i) to annihilate
*              A(i,i+2:n)
*
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               ALPHA = A( I, I+1 )
               CALL ZLARFG( N-I, ALPHA, A( I, MIN( I+2, N ) ), LDA,
     $                      TAUP( I ) )
               E( I ) = DBLE( ALPHA )
*
*              Apply G(i) to A(i+1:m,i+1:n) from the right
*
               CALL ZLARF1F( 'Right', M-I, N-I, A( I, I+1 ), LDA,
     $                     TAUP( I ), A( I+1, I+1 ), LDA, WORK )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I+1 ) = E( I )
            ELSE
               TAUP( I ) = ZERO
            END IF
   10    CONTINUE
      ELSE
*
*        Reduce to lower bidiagonal form
*
         DO 20 I = 1, M
*
*           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
*
            CALL ZLACGV( N-I+1, A( I, I ), LDA )
            ALPHA = A( I, I )
            CALL ZLARFG( N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA,
     $                   TAUP( I ) )
            D( I ) = DBLE( ALPHA )
*
*           Apply G(i) to A(i+1:m,i:n) from the right
*
            IF( I.LT.M )
     $         CALL ZLARF1F( 'Right', M-I, N-I+1, A( I, I ), LDA,
     $                     TAUP( I ), A( I+1, I ), LDA, WORK )
            CALL ZLACGV( N-I+1, A( I, I ), LDA )
            A( I, I ) = D( I )
*
            IF( I.LT.M ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:m,i)
*
               ALPHA = A( I+1, I )
               CALL ZLARFG( M-I, ALPHA, A( MIN( I+2, M ), I ), 1,
     $                      TAUQ( I ) )
               E( I ) = DBLE( ALPHA )
*
*              Apply H(i)**H to A(i+1:m,i+1:n) from the left
*
               CALL ZLARF1F( 'Left', M-I, N-I, A( I+1, I ), 1,
     $                     DCONJG( TAUQ( I ) ), A( I+1, I+1 ), LDA,
     $                     WORK )
               A( I+1, I ) = E( I )
            ELSE
               TAUQ( I ) = ZERO
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZGEBD2
*
      END
