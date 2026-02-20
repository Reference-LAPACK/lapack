*> \brief \b DLAGTFK computes an LU factorization of a matrix T-λI, where T is a general tridiagonal matrix with purely imaginary diagonal, and λ a purely imaginary scalar, using partial pivoting with row interchanges.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DLAGTFK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlagtfk.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlagtfk.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlagtfk.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAGTFK( N, A, LAMBDA, B, C, TOL, D, IN, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, N
*       DOUBLE PRECISION   LAMBDA, TOL
*       ..
*       .. Array Arguments ..
*       INTEGER            IN( * )
*       DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAGTFK factorizes the matrix (T - lambda*I), where T is an n by n
*> tridiagonal matrix purely imaginary diagonal and lambda is a purely
*> imaginary scalar, as
*>
*>    T - lambda*I = PLU,
*>
*> where P is a permutation matrix, L is a unit lower tridiagonal matrix
*> with at most one non-zero sub-diagonal elements per column and U is
*> an upper triangular matrix with at most two non-zero super-diagonal
*> elements per column.
*>
*> The factorization is obtained by Gaussian elimination with partial
*> pivoting and implicit row scaling.
*>
*> The parameter LAMBDA is included in the routine so that DLAGTFK may
*> be used, in conjunction with SLAGTS, to obtain eigenvectors of T by
*> inverse iteration.
*>
*> Purely imaginary values in this subroutine are stored with imaginary
*> part, and the algorithm guarantees that every value in the
*> factorization are either real, or purely imaginary, which can be
*> determined by the permutation matrix stored in IN.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix T.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (N)
*>          On entry, A must contain the diagonal elements of T.
*>
*>          On exit, A is overwritten by the n diagonal elements of the
*>          upper triangular matrix U of the factorization of T.
*> \endverbatim
*>
*> \param[in] LAMBDA
*> \verbatim
*>          LAMBDA is DOUBLE PRECISION
*>          On entry, the scalar lambda.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (N-1)
*>          On entry, B must contain the (n-1) super-diagonal elements of
*>          T.
*>
*>          On exit, B is overwritten by the (n-1) super-diagonal
*>          elements of the matrix U of the factorization of T.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (N-1)
*>          On entry, C must contain the (n-1) sub-diagonal elements of
*>          T.
*>
*>          On exit, C is overwritten by the (n-1) sub-diagonal elements
*>          of the matrix L of the factorization of T.
*> \endverbatim
*>
*> \param[in] TOL
*> \verbatim
*>          TOL is DOUBLE PRECISION
*>          On entry, a relative tolerance used to indicate whether or
*>          not the matrix (T - lambda*I) is nearly singular. TOL should
*>          normally be chose as approximately the largest relative error
*>          in the elements of T. For example, if the elements of T are
*>          correct to about 4 significant figures, then TOL should be
*>          set to about 5*10**(-4). If TOL is supplied as less than eps,
*>          where eps is the relative machine precision, then the value
*>          eps is used in place of TOL.
*> \endverbatim
*>
*> \param[out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N-2)
*>          On exit, D is overwritten by the (n-2) second super-diagonal
*>          elements of the matrix U of the factorization of T.
*> \endverbatim
*>
*> \param[out] IN
*> \verbatim
*>          IN is INTEGER array, dimension (N+1)
*>          On exit, IN contains details of the permutation matrix P.
*>          The value of IN(k) indicates the interchange occurred at the
*>          kth step of the elimination, and the value type of A(k), B(k),
*>          C(k) and D(k), with the following rules:
*>
*>          D(k) is always real, and
*>          if IN(k) = 0, no interchange, A(k) and C(k) is purely imaginary,
*>          B(k) is real,
*>          if IN(k) = 1, interchange, B(k) and C(k) is purely imaginary,
*>          A(k) is real,
*>          if IN(k) = 2, no interchange, B(k) is purely imaginary, A(k) and
*>          C(k) is real,
*>          if IN(k) = 3, interchange, B(k) is purely imaginary, A(k) and
*>          C(k) is real.
*>
*>          The element IN(n+1) returns the smallest positive integer j
*>          such that
*>
*>             abs( u(j,j) ) <= norm( (T - lambda*I)(j) )*TOL,
*>
*>          where norm( A(j) ) denotes the sum of the absolute values of
*>          the jth row of the matrix A. If no such j exists then IN(n+1)
*>          is returned as zero. If IN(n+1) is returned as positive, then a
*>          diagonal element of U is small, indicating that (T - lambda*I)
*>          is singular or nearly singular.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the kth argument had an illegal value
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
*> \ingroup lagtf
*
*  =====================================================================
      SUBROUTINE DLAGTFK( N, A, LAMBDA, B, C, TOL, D, IN, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
      DOUBLE PRECISION   LAMBDA, TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            K, PRMCT
      DOUBLE PRECISION   EPS, MULT, PIV1, PIV2, SCALE1, SCALE2,
     $                   TEMP, TL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DLAGTFK', -INFO )
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
*
      A( 1 ) = A( 1 ) - LAMBDA
      IN( N+1 ) = 0
      IF( N.EQ.1 ) THEN
         IF( A( 1 ).EQ.ZERO )
     $      IN( 1 ) = 1
         RETURN
      END IF
*
      EPS = DLAMCH( 'Epsilon' )
*
      PRMCT = 0
      TL = MAX( TOL, EPS )
      SCALE1 = ABS( A( 1 ) ) + ABS( B( 1 ) )
      DO 10 K = 1, N - 1
         IF ( MOD( PRMCT, 2 ) .EQ. 0 ) THEN
            IN( K ) = 0
         ELSE
            IN( K ) = 2
         END IF
         A( K+1 ) = A( K+1 ) - LAMBDA
         SCALE2 = ABS( C( K ) ) + ABS( A( K+1 ) )
         IF( K.LT.( N-1 ) )
     $      SCALE2 = SCALE2 + ABS( B( K+1 ) )
         IF( A( K ).EQ.ZERO ) THEN
            PIV1 = ZERO
         ELSE
            PIV1 = ABS( A( K ) ) / SCALE1
         END IF
         IF( C( K ).EQ.ZERO ) THEN
            PIV2 = ZERO
            SCALE1 = SCALE2
            IF( K.LT.( N-1 ) )
     $         D( K ) = ZERO
            PRMCT = 0
         ELSE
            PIV2 = ABS( C( K ) ) / SCALE2
            IF( PIV2.LE.PIV1 ) THEN
               SCALE1 = SCALE2
               IF ( MOD( PRMCT, 2 ) .EQ. 0 ) THEN
                  C( K ) = -C( K ) / A( K )
               ELSE
                  C( K ) = C( K ) / A( K )
               END IF
               A( K+1 ) = A( K+1 ) - C( K )*B( K )
               IF( K.LT.( N-1 ) )
     $            D( K ) = ZERO
               PRMCT = 0
            ELSE
               IN( K ) = IN( K ) + 1
               MULT = A( K ) / C( K )
               A( K ) = C( K )
               TEMP = A( K+1 )
               IF ( MOD( PRMCT, 2 ) .EQ. 0 ) THEN
                  A( K+1 ) = B( K ) + MULT*TEMP
               ELSE
                  A( K+1 ) = B( K ) - MULT*TEMP
               END IF
               IF( K.LT.( N-1 ) ) THEN
                  D( K ) = B( K+1 )
                  B( K+1 ) = -MULT*D( K )
               END IF
               B( K ) = TEMP
               C( K ) = MULT
               PRMCT = PRMCT+1
            END IF
         END IF
         IF( ( MAX( PIV1, PIV2 ).LE.TL ) .AND.
     $       ( IN( N+1 ) .EQ.0 ) )
     $      IN( N+1 ) = K
   10 CONTINUE
*
      IF ( MOD( PRMCT, 2 ) .EQ. 0 ) THEN
         IN( N ) = 0
      ELSE
         IN( N ) = 2
      END IF
*
      IF( ( ABS( A( N ) ).LE.SCALE1*TL ) .AND.
     $    ( IN( N+1 ).EQ.0 ) )
     $   IN( N+1 ) = N
*
      RETURN
*
*     End of DLAGTFK
*
      END
