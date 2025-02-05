*> \brief \b ZLAED8 used by ZSTEDC. Merges eigenvalues and deflates secular equation. Used when the original matrix is dense.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZLAED8 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaed8.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaed8.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaed8.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMBDA,
*                          Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR,
*                          GIVCOL, GIVNUM, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ
*       DOUBLE PRECISION   RHO
*       ..
*       .. Array Arguments ..
*       INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
*      $                   INDXQ( * ), PERM( * )
*       DOUBLE PRECISION   D( * ), DLAMBDA( * ), GIVNUM( 2, * ), W( * ),
*      $                   Z( * )
*       COMPLEX*16         Q( LDQ, * ), Q2( LDQ2, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLAED8 merges the two sets of eigenvalues together into a single
*> sorted set.  Then it tries to deflate the size of the problem.
*> There are two ways in which deflation can occur:  when two or more
*> eigenvalues are close together or if there is a tiny element in the
*> Z vector.  For each such occurrence the order of the related secular
*> equation problem is reduced by one.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[out] K
*> \verbatim
*>          K is INTEGER
*>         Contains the number of non-deflated eigenvalues.
*>         This is the order of the related secular equation.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in] QSIZ
*> \verbatim
*>          QSIZ is INTEGER
*>         The dimension of the unitary matrix used to reduce
*>         the dense or band matrix to tridiagonal form.
*>         QSIZ >= N if ICOMPQ = 1.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDQ,N)
*>         On entry, Q contains the eigenvectors of the partially solved
*>         system which has been previously updated in matrix
*>         multiplies with other partially solved eigensystems.
*>         On exit, Q contains the trailing (N-K) updated eigenvectors
*>         (those which were deflated) in its last N-K columns.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>         The leading dimension of the array Q.  LDQ >= max( 1, N ).
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>         On entry, D contains the eigenvalues of the two submatrices to
*>         be combined.  On exit, D contains the trailing (N-K) updated
*>         eigenvalues (those which were deflated) sorted into increasing
*>         order.
*> \endverbatim
*>
*> \param[in,out] RHO
*> \verbatim
*>          RHO is DOUBLE PRECISION
*>         Contains the off diagonal element associated with the rank-1
*>         cut which originally split the two submatrices which are now
*>         being recombined. RHO is modified during the computation to
*>         the value required by DLAED3.
*> \endverbatim
*>
*> \param[in] CUTPNT
*> \verbatim
*>          CUTPNT is INTEGER
*>         Contains the location of the last eigenvalue in the leading
*>         sub-matrix.  MIN(1,N) <= CUTPNT <= N.
*> \endverbatim
*>
*> \param[in] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (N)
*>         On input this vector contains the updating vector (the last
*>         row of the first sub-eigenvector matrix and the first row of
*>         the second sub-eigenvector matrix).  The contents of Z are
*>         destroyed during the updating process.
*> \endverbatim
*>
*> \param[out] DLAMBDA
*> \verbatim
*>          DLAMBDA is DOUBLE PRECISION array, dimension (N)
*>         Contains a copy of the first K eigenvalues which will be used
*>         by DLAED3 to form the secular equation.
*> \endverbatim
*>
*> \param[out] Q2
*> \verbatim
*>          Q2 is COMPLEX*16 array, dimension (LDQ2,N)
*>         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,
*>         Contains a copy of the first K eigenvectors which will be used
*>         by DLAED7 in a matrix multiply (DGEMM) to update the new
*>         eigenvectors.
*> \endverbatim
*>
*> \param[in] LDQ2
*> \verbatim
*>          LDQ2 is INTEGER
*>         The leading dimension of the array Q2.  LDQ2 >= max( 1, N ).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is DOUBLE PRECISION array, dimension (N)
*>         This will hold the first k values of the final
*>         deflation-altered z-vector and will be passed to DLAED3.
*> \endverbatim
*>
*> \param[out] INDXP
*> \verbatim
*>          INDXP is INTEGER array, dimension (N)
*>         This will contain the permutation used to place deflated
*>         values of D at the end of the array. On output INDXP(1:K)
*>         points to the nondeflated D-values and INDXP(K+1:N)
*>         points to the deflated eigenvalues.
*> \endverbatim
*>
*> \param[out] INDX
*> \verbatim
*>          INDX is INTEGER array, dimension (N)
*>         This will contain the permutation used to sort the contents of
*>         D into ascending order.
*> \endverbatim
*>
*> \param[in] INDXQ
*> \verbatim
*>          INDXQ is INTEGER array, dimension (N)
*>         This contains the permutation which separately sorts the two
*>         sub-problems in D into ascending order.  Note that elements in
*>         the second half of this permutation must first have CUTPNT
*>         added to their values in order to be accurate.
*> \endverbatim
*>
*> \param[out] PERM
*> \verbatim
*>          PERM is INTEGER array, dimension (N)
*>         Contains the permutations (from deflation and sorting) to be
*>         applied to each eigenblock.
*> \endverbatim
*>
*> \param[out] GIVPTR
*> \verbatim
*>          GIVPTR is INTEGER
*>         Contains the number of Givens rotations which took place in
*>         this subproblem.
*> \endverbatim
*>
*> \param[out] GIVCOL
*> \verbatim
*>          GIVCOL is INTEGER array, dimension (2, N)
*>         Each pair of numbers indicates a pair of columns to take place
*>         in a Givens rotation.
*> \endverbatim
*>
*> \param[out] GIVNUM
*> \verbatim
*>          GIVNUM is DOUBLE PRECISION array, dimension (2, N)
*>         Each number indicates the S value to be used in the
*>         corresponding Givens rotation.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
*> \ingroup laed8
*
*  =====================================================================
      SUBROUTINE ZLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z,
     $                   DLAMBDA,
     $                   Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR,
     $                   GIVCOL, GIVNUM, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
     $                   INDXQ( * ), PERM( * )
      DOUBLE PRECISION   D( * ), DLAMBDA( * ), GIVNUM( 2, * ), W( * ),
     $                   Z( * )
      COMPLEX*16         Q( LDQ, * ), Q2( LDQ2, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0,
     $                   TWO = 2.0D0, EIGHT = 8.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
      DOUBLE PRECISION   C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           IDAMAX, DLAMCH, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAMRG, DSCAL, XERBLA, ZCOPY,
     $                   ZDROT,
     $                   ZLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( QSIZ.LT.N ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
         INFO = -8
      ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLAED8', -INFO )
         RETURN
      END IF
*
*     Need to initialize GIVPTR to O here in case of quick exit
*     to prevent an unspecified code behavior (usually sigfault)
*     when IWORK array on entry to *stedc is not zeroed
*     (or at least some IWORK entries which used in *laed7 for GIVPTR).
*
      GIVPTR = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1
*
      T = ONE / SQRT( TWO )
      DO 10 J = 1, N
         INDX( J ) = J
   10 CONTINUE
      CALL DSCAL( N, T, Z, 1 )
      RHO = ABS( TWO*RHO )
*
*     Sort the eigenvalues into increasing order
*
      DO 20 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
   20 CONTINUE
      DO 30 I = 1, N
         DLAMBDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
   30 CONTINUE
      I = 1
      J = CUTPNT + 1
      CALL DLAMRG( N1, N2, DLAMBDA, 1, 1, INDX )
      DO 40 I = 1, N
         D( I ) = DLAMBDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
   40 CONTINUE
*
*     Calculate the allowable deflation tolerance
*
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS = DLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*ABS( D( JMAX ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     -- except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         DO 50 J = 1, N
            PERM( J ) = INDXQ( INDX( J ) )
            CALL ZCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
   50    CONTINUE
         CALL ZLACPY( 'A', QSIZ, N, Q2( 1, 1 ), LDQ2, Q( 1, 1 ),
     $                LDQ )
         RETURN
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
      K = 0
      K2 = N + 1
      DO 60 J = 1, N
         IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            INDXP( K2 ) = J
            IF( J.EQ.N )
     $         GO TO 100
         ELSE
            JLAM = J
            GO TO 70
         END IF
   60 CONTINUE
   70 CONTINUE
      J = J + 1
      IF( J.GT.N )
     $   GO TO 90
      IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         INDXP( K2 ) = J
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( JLAM )
         C = Z( J )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = DLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( J ) = TAU
            Z( JLAM ) = ZERO
*
*           Record the appropriate Givens rotation
*
            GIVPTR = GIVPTR + 1
            GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
            GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
            GIVNUM( 1, GIVPTR ) = C
            GIVNUM( 2, GIVPTR ) = S
            CALL ZDROT( QSIZ, Q( 1, INDXQ( INDX( JLAM ) ) ), 1,
     $                  Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
            T = D( JLAM )*C*C + D( J )*S*S
            D( J ) = D( JLAM )*S*S + D( J )*C*C
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
   80       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 80
               ELSE
                  INDXP( K2+I-1 ) = JLAM
               END IF
            ELSE
               INDXP( K2+I-1 ) = JLAM
            END IF
            JLAM = J
         ELSE
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMBDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         END IF
      END IF
      GO TO 70
   90 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      W( K ) = Z( JLAM )
      DLAMBDA( K ) = D( JLAM )
      INDXP( K ) = JLAM
*
  100 CONTINUE
*
*     Sort the eigenvalues and corresponding eigenvectors into DLAMBDA
*     and Q2 respectively.  The eigenvalues/vectors which were not
*     deflated go into the first K slots of DLAMBDA and Q2 respectively,
*     while those which were deflated go into the last N - K slots.
*
      DO 110 J = 1, N
         JP = INDXP( J )
         DLAMBDA( J ) = D( JP )
         PERM( J ) = INDXQ( INDX( JP ) )
         CALL ZCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  110 CONTINUE
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      IF( K.LT.N ) THEN
         CALL DCOPY( N-K, DLAMBDA( K+1 ), 1, D( K+1 ), 1 )
         CALL ZLACPY( 'A', QSIZ, N-K, Q2( 1, K+1 ), LDQ2, Q( 1,
     $                K+1 ),
     $                LDQ )
      END IF
*
      RETURN
*
*     End of ZLAED8
*
      END
