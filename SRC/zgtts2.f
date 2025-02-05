*> \brief \b ZGTTS2 solves a system of linear equations with a tridiagonal matrix using the LU factorization computed by sgttrf.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZGTTS2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgtts2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgtts2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgtts2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )
*
*       .. Scalar Arguments ..
*       INTEGER            ITRANS, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGTTS2 solves one of the systems of equations
*>    A * X = B,  A**T * X = B,  or  A**H * X = B,
*> with a tridiagonal matrix A using the LU factorization computed
*> by ZGTTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ITRANS
*> \verbatim
*>          ITRANS is INTEGER
*>          Specifies the form of the system of equations.
*>          = 0:  A * X = B     (No transpose)
*>          = 1:  A**T * X = B  (Transpose)
*>          = 2:  A**H * X = B  (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] DL
*> \verbatim
*>          DL is COMPLEX*16 array, dimension (N-1)
*>          The (n-1) multipliers that define the matrix L from the
*>          LU factorization of A.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is COMPLEX*16 array, dimension (N)
*>          The n diagonal elements of the upper triangular matrix U from
*>          the LU factorization of A.
*> \endverbatim
*>
*> \param[in] DU
*> \verbatim
*>          DU is COMPLEX*16 array, dimension (N-1)
*>          The (n-1) elements of the first super-diagonal of U.
*> \endverbatim
*>
*> \param[in] DU2
*> \verbatim
*>          DU2 is COMPLEX*16 array, dimension (N-2)
*>          The (n-2) elements of the second super-diagonal of U.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices; for 1 <= i <= n, row i of the matrix was
*>          interchanged with row IPIV(i).  IPIV(i) will always be either
*>          i or i+1; IPIV(i) = i indicates a row interchange was not
*>          required.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)
*>          On entry, the matrix of right hand side vectors B.
*>          On exit, B is overwritten by the solution vectors X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
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
*> \ingroup gtts2
*
*  =====================================================================
      SUBROUTINE ZGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B,
     $                   LDB )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            ITRANS, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
      COMPLEX*16         TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( ITRANS.EQ.0 ) THEN
*
*        Solve A*X = B using the LU factorization of A,
*        overwriting each right hand side vector with its solution.
*
         IF( NRHS.LE.1 ) THEN
            J = 1
   10       CONTINUE
*
*           Solve L*x = b.
*
            DO 20 I = 1, N - 1
               IF( IPIV( I ).EQ.I ) THEN
                  B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
               ELSE
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - DL( I )*B( I, J )
               END IF
   20       CONTINUE
*
*           Solve U*x = b.
*
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 )
     $         B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) /
     $                       D( N-1 )
            DO 30 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DU2( I )*
     $                     B( I+2, J ) ) / D( I )
   30       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 10
            END IF
         ELSE
            DO 60 J = 1, NRHS
*
*           Solve L*x = b.
*
               DO 40 I = 1, N - 1
                  IF( IPIV( I ).EQ.I ) THEN
                     B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
                  ELSE
                     TEMP = B( I, J )
                     B( I, J ) = B( I+1, J )
                     B( I+1, J ) = TEMP - DL( I )*B( I, J )
                  END IF
   40          CONTINUE
*
*           Solve U*x = b.
*
               B( N, J ) = B( N, J ) / D( N )
               IF( N.GT.1 )
     $            B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) /
     $                          D( N-1 )
               DO 50 I = N - 2, 1, -1
                  B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DU2( I )*
     $                        B( I+2, J ) ) / D( I )
   50          CONTINUE
   60       CONTINUE
         END IF
      ELSE IF( ITRANS.EQ.1 ) THEN
*
*        Solve A**T * X = B.
*
         IF( NRHS.LE.1 ) THEN
            J = 1
   70       CONTINUE
*
*           Solve U**T * x = b.
*
            B( 1, J ) = B( 1, J ) / D( 1 )
            IF( N.GT.1 )
     $         B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
            DO 80 I = 3, N
               B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )-DU2( I-2 )*
     $                     B( I-2, J ) ) / D( I )
   80       CONTINUE
*
*           Solve L**T * x = b.
*
            DO 90 I = N - 1, 1, -1
               IF( IPIV( I ).EQ.I ) THEN
                  B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
               ELSE
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                  B( I, J ) = TEMP
               END IF
   90       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 70
            END IF
         ELSE
            DO 120 J = 1, NRHS
*
*           Solve U**T * x = b.
*
               B( 1, J ) = B( 1, J ) / D( 1 )
               IF( N.GT.1 )
     $            B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
               DO 100 I = 3, N
                  B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )-
     $                        DU2( I-2 )*B( I-2, J ) ) / D( I )
  100          CONTINUE
*
*           Solve L**T * x = b.
*
               DO 110 I = N - 1, 1, -1
                  IF( IPIV( I ).EQ.I ) THEN
                     B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
                  ELSE
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                     B( I, J ) = TEMP
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
*
*        Solve A**H * X = B.
*
         IF( NRHS.LE.1 ) THEN
            J = 1
  130       CONTINUE
*
*           Solve U**H * x = b.
*
            B( 1, J ) = B( 1, J ) / DCONJG( D( 1 ) )
            IF( N.GT.1 )
     $         B( 2, J ) = ( B( 2, J )-DCONJG( DU( 1 ) )*B( 1, J ) ) /
     $                     DCONJG( D( 2 ) )
            DO 140 I = 3, N
               B( I, J ) = ( B( I, J )-DCONJG( DU( I-1 ) )*B( I-1, J )-
     $                     DCONJG( DU2( I-2 ) )*B( I-2, J ) ) /
     $                     DCONJG( D( I ) )
  140       CONTINUE
*
*           Solve L**H * x = b.
*
            DO 150 I = N - 1, 1, -1
               IF( IPIV( I ).EQ.I ) THEN
                  B( I, J ) = B( I, J ) - DCONJG( DL( I ) )*B( I+1, J )
               ELSE
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - DCONJG( DL( I ) )*TEMP
                  B( I, J ) = TEMP
               END IF
  150       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 130
            END IF
         ELSE
            DO 180 J = 1, NRHS
*
*           Solve U**H * x = b.
*
               B( 1, J ) = B( 1, J ) / DCONJG( D( 1 ) )
               IF( N.GT.1 )
     $            B( 2, J ) = ( B( 2, J )-DCONJG( DU( 1 ) )*B( 1, J ) )
     $                         / DCONJG( D( 2 ) )
               DO 160 I = 3, N
                  B( I, J ) = ( B( I, J )-DCONJG( DU( I-1 ) )*
     $                        B( I-1, J )-DCONJG( DU2( I-2 ) )*
     $                        B( I-2, J ) ) / DCONJG( D( I ) )
  160          CONTINUE
*
*           Solve L**H * x = b.
*
               DO 170 I = N - 1, 1, -1
                  IF( IPIV( I ).EQ.I ) THEN
                     B( I, J ) = B( I, J ) - DCONJG( DL( I ) )*
     $                           B( I+1, J )
                  ELSE
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - DCONJG( DL( I ) )*TEMP
                     B( I, J ) = TEMP
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      END IF
*
*     End of ZGTTS2
*
      END
