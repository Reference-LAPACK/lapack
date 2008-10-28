      SUBROUTINE CPTTS2( IUPLO, N, NRHS, D, E, B, LDB )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IUPLO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL               D( * )
      COMPLEX            B( LDB, * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  CPTTS2 solves a tridiagonal system of the form
*     A * X = B
*  using the factorization A = U'*D*U or A = L*D*L' computed by CPTTRF.
*  D is a diagonal matrix specified in the vector D, U (or L) is a unit
*  bidiagonal matrix whose superdiagonal (subdiagonal) is specified in
*  the vector E, and X and B are N by NRHS matrices.
*
*  Arguments
*  =========
*
*  IUPLO   (input) INTEGER
*          Specifies the form of the factorization and whether the
*          vector E is the superdiagonal of the upper bidiagonal factor
*          U or the subdiagonal of the lower bidiagonal factor L.
*          = 1:  A = U'*D*U, E is the superdiagonal of U
*          = 0:  A = L*D*L', E is the subdiagonal of L
*
*  N       (input) INTEGER
*          The order of the tridiagonal matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  D       (input) REAL array, dimension (N)
*          The n diagonal elements of the diagonal matrix D from the
*          factorization A = U'*D*U or A = L*D*L'.
*
*  E       (input) COMPLEX array, dimension (N-1)
*          If IUPLO = 1, the (n-1) superdiagonal elements of the unit
*          bidiagonal factor U from the factorization A = U'*D*U.
*          If IUPLO = 0, the (n-1) subdiagonal elements of the unit
*          bidiagonal factor L from the factorization A = L*D*L'.
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the right hand side vectors B for the system of
*          linear equations.
*          On exit, the solution vectors, X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL           CSSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 ) THEN
         IF( N.EQ.1 )
     $      CALL CSSCAL( NRHS, 1. / D( 1 ), B, LDB )
         RETURN
      END IF
*
      IF( IUPLO.EQ.1 ) THEN
*
*        Solve A * X = B using the factorization A = U'*D*U,
*        overwriting each right hand side vector with its solution.
*
         IF( NRHS.LE.2 ) THEN
            J = 1
    5       CONTINUE
*
*           Solve U' * x = b.
*
            DO 10 I = 2, N
               B( I, J ) = B( I, J ) - B( I-1, J )*CONJG( E( I-1 ) )
   10       CONTINUE
*
*           Solve D * U * x = b.
*
            DO 20 I = 1, N
               B( I, J ) = B( I, J ) / D( I )
   20       CONTINUE
            DO 30 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) - B( I+1, J )*E( I )
   30       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 5
            END IF
         ELSE
            DO 60 J = 1, NRHS
*
*              Solve U' * x = b.
*
               DO 40 I = 2, N
                  B( I, J ) = B( I, J ) - B( I-1, J )*CONJG( E( I-1 ) )
   40          CONTINUE
*
*              Solve D * U * x = b.
*
               B( N, J ) = B( N, J ) / D( N )
               DO 50 I = N - 1, 1, -1
                  B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
   50          CONTINUE
   60       CONTINUE
         END IF
      ELSE
*
*        Solve A * X = B using the factorization A = L*D*L',
*        overwriting each right hand side vector with its solution.
*
         IF( NRHS.LE.2 ) THEN
            J = 1
   65       CONTINUE
*
*           Solve L * x = b.
*
            DO 70 I = 2, N
               B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
   70       CONTINUE
*
*           Solve D * L' * x = b.
*
            DO 80 I = 1, N
               B( I, J ) = B( I, J ) / D( I )
   80       CONTINUE
            DO 90 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) - B( I+1, J )*CONJG( E( I ) )
   90       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 65
            END IF
         ELSE
            DO 120 J = 1, NRHS
*
*              Solve L * x = b.
*
               DO 100 I = 2, N
                  B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
  100          CONTINUE
*
*              Solve D * L' * x = b.
*
               B( N, J ) = B( N, J ) / D( N )
               DO 110 I = N - 1, 1, -1
                  B( I, J ) = B( I, J ) / D( I ) -
     $                        B( I+1, J )*CONJG( E( I ) )
  110          CONTINUE
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of CPTTS2
*
      END
