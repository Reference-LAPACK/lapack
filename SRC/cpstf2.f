*> \brief \b CPSTF2 computes the Cholesky factorization with complete pivoting of complex Hermitian positive semidefinite matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CPSTF2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpstf2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpstf2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpstf2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
*
*       .. Scalar Arguments ..
*       REAL               TOL
*       INTEGER            INFO, LDA, N, RANK
*       CHARACTER          UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX            A( LDA, * )
*       REAL               WORK( 2*N )
*       INTEGER            PIV( N )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CPSTF2 computes the Cholesky factorization with complete
*> pivoting of a complex Hermitian positive semidefinite matrix A.
*>
*> The factorization has the form
*>    P**T * A * P = U**H * U ,  if UPLO = 'U',
*>    P**T * A * P = L  * L**H,  if UPLO = 'L',
*> where U is an upper triangular matrix and L is lower triangular, and
*> P is stored as vector PIV.
*>
*> This algorithm does not attempt to check that A is positive
*> semidefinite. This version of the algorithm calls level 2 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is stored.
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
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          n by n upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading n by n lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky
*>          factorization as above.
*> \endverbatim
*>
*> \param[out] PIV
*> \verbatim
*>          PIV is INTEGER array, dimension (N)
*>          PIV is such that the nonzero entries are P( PIV(K), K ) = 1.
*> \endverbatim
*>
*> \param[out] RANK
*> \verbatim
*>          RANK is INTEGER
*>          The rank of A given by the number of steps the algorithm
*>          completed.
*> \endverbatim
*>
*> \param[in] TOL
*> \verbatim
*>          TOL is REAL
*>          User defined tolerance. If TOL < 0, then N*U*MAX( A( K,K ) )
*>          will be used. The algorithm terminates at the (K-1)st step
*>          if the pivot <= TOL.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (2*N)
*>          Work space.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          < 0: If INFO = -K, the K-th argument had an illegal value,
*>          = 0: algorithm completed successfully, and
*>          > 0: the matrix A is either rank deficient with computed rank
*>               as returned in RANK, or is not positive semidefinite. See
*>               Section 7 of LAPACK Working Note #161 for further
*>               information.
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
*> \ingroup pstf2
*
*  =====================================================================
      SUBROUTINE CPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               TOL
      INTEGER            INFO, LDA, N, RANK
      CHARACTER          UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * )
      REAL               WORK( 2*N )
      INTEGER            PIV( N )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX            CTEMP
      REAL               AJJ, SSTOP, STEMP
      INTEGER            I, ITEMP, J, PVT
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      LOGICAL            LSAME, SISNAN
      EXTERNAL           SLAMCH, LSAME, SISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV, CLACGV, CSSCAL, CSWAP,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPSTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Initialize PIV
*
      DO 100 I = 1, N
         PIV( I ) = I
  100 CONTINUE
*
*     Compute stopping value
*
      DO 110 I = 1, N
         WORK( I ) = REAL( A( I, I ) )
  110 CONTINUE
      PVT = MAXLOC( WORK( 1:N ), 1 )
      AJJ = REAL ( A( PVT, PVT ) )
      IF( AJJ.LE.ZERO.OR.SISNAN( AJJ ) ) THEN
         RANK = 0
         INFO = 1
         GO TO 200
      END IF
*
*     Compute stopping value if not supplied
*
      IF( TOL.LT.ZERO ) THEN
         SSTOP = REAL( N ) * SLAMCH( 'Epsilon' ) * AJJ
      ELSE
         SSTOP = TOL
      END IF
*
*     Set first half of WORK to zero, holds dot products
*
      DO 120 I = 1, N
         WORK( I ) = 0
  120 CONTINUE
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization P**T * A * P = U**H * U
*
         DO 150 J = 1, N
*
*        Find pivot, test for exit, else swap rows and columns
*        Update dot products, compute possible pivots which are
*        stored in the second half of WORK
*
            DO 130 I = J, N
*
               IF( J.GT.1 ) THEN
                  WORK( I ) = WORK( I ) +
     $                        REAL( CONJG( A( J-1, I ) )*
     $                              A( J-1, I ) )
               END IF
               WORK( N+I ) = REAL( A( I, I ) ) - WORK( I )
*
  130       CONTINUE
*
            IF( J.GT.1 ) THEN
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               IF( AJJ.LE.SSTOP.OR.SISNAN( AJJ ) ) THEN
                  A( J, J ) = AJJ
                  GO TO 190
               END IF
            END IF
*
            IF( J.NE.PVT ) THEN
*
*              Pivot OK, so can now swap pivot rows and columns
*
               A( PVT, PVT ) = A( J, J )
               CALL CSWAP( J-1, A( 1, J ), 1, A( 1, PVT ), 1 )
               IF( PVT.LT.N )
     $            CALL CSWAP( N-PVT, A( J, PVT+1 ), LDA,
     $                        A( PVT, PVT+1 ), LDA )
               DO 140 I = J + 1, PVT - 1
                  CTEMP = CONJG( A( J, I ) )
                  A( J, I ) = CONJG( A( I, PVT ) )
                  A( I, PVT ) = CTEMP
  140          CONTINUE
               A( J, PVT ) = CONJG( A( J, PVT ) )
*
*              Swap dot products and PIV
*
               STEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = STEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            END IF
*
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of row J
*
            IF( J.LT.N ) THEN
               CALL CLACGV( J-1, A( 1, J ), 1 )
               CALL CGEMV( 'Trans', J-1, N-J, -CONE, A( 1, J+1 ),
     $                     LDA,
     $                     A( 1, J ), 1, CONE, A( J, J+1 ), LDA )
               CALL CLACGV( J-1, A( 1, J ), 1 )
               CALL CSSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
*
  150    CONTINUE
*
      ELSE
*
*        Compute the Cholesky factorization P**T * A * P = L * L**H
*
         DO 180 J = 1, N
*
*        Find pivot, test for exit, else swap rows and columns
*        Update dot products, compute possible pivots which are
*        stored in the second half of WORK
*
            DO 160 I = J, N
*
               IF( J.GT.1 ) THEN
                  WORK( I ) = WORK( I ) +
     $                        REAL( CONJG( A( I, J-1 ) )*
     $                              A( I, J-1 ) )
               END IF
               WORK( N+I ) = REAL( A( I, I ) ) - WORK( I )
*
  160       CONTINUE
*
            IF( J.GT.1 ) THEN
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               IF( AJJ.LE.SSTOP.OR.SISNAN( AJJ ) ) THEN
                  A( J, J ) = AJJ
                  GO TO 190
               END IF
            END IF
*
            IF( J.NE.PVT ) THEN
*
*              Pivot OK, so can now swap pivot rows and columns
*
               A( PVT, PVT ) = A( J, J )
               CALL CSWAP( J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA )
               IF( PVT.LT.N )
     $            CALL CSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1,
     $                        PVT ),
     $                        1 )
               DO 170 I = J + 1, PVT - 1
                  CTEMP = CONJG( A( I, J ) )
                  A( I, J ) = CONJG( A( PVT, I ) )
                  A( PVT, I ) = CTEMP
  170          CONTINUE
               A( PVT, J ) = CONJG( A( PVT, J ) )
*
*              Swap dot products and PIV
*
               STEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = STEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            END IF
*
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of column J
*
            IF( J.LT.N ) THEN
               CALL CLACGV( J-1, A( J, 1 ), LDA )
               CALL CGEMV( 'No Trans', N-J, J-1, -CONE, A( J+1, 1 ),
     $                     LDA, A( J, 1 ), LDA, CONE, A( J+1, J ), 1 )
               CALL CLACGV( J-1, A( J, 1 ), LDA )
               CALL CSSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
*
  180    CONTINUE
*
      END IF
*
*     Ran to completion, A has full rank
*
      RANK = N
*
      GO TO 200
  190 CONTINUE
*
*     Rank is number of steps completed.  Set INFO = 1 to signal
*     that the factorization cannot be used to solve a system.
*
      RANK = J - 1
      INFO = 1
*
  200 CONTINUE
      RETURN
*
*     End of CPSTF2
*
      END
