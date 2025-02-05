*> \brief \b SLAED7 used by SSTEDC. Computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is dense.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SLAED7 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed7.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed7.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed7.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,
*                          LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR,
*                          PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N,
*      $                   QSIZ, TLVLS
*       REAL               RHO
*       ..
*       .. Array Arguments ..
*       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),
*      $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
*       REAL               D( * ), GIVNUM( 2, * ), Q( LDQ, * ),
*      $                   QSTORE( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SLAED7 computes the updated eigensystem of a diagonal
*> matrix after modification by a rank-one symmetric matrix. This
*> routine is used only for the eigenproblem which requires all
*> eigenvalues and optionally eigenvectors of a dense symmetric matrix
*> that has been reduced to tridiagonal form.  SLAED1 handles
*> the case in which all eigenvalues and eigenvectors of a symmetric
*> tridiagonal matrix are desired.
*>
*>   T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
*>
*>    where Z = Q**Tu, u is a vector of length N with ones in the
*>    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
*>
*>    The eigenvectors of the original matrix are stored in Q, and the
*>    eigenvalues are in D.  The algorithm consists of three stages:
*>
*>       The first stage consists of deflating the size of the problem
*>       when there are multiple eigenvalues or if there is a zero in
*>       the Z vector.  For each such occurrence the dimension of the
*>       secular equation problem is reduced by one.  This stage is
*>       performed by the routine SLAED8.
*>
*>       The second stage consists of calculating the updated
*>       eigenvalues. This is done by finding the roots of the secular
*>       equation via the routine SLAED4 (as called by SLAED9).
*>       This routine also calculates the eigenvectors of the current
*>       problem.
*>
*>       The final stage consists of computing the updated eigenvectors
*>       directly using the updated eigenvalues.  The eigenvectors for
*>       the current problem are multiplied with the eigenvectors from
*>       the overall problem.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ICOMPQ
*> \verbatim
*>          ICOMPQ is INTEGER
*>          = 0:  Compute eigenvalues only.
*>          = 1:  Compute eigenvectors of original dense symmetric matrix
*>                also.  On entry, Q contains the orthogonal matrix used
*>                to reduce the original matrix to tridiagonal form.
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
*>         The dimension of the orthogonal matrix used to reduce
*>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
*> \endverbatim
*>
*> \param[in] TLVLS
*> \verbatim
*>          TLVLS is INTEGER
*>         The total number of merging levels in the overall divide and
*>         conquer tree.
*> \endverbatim
*>
*> \param[in] CURLVL
*> \verbatim
*>          CURLVL is INTEGER
*>         The current level in the overall merge routine,
*>         0 <= CURLVL <= TLVLS.
*> \endverbatim
*>
*> \param[in] CURPBM
*> \verbatim
*>          CURPBM is INTEGER
*>         The current problem in the current level in the overall
*>         merge routine (counting from upper left to lower right).
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is REAL array, dimension (N)
*>         On entry, the eigenvalues of the rank-1-perturbed matrix.
*>         On exit, the eigenvalues of the repaired matrix.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is REAL array, dimension (LDQ, N)
*>         On entry, the eigenvectors of the rank-1-perturbed matrix.
*>         On exit, the eigenvectors of the repaired tridiagonal matrix.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>         The leading dimension of the array Q.  LDQ >= max(1,N).
*> \endverbatim
*>
*> \param[out] INDXQ
*> \verbatim
*>          INDXQ is INTEGER array, dimension (N)
*>         The permutation which will reintegrate the subproblem just
*>         solved back into sorted order, i.e., D( INDXQ( I = 1, N ) )
*>         will be in ascending order.
*> \endverbatim
*>
*> \param[in] RHO
*> \verbatim
*>          RHO is REAL
*>         The subdiagonal element used to create the rank-1
*>         modification.
*> \endverbatim
*>
*> \param[in] CUTPNT
*> \verbatim
*>          CUTPNT is INTEGER
*>         Contains the location of the last eigenvalue in the leading
*>         sub-matrix.  min(1,N) <= CUTPNT <= N.
*> \endverbatim
*>
*> \param[in,out] QSTORE
*> \verbatim
*>          QSTORE is REAL array, dimension (N**2+1)
*>         Stores eigenvectors of submatrices encountered during
*>         divide and conquer, packed together. QPTR points to
*>         beginning of the submatrices.
*> \endverbatim
*>
*> \param[in,out] QPTR
*> \verbatim
*>          QPTR is INTEGER array, dimension (N+2)
*>         List of indices pointing to beginning of submatrices stored
*>         in QSTORE. The submatrices are numbered starting at the
*>         bottom left of the divide and conquer tree, from left to
*>         right and bottom to top.
*> \endverbatim
*>
*> \param[in] PRMPTR
*> \verbatim
*>          PRMPTR is INTEGER array, dimension (N lg N)
*>         Contains a list of pointers which indicate where in PERM a
*>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
*>         indicates the size of the permutation and also the size of
*>         the full, non-deflated problem.
*> \endverbatim
*>
*> \param[in] PERM
*> \verbatim
*>          PERM is INTEGER array, dimension (N lg N)
*>         Contains the permutations (from deflation and sorting) to be
*>         applied to each eigenblock.
*> \endverbatim
*>
*> \param[in] GIVPTR
*> \verbatim
*>          GIVPTR is INTEGER array, dimension (N lg N)
*>         Contains a list of pointers which indicate where in GIVCOL a
*>         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
*>         indicates the number of Givens rotations.
*> \endverbatim
*>
*> \param[in] GIVCOL
*> \verbatim
*>          GIVCOL is INTEGER array, dimension (2, N lg N)
*>         Each pair of numbers indicates a pair of columns to take place
*>         in a Givens rotation.
*> \endverbatim
*>
*> \param[in] GIVNUM
*> \verbatim
*>          GIVNUM is REAL array, dimension (2, N lg N)
*>         Each number indicates the S value to be used in the
*>         corresponding Givens rotation.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (3*N+2*QSIZ*N)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (4*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  if INFO = 1, an eigenvalue did not converge
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
*> \ingroup laed7
*
*> \par Contributors:
*  ==================
*>
*> Jeff Rutter, Computer Science Division, University of California
*> at Berkeley, USA
*
*  =====================================================================
      SUBROUTINE SLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D,
     $                   Q,
     $                   LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR,
     $                   PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N,
     $                   QSIZ, TLVLS
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),
     $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
      REAL               D( * ), GIVNUM( 2, * ), Q( LDQ, * ),
     $                   QSTORE( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            COLTYP, CURR, I, IDLMDA, INDX, INDXC, INDXP,
     $                   IQ2, IS, IW, IZ, K, LDQ2, N1, N2, PTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SLAED8, SLAED9, SLAEDA, SLAMRG,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( MIN( 1, N ).GT.CUTPNT .OR. N.LT.CUTPNT ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED7', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     The following values are for bookkeeping purposes only.  They are
*     integer pointers which indicate the portion of the workspace
*     used by a particular array in SLAED8 and SLAED9.
*
      IF( ICOMPQ.EQ.1 ) THEN
         LDQ2 = QSIZ
      ELSE
         LDQ2 = N
      END IF
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ2 = IW + N
      IS = IQ2 + N*LDQ2
*
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
*
*     Form the z-vector which consists of the last row of Q_1 and the
*     first row of Q_2.
*
      PTR = 1 + 2**TLVLS
      DO 10 I = 1, CURLVL - 1
         PTR = PTR + 2**( TLVLS-I )
   10 CONTINUE
      CURR = PTR + CURPBM
      CALL SLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,
     $             GIVCOL, GIVNUM, QSTORE, QPTR, WORK( IZ ),
     $             WORK( IZ+N ), INFO )
*
*     When solving the final problem, we no longer need the stored data,
*     so we will overwrite the data from this level onto the previously
*     used storage space.
*
      IF( CURLVL.EQ.TLVLS ) THEN
         QPTR( CURR ) = 1
         PRMPTR( CURR ) = 1
         GIVPTR( CURR ) = 1
      END IF
*
*     Sort and Deflate eigenvalues.
*
      CALL SLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, CUTPNT,
     $             WORK( IZ ), WORK( IDLMDA ), WORK( IQ2 ), LDQ2,
     $             WORK( IW ), PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ),
     $             GIVCOL( 1, GIVPTR( CURR ) ),
     $             GIVNUM( 1, GIVPTR( CURR ) ), IWORK( INDXP ),
     $             IWORK( INDX ), INFO )
      PRMPTR( CURR+1 ) = PRMPTR( CURR ) + N
      GIVPTR( CURR+1 ) = GIVPTR( CURR+1 ) + GIVPTR( CURR )
*
*     Solve Secular Equation.
*
      IF( K.NE.0 ) THEN
         CALL SLAED9( K, 1, K, N, D, WORK( IS ), K, RHO,
     $                WORK( IDLMDA ),
     $                WORK( IW ), QSTORE( QPTR( CURR ) ), K, INFO )
         IF( INFO.NE.0 )
     $      GO TO 30
         IF( ICOMPQ.EQ.1 ) THEN
            CALL SGEMM( 'N', 'N', QSIZ, K, K, ONE, WORK( IQ2 ), LDQ2,
     $                  QSTORE( QPTR( CURR ) ), K, ZERO, Q, LDQ )
         END IF
         QPTR( CURR+1 ) = QPTR( CURR ) + K**2
*
*     Prepare the INDXQ sorting permutation.
*
         N1 = K
         N2 = N - K
         CALL SLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         QPTR( CURR+1 ) = QPTR( CURR )
         DO 20 I = 1, N
            INDXQ( I ) = I
   20    CONTINUE
      END IF
*
   30 CONTINUE
      RETURN
*
*     End of SLAED7
*
      END
