*> \brief \b DLAEDK0 used by DKTEDC. Computes all eigenvalues and corresponding eigenvectors of an unreduced skew-symmetric tridiagonal matrix using the divide and conquer method.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DLAEDK0 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaedk0.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaedk0.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaedk0.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAEDK0( ICOMPQ, QSIZ, N, E, Q, LDQ, QSTORE, LDQS,
*                          WORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   E( * ), Q( LDQ, * ), QSTORE( LDQS, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAEDK0 computes all eigenvalues and corresponding eigenvectors of a
*> skew-symmetric tridiagonal matrix using the divide and conquer method.
*> This subroutine firstly perform an odd-even reordering and then use
*> divide and conquer method of bidiagonal SVD.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ICOMPQ
*> \verbatim
*>          ICOMPQ is INTEGER
*>          = 0:  Compute eigenvalues only.
*>          = 1:  Compute eigenvectors of original dense skew-symmetric matrix
*>                also.  On entry, Q contains the orthogonal matrix used
*>                to reduce the original matrix to tridiagonal form.
*>          = 2:  Compute eigenvalues and eigenvectors of tridiagonal
*>                matrix.
*> \endverbatim
*>
*> \param[in] QSIZ
*> \verbatim
*>          QSIZ is INTEGER
*>         The dimension of the orthogonal matrix used to reduce
*>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         The dimension of the skew-symmetric tridiagonal matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>         The off-diagonal elements of the tridiagonal matrix.
*>         On exit, if INFO = 0, the block form of eigenvalues in absolute
*>         descending order.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
*>         On entry, Q must contain an N-by-N orthogonal matrix.
*>         If ICOMPQ = 0    Q is not referenced.
*>         If ICOMPQ = 1    On entry, Q is a subset of the columns of the
*>                          orthogonal matrix used to reduce the full
*>                          matrix to tridiagonal form corresponding to
*>                          the subset of the full matrix which is being
*>                          decomposed at this time.
*>         If ICOMPQ = 2    On entry, Q will be the identity matrix.
*>                          On exit, Q contains the eigenvectors of the
*>                          tridiagonal matrix.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>         The leading dimension of the array Q.  If eigenvectors are
*>         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.
*> \endverbatim
*>
*> \param[out] QSTORE
*> \verbatim
*>          QSTORE is DOUBLE PRECISION array, dimension (LDQS, N)
*>         Referenced only when ICOMPQ = 1.  Used to store parts of
*>         the eigenvector matrix when the updating matrix multiplies
*>         take place.
*> \endverbatim
*>
*> \param[in] LDQS
*> \verbatim
*>          LDQS is INTEGER
*>         The leading dimension of the array QSTORE.  If ICOMPQ = 1,
*>         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array,
*>         If ICOMPQ = 0, WORK is not referenced.
*>         If ICOMPQ = 1 or 2, the dimension of WORK must be at least
*>                     N**2 + 2*N + floor(N/2) - 4.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array,
*>         If ICOMPQ = 0, IWORK is not referenced.
*>         If ICOMPQ = 1 or 2, the dimension of IWORK must be at least 4*N.
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
*> \ingroup laedk0
*
*> \par Contributors:
*  ==================
*>
*> Jeff Rutter, Computer Science Division, University of California
*> at Berkeley, USA
*
*  =====================================================================
      SUBROUTINE DLAEDK0( ICOMPQ, QSIZ, N, E, Q, LDQ, QSTORE, LDQS,
     $                    WORK, IWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   E( * ), Q( LDQ, * ), QSTORE( LDQS, * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.D0, ONE = 1.D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INODE, NDIML, NDIMR, SMLSIZ, SQRE, NSVD,
     $                   LVL, ND, IDXQ, IDXQC, IWK, I, J, I1,
     $                   IC, NDB1, NL, NLF, NLP1, NLVL, NR, NRF,
     $                   NRP1, LF, LL, SQREI, IM1, IDXQR, IDXD,
     $                   ITEMP
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DCOPY, DGEMM, DLACPY, DLASDT,
     $                   DLASET, DKTEQR, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT, ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.2 ) THEN
         INFO = -1
      ELSE IF( ( ICOMPQ.EQ.1 ) .AND. ( QSIZ.LT.MAX( 0, N ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDQS.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAEDK0', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( ICOMPQ.EQ.0 ) THEN
         CALL DKTEQR( 'N', N, E, Q, LDQ, WORK, INFO )
         RETURN
      END IF
      IF( ICOMPQ.EQ.1 ) THEN
         CALL DLACPY( 'A', QSIZ, N, Q, LDQ, QSTORE, LDQS )
         CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
      END IF
*
      SMLSIZ = ILAENV( 9, 'DLAEDK0', ' ', 0, 0, 0, 0 )
      NSVD = N/2
      SQRE = MOD(N, 2)
      IDXQR = N*N+1
      IDXD = IDXQR+2*N-4
*
      IF( NSVD.LE.SMLSIZ ) THEN
         CALL DKTEQR( 'I', N, E, Q, LDQ, WORK, INFO )
         GO TO 60
      END IF
*
*     Determine the size and placement of the submatrices, and save in
*     the leading elements of IWORK. Reuse DLASDT to set up the computation
*     tree.
*
      INODE = 1
      NDIML = INODE + NSVD
      NDIMR = NDIML + NSVD
      IDXQ = NDIMR + NSVD
      IWK = IDXQ + NSVD
      CALL DLASDT( NSVD, NLVL, ND, IWORK( INODE ), IWORK( NDIML ),
     $             IWORK( NDIMR ), SMLSIZ )
*
*     For the nodes on bottom level of the tree, solve
*     their subproblems by SSTEQR.
*
      NDB1 = ( ND+1 ) / 2
      DO 30 I = NDB1, ND
*
*        IC : center row of each node
*        NL : number of rows of left  subproblem
*        NR : number of rows of right subproblem
*        NLF: starting row of the left   subproblem
*        NRF: starting row of the right  subproblem
*
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NLP1 = NL + 1
         NR = IWORK( NDIMR+I1 )
         NRP1 = NR + 1
         NLF = IC - NL
         NRF = IC + 1
         SQREI = 1
         CALL DKTEQR( 'I', 2*NL+1, E( 2*NLF-1 ),
     $                WORK( (2*NLF-2)*N+2*NLF-1 ), N,
     $                WORK( IDXQR ), INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
*
*        Reorganize the output of DKTEQR into DLASD1's layout
*
         DO J = 1, NL
            IF( DDOT( NL+1, WORK( 2*NLF-1+(2*(NLF+J-1)-2)*N ),
     $          2, WORK( 2*NLF-1+(2*(NLF+J-1)-2)*N ), 2 ) .GT.
     $          DDOT( NL, WORK( 2*NLF+(2*(NLF+J-1)-2)*N ),
     $          2, WORK( 2*NLF+(2*(NLF+J-1)-2)*N ), 2 ) ) THEN
               WORK( IDXD+NLF+J-2 ) = E( 2*(NLF+J-1)-1 )
               CALL DCOPY( NL, WORK( 2*NLF+(2*(NLF+J-1)-1)*N ),
     $                     2, Q( NLF, NLF+J-1 ), 1 )
               CALL DCOPY( NL+1,
     $                     WORK( 2*NLF-1+(2*(NLF+J-1)-2)*N ),
     $                     2, Q( NSVD+NLF+J-1, NSVD+NLF ), LDQ )
            ELSE
               WORK( IDXD+NLF+J-2 ) = E( 2*(NLF+J-1)-1 )
               CALL DCOPY( NL, WORK( 2*NLF+(2*(NLF+J-1)-2)*N ),
     $                     2, Q( NLF, NLF+J-1 ), 1 )
               CALL DSCAL( NL, -ONE, Q( NLF, NLF+J-1 ), 1 )
               CALL DCOPY( NL+1,
     $                     WORK( 2*NLF-1+(2*(NLF+J-1)-1)*N ),
     $                     2, Q( NSVD+NLF+J-1, NSVD+NLF ), LDQ )
            END IF
         END DO
         CALL DCOPY( NL+1, WORK( 2*NLF-1+(2*(NLF+NL)-2)*N ),
     $               2, Q( NSVD+NLF+NL, NSVD+NLF ), LDQ )
*
*        End of reorganization
*
         ITEMP = IDXQ + NLF - 2
         DO 10 J = 1, NL
            IWORK( ITEMP+J ) = NL-J+1
   10    CONTINUE
         IF( I.EQ.ND ) THEN
            SQREI = SQRE
         ELSE
            SQREI = 1
         END IF
         NRP1 = NR + SQREI
         CALL DKTEQR( 'I', 2*NR+SQREI, E( 2*NRF-1 ),
     $                WORK( (2*NRF-2)*N+2*NRF-1 ), N,
     $                WORK( IDXQR ), INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
*
*        Reorganize the output of DKTEQR into DLASD1's layout
*
         DO J = 1, NR
            IF( DDOT( NR+SQREI, WORK( 2*NRF-1+(2*(NRF+J-1)-2)*N ),
     $          2, WORK( 2*NRF-1+(2*(NRF+J-1)-2)*N ), 2 ) .GT.
     $          DDOT( NR, WORK( 2*NRF+(2*(NRF+J-1)-2)*N ),
     $          2, WORK( 2*NRF+(2*(NRF+J-1)-2)*N ), 2 ) ) THEN
               WORK( IDXD+NRF+J-2 ) = E( 2*(NRF+J-1)-1 )
               CALL DCOPY( NR, WORK( 2*NRF+(2*(NRF+J-1)-1)*N ),
     $                     2, Q( NRF, NRF+J-1 ), 1 )
               CALL DCOPY( NR+SQREI,
     $                     WORK( 2*NRF-1+(2*(NRF+J-1)-2)*N ),
     $                     2, Q( NSVD+NRF+J-1, NSVD+NRF ), LDQ )
            ELSE
               WORK( IDXD+NRF+J-2 ) = E( 2*(NRF+J-1)-1 )
               CALL DCOPY( NR, WORK( 2*NRF+(2*(NRF+J-1)-2)*N ),
     $                     2, Q( NRF, NRF+J-1 ), 1 )
               CALL DSCAL( NR, -ONE, Q( NRF, NRF+J-1 ), 1 )
               CALL DCOPY( NR+SQREI,
     $                     WORK( 2*NRF-1+(2*(NRF+J-1)-1)*N ),
     $                     2, Q( NSVD+NRF+J-1, NSVD+NRF ), LDQ )
            END IF
         END DO
         IF( SQREI.EQ.1 ) THEN
            CALL DCOPY( NR+1, WORK( 2*NRF-1+(2*(NRF+NR)-2)*N ),
     $                  2, Q( NSVD+NRF+NR, NSVD+NRF ), LDQ )
         END IF
*
*        End of reorganization
*
         ITEMP = IDXQ + IC
         DO 20 J = 1, NR
            IWORK( ITEMP+J-1 ) = NR-J+1
   20    CONTINUE
   30 CONTINUE
*
*     Now conquer each subproblem bottom-up.
*
      DO 50 LVL = NLVL, 1, -1
*
*        Find the first node LF and last node LL on the
*        current level LVL.
*
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO 40 I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            IF( ( SQRE.EQ.0 ) .AND. ( I.EQ.LL ) ) THEN
               SQREI = SQRE
            ELSE
               SQREI = 1
            END IF
            IDXQC = IDXQ + NLF - 1
            ALPHA = E( 2*IC-1 )
            BETA = -E( 2*IC )
            CALL DLASD1( NL, NR, SQREI, WORK( IDXD+NLF-1 ), ALPHA,
     $                   BETA, Q( NLF, NLF ), LDQ,
     $                   Q( NSVD+NLF, NSVD+NLF ), LDQ,
     $                   IWORK( IDXQC ), IWORK( IWK ),
     $                   WORK, INFO )
*
*           Report the possible convergence failure.
*
            IF( INFO.NE.0 ) THEN
               RETURN
            END IF
   40    CONTINUE
   50 CONTINUE
*
*     Copy results in Q to WORK
*
      CALL DLACPY( 'A', NSVD, NSVD, Q, LDQ, WORK, NSVD )
      CALL DLASET( 'Full', NSVD, NSVD, ZERO, ZERO, Q, LDQ )
      CALL DLACPY( 'A', NSVD+SQRE, NSVD+SQRE, Q(NSVD+1, NSVD+1),
     $             LDQ, WORK(NSVD*NSVD+1), NSVD+SQRE )
      CALL DLASET( 'Full', NSVD+SQRE, NSVD+SQRE, ZERO, ZERO,
     $             Q(NSVD+1, NSVD+1), LDQ )
*
*     Reorganize the output of DLASD1 into DKTEQR's layout
*
      DO J = 1, NSVD
         E( 2*J-1 ) = WORK( IDXD+IWORK( IDXQC+NSVD-J )-1 )
         IF ( SQRE.EQ.1 .OR. J.NE.NSVD ) THEN
            E( 2*J ) = ZERO
         END IF
         CALL DCOPY( NSVD+SQRE,
     $               WORK( NSVD*NSVD+IWORK( IDXQC+NSVD-J ) ),
     $               NSVD+SQRE, Q( 1, 2*J-1 ), 2 )
         CALL DCOPY( NSVD,
     $               WORK( (IWORK( IDXQC+NSVD-J )-1)*NSVD+1 ),
     $               1, Q( 2, 2*J ), 2 )
      END DO
      IF( SQRE.EQ.1 ) THEN
         CALL DCOPY( NSVD+1, WORK( NSVD*NSVD+NSVD+1 ), NSVD+SQRE,
     $               Q( 1, N ), 2 )
      END IF
*
*     End of reorganization
*
*     If ICOMPQ = 1, build the final matrix
*
   60 CONTINUE
      IF( ICOMPQ.EQ.1 ) THEN
         CALL DLACPY( 'A', N, N, Q, LDQ, WORK, N )
         CALL DGEMM( 'N', 'N', QSIZ, N, N, ONE, QSTORE, LDQS, WORK,
     $               N, ZERO, Q, LDQ )
      END IF
*
      RETURN
*
*     End of DLAEDK0
*
      END
