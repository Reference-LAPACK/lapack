*> \brief \b SLASDT creates a tree of subproblems for bidiagonal divide and conquer. Used by sbdsdc.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SLASDT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasdt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasdt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasdt.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )
*
*       .. Scalar Arguments ..
*       INTEGER            LVL, MSUB, N, ND
*       ..
*       .. Array Arguments ..
*       INTEGER            INODE( * ), NDIML( * ), NDIMR( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SLASDT creates a tree of subproblems for bidiagonal divide and
*> conquer.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          On entry, the number of diagonal elements of the
*>          bidiagonal matrix.
*> \endverbatim
*>
*> \param[out] LVL
*> \verbatim
*>          LVL is INTEGER
*>          On exit, the number of levels on the computation tree.
*> \endverbatim
*>
*> \param[out] ND
*> \verbatim
*>          ND is INTEGER
*>          On exit, the number of nodes on the tree.
*> \endverbatim
*>
*> \param[out] INODE
*> \verbatim
*>          INODE is INTEGER array, dimension ( N )
*>          On exit, centers of subproblems.
*> \endverbatim
*>
*> \param[out] NDIML
*> \verbatim
*>          NDIML is INTEGER array, dimension ( N )
*>          On exit, row dimensions of left children.
*> \endverbatim
*>
*> \param[out] NDIMR
*> \verbatim
*>          NDIMR is INTEGER array, dimension ( N )
*>          On exit, row dimensions of right children.
*> \endverbatim
*>
*> \param[in] MSUB
*> \verbatim
*>          MSUB is INTEGER
*>          On entry, the maximum row dimension each subproblem at the
*>          bottom of the tree can be of.
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
*> \ingroup lasdt
*
*> \par Contributors:
*  ==================
*>
*>     Ming Gu and Huan Ren, Computer Science Division, University of
*>     California at Berkeley, USA
*>
*  =====================================================================
      SUBROUTINE SLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LVL, MSUB, N, ND
*     ..
*     .. Array Arguments ..
      INTEGER            INODE( * ), NDIML( * ), NDIMR( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               TWO
      PARAMETER          ( TWO = 2.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IL, IR, LLST, MAXN, NCRNT, NLVL
      REAL               TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, LOG, MAX, REAL
*     ..
*     .. Executable Statements ..
*
*     Find the number of levels on the tree.
*
      MAXN = MAX( 1, N )
      TEMP = LOG( REAL( MAXN ) / REAL( MSUB+1 ) ) / LOG( TWO )
      LVL = INT( TEMP ) + 1
*
      I = N / 2
      INODE( 1 ) = I + 1
      NDIML( 1 ) = I
      NDIMR( 1 ) = N - I - 1
      IL = 0
      IR = 1
      LLST = 1
      DO 20 NLVL = 1, LVL - 1
*
*        Constructing the tree at (NLVL+1)-st level. The number of
*        nodes created on this level is LLST * 2.
*
         DO 10 I = 0, LLST - 1
            IL = IL + 2
            IR = IR + 2
            NCRNT = LLST + I
            NDIML( IL ) = NDIML( NCRNT ) / 2
            NDIMR( IL ) = NDIML( NCRNT ) - NDIML( IL ) - 1
            INODE( IL ) = INODE( NCRNT ) - NDIMR( IL ) - 1
            NDIML( IR ) = NDIMR( NCRNT ) / 2
            NDIMR( IR ) = NDIMR( NCRNT ) - NDIML( IR ) - 1
            INODE( IR ) = INODE( NCRNT ) + NDIML( IR ) + 1
   10    CONTINUE
         LLST = LLST*2
   20 CONTINUE
      ND = LLST*2 - 1
*
      RETURN
*
*     End of SLASDT
*
      END
