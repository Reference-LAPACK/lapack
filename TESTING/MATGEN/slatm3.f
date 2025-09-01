*> \brief \b SLATM3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       REAL             FUNCTION SLATM3( M, N, I, J, ISUB, JSUB, KL, KU,
*                        IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK,
*                        SPARSE )
*
*       .. Scalar Arguments ..
*
*       INTEGER            I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL,
*      $                   KU, M, N
*       REAL               SPARSE
*       ..
*
*       .. Array Arguments ..
*
*       INTEGER            ISEED( 4 ), IWORK( * )
*       REAL               D( * ), DL( * ), DR( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    SLATM3 returns the (ISUB,JSUB) entry of a random matrix of
*>    dimension (M, N) described by the other parameters. (ISUB,JSUB)
*>    is the final position of the (I,J) entry after pivoting
*>    according to IPVTNG and IWORK. SLATM3 is called by the
*>    SLATMR routine in order to build random test matrices. No error
*>    checking on parameters is done, because this routine is called in
*>    a tight loop by SLATMR which has already checked the parameters.
*>
*>    Use of SLATM3 differs from SLATM2 in the order in which the random
*>    number generator is called to fill in random matrix entries.
*>    With SLATM2, the generator is called to fill in the pivoted matrix
*>    columnwise. With SLATM3, the generator is called to fill in the
*>    matrix columnwise, after which it is pivoted. Thus, SLATM3 can
*>    be used to construct random matrices which differ only in their
*>    order of rows and/or columns. SLATM2 is used to construct band
*>    matrices while avoiding calling the random number generator for
*>    entries outside the band (and therefore generating random numbers
*>    in different orders for different pivot orders).
*>
*>    The matrix whose (ISUB,JSUB) entry is returned is constructed as
*>    follows (this routine only computes one entry):
*>
*>      If ISUB is outside (1..M) or JSUB is outside (1..N), return zero
*>         (this is convenient for generating matrices in band format).
*>
*>      Generate a matrix A with random entries of distribution IDIST.
*>
*>      Set the diagonal to D.
*>
*>      Grade the matrix, if desired, from the left (by DL) and/or
*>         from the right (by DR or DL) as specified by IGRADE.
*>
*>      Permute, if desired, the rows and/or columns as specified by
*>         IPVTNG and IWORK.
*>
*>      Band the matrix to have lower bandwidth KL and upper
*>         bandwidth KU.
*>
*>      Set random entries to zero as specified by SPARSE.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           Number of rows of matrix. Not modified.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           Number of columns of matrix. Not modified.
*> \endverbatim
*>
*> \param[in] I
*> \verbatim
*>          I is INTEGER
*>           Row of unpivoted entry to be returned. Not modified.
*> \endverbatim
*>
*> \param[in] J
*> \verbatim
*>          J is INTEGER
*>           Column of unpivoted entry to be returned. Not modified.
*> \endverbatim
*>
*> \param[in,out] ISUB
*> \verbatim
*>          ISUB is INTEGER
*>           Row of pivoted entry to be returned. Changed on exit.
*> \endverbatim
*>
*> \param[in,out] JSUB
*> \verbatim
*>          JSUB is INTEGER
*>           Column of pivoted entry to be returned. Changed on exit.
*> \endverbatim
*>
*> \param[in] KL
*> \verbatim
*>          KL is INTEGER
*>           Lower bandwidth. Not modified.
*> \endverbatim
*>
*> \param[in] KU
*> \verbatim
*>          KU is INTEGER
*>           Upper bandwidth. Not modified.
*> \endverbatim
*>
*> \param[in] IDIST
*> \verbatim
*>          IDIST is INTEGER
*>           On entry, IDIST specifies the type of distribution to be
*>           used to generate a random matrix .
*>           1 => UNIFORM( 0, 1 )
*>           2 => UNIFORM( -1, 1 )
*>           3 => NORMAL( 0, 1 )
*>           Not modified.
*> \endverbatim
*>
*> \param[in,out] ISEED
*> \verbatim
*>          ISEED is INTEGER array of dimension ( 4 )
*>           Seed for random number generator.
*>           Changed on exit.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is REAL array of dimension ( MIN( I , J ) )
*>           Diagonal entries of matrix. Not modified.
*> \endverbatim
*>
*> \param[in] IGRADE
*> \verbatim
*>          IGRADE is INTEGER
*>           Specifies grading of matrix as follows:
*>           0  => no grading
*>           1  => matrix premultiplied by diag( DL )
*>           2  => matrix postmultiplied by diag( DR )
*>           3  => matrix premultiplied by diag( DL ) and
*>                         postmultiplied by diag( DR )
*>           4  => matrix premultiplied by diag( DL ) and
*>                         postmultiplied by inv( diag( DL ) )
*>           5  => matrix premultiplied by diag( DL ) and
*>                         postmultiplied by diag( DL )
*>           Not modified.
*> \endverbatim
*>
*> \param[in] DL
*> \verbatim
*>          DL is REAL array ( I or J, as appropriate )
*>           Left scale factors for grading matrix.  Not modified.
*> \endverbatim
*>
*> \param[in] DR
*> \verbatim
*>          DR is REAL array ( I or J, as appropriate )
*>           Right scale factors for grading matrix.  Not modified.
*> \endverbatim
*>
*> \param[in] IPVTNG
*> \verbatim
*>          IPVTNG is INTEGER
*>           On entry specifies pivoting permutations as follows:
*>           0 => none.
*>           1 => row pivoting.
*>           2 => column pivoting.
*>           3 => full pivoting, i.e., on both sides.
*>           Not modified.
*> \endverbatim
*>
*> \param[in] IWORK
*> \verbatim
*>          IWORK is INTEGER array ( I or J, as appropriate )
*>           This array specifies the permutation used. The
*>           row (or column) originally in position K is in
*>           position IWORK( K ) after pivoting.
*>           This differs from IWORK for SLATM2. Not modified.
*> \endverbatim
*>
*> \param[in] SPARSE
*> \verbatim
*>          SPARSE is REAL between 0. and 1.
*>           On entry specifies the sparsity of the matrix
*>           if sparse matrix is to be generated.
*>           SPARSE should lie between 0 and 1.
*>           A uniform ( 0, 1 ) random number x is generated and
*>           compared to SPARSE; if x is larger the matrix entry
*>           is unchanged and if x is smaller the entry is set
*>           to zero. Thus on the average a fraction SPARSE of the
*>           entries will be set to zero.
*>           Not modified.
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
*> \ingroup real_matgen
*
*  =====================================================================
      REAL             FUNCTION SLATM3( M, N, I, J, ISUB, JSUB, KL,
     $                 KU,
     $                 IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK,
     $                 SPARSE )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
*
      INTEGER            I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL,
     $                   KU, M, N
      REAL               SPARSE
*     ..
*
*     .. Array Arguments ..
*
      INTEGER            ISEED( 4 ), IWORK( * )
      REAL               D( * ), DL( * ), DR( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*
*     .. Local Scalars ..
*
      REAL               TEMP
*     ..
*
*     .. External Functions ..
*
      REAL               SLARAN, SLARND
      EXTERNAL           SLARAN, SLARND
*     ..
*
*-----------------------------------------------------------------------
*
*     .. Executable Statements ..
*
*
*     Check for I and J in range
*
      IF( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) THEN
         ISUB = I
         JSUB = J
         SLATM3 = ZERO
         RETURN
      END IF
*
*     Compute subscripts depending on IPVTNG
*
      IF( IPVTNG.EQ.0 ) THEN
         ISUB = I
         JSUB = J
      ELSE IF( IPVTNG.EQ.1 ) THEN
         ISUB = IWORK( I )
         JSUB = J
      ELSE IF( IPVTNG.EQ.2 ) THEN
         ISUB = I
         JSUB = IWORK( J )
      ELSE IF( IPVTNG.EQ.3 ) THEN
         ISUB = IWORK( I )
         JSUB = IWORK( J )
      END IF
*
*     Check for banding
*
      IF( JSUB.GT.ISUB+KU .OR. JSUB.LT.ISUB-KL ) THEN
         SLATM3 = ZERO
         RETURN
      END IF
*
*     Check for sparsity
*
      IF( SPARSE.GT.ZERO ) THEN
         IF( SLARAN( ISEED ).LT.SPARSE ) THEN
            SLATM3 = ZERO
            RETURN
         END IF
      END IF
*
*     Compute entry and grade it according to IGRADE
*
      IF( I.EQ.J ) THEN
         TEMP = D( I )
      ELSE
         TEMP = SLARND( IDIST, ISEED )
      END IF
      IF( IGRADE.EQ.1 ) THEN
         TEMP = TEMP*DL( I )
      ELSE IF( IGRADE.EQ.2 ) THEN
         TEMP = TEMP*DR( J )
      ELSE IF( IGRADE.EQ.3 ) THEN
         TEMP = TEMP*DL( I )*DR( J )
      ELSE IF( IGRADE.EQ.4 .AND. I.NE.J ) THEN
         TEMP = TEMP*DL( I ) / DL( J )
      ELSE IF( IGRADE.EQ.5 ) THEN
         TEMP = TEMP*DL( I )*DL( J )
      END IF
      SLATM3 = TEMP
      RETURN
*
*     End of SLATM3
*
      END
