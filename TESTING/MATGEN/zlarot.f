*> \brief \b ZLAROT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAROT( LROWS, LLEFT, LRIGHT, NL, C, S, A, LDA, XLEFT,
*                          XRIGHT )
*
*       .. Scalar Arguments ..
*       LOGICAL            LLEFT, LRIGHT, LROWS
*       INTEGER            LDA, NL
*       COMPLEX*16         C, S, XLEFT, XRIGHT
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZLAROT applies a (Givens) rotation to two adjacent rows or
*>    columns, where one element of the first and/or last column/row
*>    for use on matrices stored in some format other than GE, so
*>    that elements of the matrix may be used or modified for which
*>    no array element is provided.
*>
*>    One example is a symmetric matrix in SB format (bandwidth=4), for
*>    which UPLO='L':  Two adjacent rows will have the format:
*>
*>    row j:     C> C> C> C> C> .  .  .  .
*>    row j+1:      C> C> C> C> C> .  .  .  .
*>
*>    '*' indicates elements for which storage is provided,
*>    '.' indicates elements for which no storage is provided, but
*>    are not necessarily zero; their values are determined by
*>    symmetry.  ' ' indicates elements which are necessarily zero,
*>     and have no storage provided.
*>
*>    Those columns which have two '*'s can be handled by DROT.
*>    Those columns which have no '*'s can be ignored, since as long
*>    as the Givens rotations are carefully applied to preserve
*>    symmetry, their values are determined.
*>    Those columns which have one '*' have to be handled separately,
*>    by using separate variables "p" and "q":
*>
*>    row j:     C> C> C> C> C> p  .  .  .
*>    row j+1:   q  C> C> C> C> C> .  .  .  .
*>
*>    The element p would have to be set correctly, then that column
*>    is rotated, setting p to its new value.  The next call to
*>    ZLAROT would rotate columns j and j+1, using p, and restore
*>    symmetry.  The element q would start out being zero, and be
*>    made non-zero by the rotation.  Later, rotations would presumably
*>    be chosen to zero q out.
*>
*>    Typical Calling Sequences: rotating the i-th and (i+1)-st rows.
*>    ------- ------- ---------
*>
*>      General dense matrix:
*>
*>              CALL ZLAROT(.TRUE.,.FALSE.,.FALSE., N, C,S,
*>                      A(i,1),LDA, DUMMY, DUMMY)
*>
*>      General banded matrix in GB format:
*>
*>              j = MAX(1, i-KL )
*>              NL = MIN( N, i+KU+1 ) + 1-j
*>              CALL ZLAROT( .TRUE., i-KL.GE.1, i+KU.LT.N, NL, C,S,
*>                      A(KU+i+1-j,j),LDA-1, XLEFT, XRIGHT )
*>
*>              [ note that i+1-j is just MIN(i,KL+1) ]
*>
*>      Symmetric banded matrix in SY format, bandwidth K,
*>      lower triangle only:
*>
*>              j = MAX(1, i-K )
*>              NL = MIN( K+1, i ) + 1
*>              CALL ZLAROT( .TRUE., i-K.GE.1, .TRUE., NL, C,S,
*>                      A(i,j), LDA, XLEFT, XRIGHT )
*>
*>      Same, but upper triangle only:
*>
*>              NL = MIN( K+1, N-i ) + 1
*>              CALL ZLAROT( .TRUE., .TRUE., i+K.LT.N, NL, C,S,
*>                      A(i,i), LDA, XLEFT, XRIGHT )
*>
*>      Symmetric banded matrix in SB format, bandwidth K,
*>      lower triangle only:
*>
*>              [ same as for SY, except:]
*>                  . . . .
*>                      A(i+1-j,j), LDA-1, XLEFT, XRIGHT )
*>
*>              [ note that i+1-j is just MIN(i,K+1) ]
*>
*>      Same, but upper triangle only:
*>                  . . .
*>                      A(K+1,i), LDA-1, XLEFT, XRIGHT )
*>
*>      Rotating columns is just the transpose of rotating rows, except
*>      for GB and SB: (rotating columns i and i+1)
*>
*>      GB:
*>              j = MAX(1, i-KU )
*>              NL = MIN( N, i+KL+1 ) + 1-j
*>              CALL ZLAROT( .TRUE., i-KU.GE.1, i+KL.LT.N, NL, C,S,
*>                      A(KU+j+1-i,i),LDA-1, XTOP, XBOTTM )
*>
*>              [note that KU+j+1-i is just MAX(1,KU+2-i)]
*>
*>      SB: (upper triangle)
*>
*>                   . . . . . .
*>                      A(K+j+1-i,i),LDA-1, XTOP, XBOTTM )
*>
*>      SB: (lower triangle)
*>
*>                   . . . . . .
*>                      A(1,i),LDA-1, XTOP, XBOTTM )
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \verbatim
*>  LROWS  - LOGICAL
*>           If .TRUE., then ZLAROT will rotate two rows.  If .FALSE.,
*>           then it will rotate two columns.
*>           Not modified.
*>
*>  LLEFT  - LOGICAL
*>           If .TRUE., then XLEFT will be used instead of the
*>           corresponding element of A for the first element in the
*>           second row (if LROWS=.FALSE.) or column (if LROWS=.TRUE.)
*>           If .FALSE., then the corresponding element of A will be
*>           used.
*>           Not modified.
*>
*>  LRIGHT - LOGICAL
*>           If .TRUE., then XRIGHT will be used instead of the
*>           corresponding element of A for the last element in the
*>           first row (if LROWS=.FALSE.) or column (if LROWS=.TRUE.) If
*>           .FALSE., then the corresponding element of A will be used.
*>           Not modified.
*>
*>  NL     - INTEGER
*>           The length of the rows (if LROWS=.TRUE.) or columns (if
*>           LROWS=.FALSE.) to be rotated.  If XLEFT and/or XRIGHT are
*>           used, the columns/rows they are in should be included in
*>           NL, e.g., if LLEFT = LRIGHT = .TRUE., then NL must be at
*>           least 2.  The number of rows/columns to be rotated
*>           exclusive of those involving XLEFT and/or XRIGHT may
*>           not be negative, i.e., NL minus how many of LLEFT and
*>           LRIGHT are .TRUE. must be at least zero; if not, XERBLA
*>           will be called.
*>           Not modified.
*>
*>  C, S   - COMPLEX*16
*>           Specify the Givens rotation to be applied.  If LROWS is
*>           true, then the matrix ( c  s )
*>                                 ( _  _ )
*>                                 (-s  c )  is applied from the left;
*>           if false, then the transpose (not conjugated) thereof is
*>           applied from the right.  Note that in contrast to the
*>           output of ZROTG or to most versions of ZROT, both C and S
*>           are complex.  For a Givens rotation, |C|**2 + |S|**2 should
*>           be 1, but this is not checked.
*>           Not modified.
*>
*>  A      - COMPLEX*16 array.
*>           The array containing the rows/columns to be rotated.  The
*>           first element of A should be the upper left element to
*>           be rotated.
*>           Read and modified.
*>
*>  LDA    - INTEGER
*>           The "effective" leading dimension of A.  If A contains
*>           a matrix stored in GE, HE, or SY format, then this is just
*>           the leading dimension of A as dimensioned in the calling
*>           routine.  If A contains a matrix stored in band (GB, HB, or
*>           SB) format, then this should be *one less* than the leading
*>           dimension used in the calling routine.  Thus, if A were
*>           dimensioned A(LDA,*) in ZLAROT, then A(1,j) would be the
*>           j-th element in the first of the two rows to be rotated,
*>           and A(2,j) would be the j-th in the second, regardless of
*>           how the array may be stored in the calling routine.  [A
*>           cannot, however, actually be dimensioned thus, since for
*>           band format, the row number may exceed LDA, which is not
*>           legal FORTRAN.]
*>           If LROWS=.TRUE., then LDA must be at least 1, otherwise
*>           it must be at least NL minus the number of .TRUE. values
*>           in XLEFT and XRIGHT.
*>           Not modified.
*>
*>  XLEFT  - COMPLEX*16
*>           If LLEFT is .TRUE., then XLEFT will be used and modified
*>           instead of A(2,1) (if LROWS=.TRUE.) or A(1,2)
*>           (if LROWS=.FALSE.).
*>           Read and modified.
*>
*>  XRIGHT - COMPLEX*16
*>           If LRIGHT is .TRUE., then XRIGHT will be used and modified
*>           instead of A(1,NL) (if LROWS=.TRUE.) or A(NL,1)
*>           (if LROWS=.FALSE.).
*>           Read and modified.
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
*> \ingroup complex16_matgen
*
*  =====================================================================
      SUBROUTINE ZLAROT( LROWS, LLEFT, LRIGHT, NL, C, S, A, LDA,
     $                   XLEFT,
     $                   XRIGHT )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            LLEFT, LRIGHT, LROWS
      INTEGER            LDA, NL
      COMPLEX*16         C, S, XLEFT, XRIGHT
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            IINC, INEXT, IX, IY, IYT, J, NT
      COMPLEX*16         TEMPX
*     ..
*     .. Local Arrays ..
      COMPLEX*16         XT( 2 ), YT( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. Executable Statements ..
*
*     Set up indices, arrays for ends
*
      IF( LROWS ) THEN
         IINC = LDA
         INEXT = 1
      ELSE
         IINC = 1
         INEXT = LDA
      END IF
*
      IF( LLEFT ) THEN
         NT = 1
         IX = 1 + IINC
         IY = 2 + LDA
         XT( 1 ) = A( 1 )
         YT( 1 ) = XLEFT
      ELSE
         NT = 0
         IX = 1
         IY = 1 + INEXT
      END IF
*
      IF( LRIGHT ) THEN
         IYT = 1 + INEXT + ( NL-1 )*IINC
         NT = NT + 1
         XT( NT ) = XRIGHT
         YT( NT ) = A( IYT )
      END IF
*
*     Check for errors
*
      IF( NL.LT.NT ) THEN
         CALL XERBLA( 'ZLAROT', 4 )
         RETURN
      END IF
      IF( LDA.LE.0 .OR. ( .NOT.LROWS .AND. LDA.LT.NL-NT ) ) THEN
         CALL XERBLA( 'ZLAROT', 8 )
         RETURN
      END IF
*
*     Rotate
*
*     ZROT( NL-NT, A(IX),IINC, A(IY),IINC, C, S ) with complex C, S
*
      DO 10 J = 0, NL - NT - 1
         TEMPX = C*A( IX+J*IINC ) + S*A( IY+J*IINC )
         A( IY+J*IINC ) = -DCONJG( S )*A( IX+J*IINC ) +
     $                    DCONJG( C )*A( IY+J*IINC )
         A( IX+J*IINC ) = TEMPX
   10 CONTINUE
*
*     ZROT( NT, XT,1, YT,1, C, S ) with complex C, S
*
      DO 20 J = 1, NT
         TEMPX = C*XT( J ) + S*YT( J )
         YT( J ) = -DCONJG( S )*XT( J ) + DCONJG( C )*YT( J )
         XT( J ) = TEMPX
   20 CONTINUE
*
*     Stuff values back into XLEFT, XRIGHT, etc.
*
      IF( LLEFT ) THEN
         A( 1 ) = XT( 1 )
         XLEFT = YT( 1 )
      END IF
*
      IF( LRIGHT ) THEN
         XRIGHT = XT( NT )
         A( IYT ) = YT( NT )
      END IF
*
      RETURN
*
*     End of ZLAROT
*
      END
