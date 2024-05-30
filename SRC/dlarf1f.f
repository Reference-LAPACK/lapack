*> \brief \b DLARF1F applies an elementary reflector to a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLARF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE
*       INTEGER            INCV, LDC, M, N
*       DOUBLE PRECISION   TAU
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARF applies a real elementary reflector H to a real m by n matrix
*> C, from either the left or the right. H is represented in the form
*>
*>       H = I - tau * v * v**T
*>
*> where tau is a real scalar and v is a real vector.
*>
*> If tau = 0, then H is taken to be the unit matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': form  H * C
*>          = 'R': form  C * H
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension
*>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*>          The vector v in the representation of H. V is not used if
*>          TAU = 0.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          The increment between elements of v. INCV <> 0.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION
*>          The value tau in the representation of H.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
*>          On entry, the m by n matrix C.
*>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*>          or C * H if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension
*>                         (N) if SIDE = 'L'
*>                      or (M) if SIDE = 'R'
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
*> \ingroup larf
*
*  =====================================================================
      SUBROUTINE DLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
      DOUBLE PRECISION   C11, DOT1, DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DDOT, DAXPY, DCOPY, DSCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILADLR, ILADLC
      EXTERNAL           LSAME, ILADLR, ILADLC
*     ..
*     .. Executable Statements ..
*
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILADLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILADLR(M, LASTV, C, LDC)
         END IF
      END IF

      IF( LASTC.EQ.0 .OR. LASTV.EQ.0 ) THEN
         RETURN
      END IF

      IF( APPLYLEFT ) THEN
*
*        Form  H * C
*
         IF( LASTV.EQ.1 ) THEN
            CALL DSCAL(LASTC, ONE - TAU, C, LDC)
         ELSE
            DOT1 = - TAU * DDOT( LASTV - 1, V( 1 + INCV ), INCV,
     $            C( 2, 1 ), 1 )

            C11 = (ONE - TAU) * C( 1, 1 ) + DOT1
*
*           Prepare WORK
*
            CALL DCOPY( LASTC - 1, C( 1, 2 ), LDC, WORK, 1 )

            CALL DGEMV( 'Transpose', LASTV - 1, LASTC - 1, -TAU,
     $            C( 2, 2 ), LDC, V( 1 + INCV ), INCV, -TAU, WORK, 1 )
*
*           Update C12
*
            CALL DAXPY( LASTC - 1, ONE, WORK, 1, C( 1, 2 ), LDC )
*
*           Update C21
*
            CALL DAXPY( LASTV - 1, -TAU * C( 1, 1 ) + DOT1,
     $            V( 1 + INCV ), INCV, C( 2, 1 ), 1 )
*
*           Update C11
*
            C( 1, 1 ) = C11
*
*           Update C22
*
            CALL DGER( LASTV - 1, LASTC - 1, ONE, V( 1 + INCV ),
     $            INCV, WORK, 1, C( 2, 2 ), LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( LASTV.EQ.1 ) THEN
            CALL DSCAL(LASTC, ONE - TAU, C, 1)
         ELSE
            DOT1 = - TAU * DDOT( LASTV - 1, V( 1 + INCV ), INCV,
     $            C( 1, 2 ), LDC )

            C11 = (ONE - TAU) * C( 1, 1 ) + DOT1
*
*           Prepare WORK
*
            CALL DCOPY( LASTC - 1, C( 2, 1 ), 1, WORK, 1 )

            CALL DGEMV( 'No transpose', LASTC - 1, LASTV - 1, -TAU,
     $            C( 2, 2 ), LDC, V( 1 + INCV ), INCV, -TAU, WORK, 1 )
*
*           Update C12
*
            CALL DAXPY( LASTV - 1, -TAU * C( 1, 1 ) + DOT1,
     $            V( 1 + INCV ), INCV, C( 1, 2 ), LDC )
*
*           Update C11
*
            C( 1, 1 ) = C11
*
*           Update C21
*
            CALL DAXPY( LASTC - 1, ONE, WORK, 1, C( 2, 1 ), 1 )
*
*           Update C22
*
            CALL DGER( LASTC - 1, LASTV - 1, ONE, WORK, 1,
     $            V( 1 + INCV ), INCV, C( 2, 2 ), LDC )
         END IF
      END IF
      RETURN
*
*     End of DLARF1F
*
      END
