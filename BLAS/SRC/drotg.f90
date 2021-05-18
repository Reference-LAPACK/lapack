!> \brief \b DROTG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!  DROTG constructs a plane rotation
!     [  c  s ] [ a ] = [ r ]
!     [ -s  c ] [ b ]   [ 0 ]
!  satisfying c**2 + s**2 = 1.
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> The computation uses the formulas
!>    sigma = sgn(a)    if |a| >  |b|
!>          = sgn(b)    if |b| >= |a|
!>    r = sigma*sqrt( a**2 + b**2 )
!>    c = 1; s = 0      if r = 0
!>    c = a/r; s = b/r  if r != 0
!> The subroutine also computes
!>    z = s    if |a| > |b|,
!>      = 1/c  if |b| >= |a| and c != 0
!>      = 1    if c = 0
!> This allows c and s to be reconstructed from z as follows:
!>    If z = 1, set c = 0, s = 1.
!>    If |z| < 1, set c = sqrt(1 - z**2) and s = z.
!>    If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION
!>          On entry, the scalar a.
!>          On exit, the scalar r.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION
!>          On entry, the scalar b.
!>          On exit, the scalar z.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION
!>          The scalar c.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION
!>          The scalar s.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \ingroup single_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!> \endverbatim
!
!  =====================================================================
subroutine DROTG( a, b, c, s )
   integer, parameter :: wp = kind(1.d0)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  .. Constants ..
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one  = 1.0_wp
!  ..
!  .. Scaling constants ..
   real(wp), parameter :: safmin = real(radix(0._wp),wp)**max( &
      minexponent(0._wp)-1, &
      1-maxexponent(0._wp) &
   )
   real(wp), parameter :: safmax = real(radix(0._wp),wp)**max( &
      1-minexponent(0._wp), &
      maxexponent(0._wp)-1 &
   )
!  ..
!  .. Scalar Arguments ..
   real(wp) :: a, b, c, s
!  ..
!  .. Local Scalars ..
   real(wp) :: anorm, bnorm, scl, sigma, r, z
!  ..
   anorm = abs(a)
   bnorm = abs(b)
   if( bnorm == zero ) then
      c = one
      s = zero
      b = zero
   else if( anorm == zero ) then
      c = zero
      s = one
      a = b
      b = one
   else
      scl = min( safmax, max( safmin, anorm, bnorm ) )
      if( anorm > bnorm ) then
         sigma = sign(one,a)
      else
         sigma = sign(one,b)
      end if
      r = sigma*( scl*sqrt((a/scl)**2 + (b/scl)**2) )
      c = a/r
      s = b/r
      if( anorm > bnorm ) then
         z = s
      else if( c /= zero ) then
         z = one/c
      else
         z = one
      end if
      a = r
      b = z
   end if
   return
end subroutine
