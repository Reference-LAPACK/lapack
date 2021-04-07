!> \brief \b CROTG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!  CROTG constructs a plane rotation
!     [  c         s ] [ a ] = [ r ]
!     [ -conjg(s)  c ] [ b ]   [ 0 ]
!  where c is real, s ic complex, and c**2 + conjg(s)*s = 1.
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> The computation uses the formulas
!>    |x| = sqrt( Re(x)**2 + Im(x)**2 )
!>    sgn(x) = x / |x|  if x /= 0
!>           = 1        if x  = 0
!>    c = |a| / sqrt(|a|**2 + |b|**2)
!>    s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
!> When a and b are real and r /= 0, the formulas simplify to
!>    r = sgn(a)*sqrt(|a|**2 + |b|**2)
!>    c = a / r
!>    s = b / r
!> the same as in CROTG when |a| > |b|.  When |b| >= |a|, the
!> sign of c and s will be different from those computed by CROTG
!> if the signs of a and b are not the same.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX
!>          On entry, the scalar a.
!>          On exit, the scalar r.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX
!>          The scalar b.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL
!>          The scalar c.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL
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
subroutine CROTG( a, b, c, s )
   integer, parameter :: wp = kind(1.e0)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  .. Constants ..
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one  = 1.0_wp
   complex(wp), parameter :: czero  = 0.0_wp
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
   real(wp), parameter :: rtmin = sqrt( real(radix(0._wp),wp)**max( &
      minexponent(0._wp)-1, &
      1-maxexponent(0._wp) &
   ) / epsilon(0._wp) )
   real(wp), parameter :: rtmax = sqrt( real(radix(0._wp),wp)**max( &
      1-minexponent(0._wp), &
      maxexponent(0._wp)-1 &
   ) * epsilon(0._wp) )
!  ..
!  .. Scalar Arguments ..
   real(wp) :: c
   complex(wp) :: a, b, s
!  ..
!  .. Local Scalars ..
   real(wp) :: d, f1, f2, g1, g2, h2, p, u, uu, v, vv, w
   complex(wp) :: f, fs, g, gs, r, t
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, aimag, conjg, max, min, real, sqrt
!  ..
!  .. Statement Functions ..
   real(wp) :: ABSSQ
!  ..
!  .. Statement Function definitions ..
   ABSSQ( t ) = real( t )**2 + aimag( t )**2
!  ..
!  .. Executable Statements ..
!
   f = a
   g = b
   if( g == czero ) then
      c = one
      s = czero
      r = f
   else if( f == czero ) then
      c = zero
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         g2 = ABSSQ( g )
         d = sqrt( g2 )
         s = conjg( g ) / d
         r = d
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, g1 ) )
         uu = one / u
         gs = g*uu
         g2 = ABSSQ( gs )
         d = sqrt( g2 )
         s = conjg( gs ) / d
         r = d*u
      end if
   else
      f1 = max( abs(real(f)), abs(aimag(f)) )
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( f1 > rtmin .and. f1 < rtmax .and. &
          g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         f2 = ABSSQ( f )
         g2 = ABSSQ( g )
         h2 = f2 + g2
         if( f2 > rtmin .and. h2 < rtmax ) then
            d = sqrt( f2*h2 )
         else
            d = sqrt( f2 )*sqrt( h2 )
         end if
         p = 1 / d
         c = f2*p
         s = conjg( g )*( f*p )
         r = f*( h2*p )
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, f1, g1 ) )
         uu = one / u
         gs = g*uu
         g2 = ABSSQ( gs )
         if( f1*uu < rtmin ) then
!
!           f is not well-scaled when scaled by g1.
!           Use a different scaling for f.
!
            v = min( safmax, max( safmin, f1 ) )
            vv = one / v
            w = v * uu
            fs = f*vv
            f2 = ABSSQ( fs )
            h2 = f2*w**2 + g2
         else
!
!           Otherwise use the same scaling for f and g.
!
            w = one
            fs = f*uu
            f2 = ABSSQ( fs )
            h2 = f2 + g2
         end if
         if( f2 > rtmin .and. h2 < rtmax ) then
            d = sqrt( f2*h2 )
         else
            d = sqrt( f2 )*sqrt( h2 )
         end if
         p = 1 / d
         c = ( f2*p )*w
         s = conjg( gs )*( fs*p )
         r = ( fs*( h2*p ) )*u
      end if
   end if
   a = r
   return
end subroutine
