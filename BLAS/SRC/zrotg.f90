!> \brief \b ZROTG  generates a Givens rotation with real cosine and complex sine.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZROTG constructs a plane rotation
!>    [  c         s ] [ a ] = [ r ]
!>    [ -conjg(s)  c ] [ b ]   [ 0 ]
!> where c is real, s is complex, and c**2 + conjg(s)*s = 1.
!>
!> The computation uses the formulas
!>    |x| = sqrt( Re(x)**2 + Im(x)**2 )
!>    sgn(x) = x / |x|  if x /= 0
!>           = 1        if x  = 0
!>    c = |a| / sqrt(|a|**2 + |b|**2)
!>    s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
!>    r = sgn(a)*sqrt(|a|**2 + |b|**2)
!> When a and b are real and r /= 0, the formulas simplify to
!>    c = a / r
!>    s = b / r
!> the same as in DROTG when |a| > |b|.  When |b| >= |a|, the
!> sign of c and s will be different from those computed by DROTG
!> if the signs of a and b are not the same.
!>
!> \endverbatim
!>
!> @see lartg, @see lartgp
!
!  Arguments:
!  ==========
!
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE COMPLEX
!>          On entry, the scalar a.
!>          On exit, the scalar r.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE COMPLEX
!>          The scalar b.
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
!>          S is DOUBLE COMPLEX
!>          The scalar s.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Weslley Pereira, University of Colorado Denver, USA
!
!> \date December 2021
!
!> \ingroup rotg
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!> Based on the algorithm from
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!> \endverbatim
!
!  =====================================================================
subroutine ZROTG( a, b, c, s )
   integer, parameter :: wp = kind(1.d0)
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
   real(wp), parameter :: rtmin = sqrt( safmin )
!  ..
!  .. Scalar Arguments ..
   real(wp) :: c
   complex(wp) :: a, b, s
!  ..
!  .. Local Scalars ..
   real(wp) :: d, f1, f2, g1, g2, h2, u, v, w, rtmax
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
      if( real(g) == zero ) then
         r = abs(aimag(g))
         s = conjg( g ) / r
      elseif( aimag(g) == zero ) then
         r = abs(real(g))
         s = conjg( g ) / r
      else
         g1 = max( abs(real(g)), abs(aimag(g)) )
         rtmax = sqrt( safmax/2 )
         if( g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
!           The following two lines can be replaced by `d = abs( g )`.
!           This algorithm do not use the intrinsic complex abs.
            g2 = ABSSQ( g )
            d = sqrt( g2 )
            s = conjg( g ) / d
            r = d
         else
!
!        Use scaled algorithm
!
            u = min( safmax, max( safmin, g1 ) )
            gs = g / u
!           The following two lines can be replaced by `d = abs( gs )`.
!           This algorithm do not use the intrinsic complex abs.
            g2 = ABSSQ( gs )
            d = sqrt( g2 )
            s = conjg( gs ) / d
            r = d*u
         end if
      end if
   else
      f1 = max( abs(real(f)), abs(aimag(f)) )
      g1 = max( abs(real(g)), abs(aimag(g)) )
      rtmax = sqrt( safmax/4 )
      if( f1 > rtmin .and. f1 < rtmax .and. &
          g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         f2 = ABSSQ( f )
         g2 = ABSSQ( g )
         h2 = f2 + g2
         ! safmin <= f2 <= h2 <= safmax
         if( f2 >= h2 * safmin ) then
            ! safmin <= f2/h2 <= 1, and h2/f2 is finite
            c = sqrt( f2 / h2 )
            r = f / c
            rtmax = rtmax * 2
            if( f2 > rtmin .and. h2 < rtmax ) then
               ! safmin <= sqrt( f2*h2 ) <= safmax
               s = conjg( g ) * ( f / sqrt( f2*h2 ) )
            else
               s = conjg( g ) * ( r / h2 )
            end if
         else
            ! f2/h2 <= safmin may be subnormal, and h2/f2 may overflow.
            ! Moreover,
            !  safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= safmax,
            !  sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax).
            ! Also,
            !  g2 >> f2, which means that h2 = g2.
            d = sqrt( f2 * h2 )
            c = f2 / d
            if( c >= safmin ) then
               r = f / c
            else
               ! f2 / sqrt(f2 * h2) < safmin, then
               !  sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2) <= h2 * (safmin / f2) <= h2 <= safmax
               r = f * ( h2 / d )
            end if
            s = conjg( g ) * ( f / d )
         end if
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, f1, g1 ) )
         gs = g / u
         g2 = ABSSQ( gs )
         if( f1 / u < rtmin ) then
!
!           f is not well-scaled when scaled by g1.
!           Use a different scaling for f.
!
            v = min( safmax, max( safmin, f1 ) )
            w = v / u
            fs = f / v
            f2 = ABSSQ( fs )
            h2 = f2*w**2 + g2
         else
!
!           Otherwise use the same scaling for f and g.
!
            w = one
            fs = f / u
            f2 = ABSSQ( fs )
            h2 = f2 + g2
         end if
         ! safmin <= f2 <= h2 <= safmax
         if( f2 >= h2 * safmin ) then
            ! safmin <= f2/h2 <= 1, and h2/f2 is finite
            c = sqrt( f2 / h2 )
            r = fs / c
            rtmax = rtmax * 2
            if( f2 > rtmin .and. h2 < rtmax ) then
               ! safmin <= sqrt( f2*h2 ) <= safmax
               s = conjg( gs ) * ( fs / sqrt( f2*h2 ) )
            else
               s = conjg( gs ) * ( r / h2 )
            end if
         else
            ! f2/h2 <= safmin may be subnormal, and h2/f2 may overflow.
            ! Moreover,
            !  safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= safmax,
            !  sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax).
            ! Also,
            !  g2 >> f2, which means that h2 = g2.
            d = sqrt( f2 * h2 )
            c = f2 / d
            if( c >= safmin ) then
               r = fs / c
            else
               ! f2 / sqrt(f2 * h2) < safmin, then
               !  sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2) <= h2 * (safmin / f2) <= h2 <= safmax
               r = fs * ( h2 / d )
            end if
            s = conjg( gs ) * ( fs / d )
         end if
         ! Rescale c and r
         c = c * w
         r = r * u
      end if
   end if
   a = r
   return
end subroutine
