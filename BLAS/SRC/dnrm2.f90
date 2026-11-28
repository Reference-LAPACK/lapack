!> \brief \b DNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    DNRM2 := sqrt( x'*x )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER, storage spacing between elements of X
!>          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!>          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!>          If INCX = 0, x isn't a vector so there is no need to call
!>          this subroutine.  If you call it anyway, it will count x(1)
!>          in the vector norm N times.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup nrm2
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
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
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!>
!  =====================================================================
function DNRM2( n, x, incx )
   implicit none
   integer, parameter :: wp = kind(1.d0)
   real(wp) :: DNRM2
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  .. Constants ..
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one  = 1.0_wp
   real(wp), parameter :: maxN = huge(0.0_wp)
!  ..
!  .. Blue's scaling constants ..
   real(wp), parameter :: tsml = real(radix(0._wp), wp)**ceiling( &
       (minexponent(0._wp) - 1) * 0.5_wp)
   real(wp), parameter :: tbig = real(radix(0._wp), wp)**floor( &
       (maxexponent(0._wp) - digits(0._wp) + 1) * 0.5_wp)
   real(wp), parameter :: ssml = real(radix(0._wp), wp)**( - floor( &
       (minexponent(0._wp) - digits(0._wp)) * 0.5_wp))
   real(wp), parameter :: sbig = real(radix(0._wp), wp)**( - ceiling( &
       (maxexponent(0._wp) + digits(0._wp) - 1) * 0.5_wp))
!  ..
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*)
!  ..
!
!  A conservative exact threshold for the exceptional path.  Unlike
!  (maxN/scl)**2, this parameter cannot raise inexact on an ordinary call.
   real(wp), parameter :: overflow_bound = scale( one-epsilon(one), &
                          2*(maxexponent(zero)+exponent(sbig)-1) )
!
!  Exact-bin parameters for the exceptional overflow boundary only.
   integer, parameter :: bin_bits = max( 1, min( digits(zero)/3, &
                                      (maxexponent(zero)-1)/3 ) )
   integer, parameter :: bin_count = (digits(zero)+bin_bits-1)/bin_bits
   integer, parameter :: bin_lo = 2*(minexponent(zero)-digits(zero))/bin_bits &
                                 -2*bin_count-2
   real(wp), parameter :: bin_base = real(radix(zero),wp)**bin_bits
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin, entry
!
!  Quick return if possible
!
   DNRM2 = zero
   if( n <= 0 ) return
!
   scl = one
   sumsq = zero
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(x(ix))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         abig = abig + (amed*sbig)*sbig
      end if
      scl = one / sbig
      sumsq = abig
!
!     A rounded scaled sum can make the final product overflow even when
!     the exact norm fits.  Recompute only this exceptional boundary; the
!     ordinary Blue accumulation and all other finishing paths are unchanged.
!     First exclude sums too small to reach the boundary even at the
!     largest index value.  This exact kind-dependent scale normally folds
!     to a constant and avoids the integer logarithm on ordinary calls.
      entry = zero
      if (bit_size(n)-1 <= exponent(overflow_bound)-minexponent(zero)) &
         entry = scale( overflow_bound, 1-bit_size(n) )
      if (sumsq >= entry) then
         if (sumsq >= entry_threshold(n)) then
            DNRM2 = finite_nrm2( n, x, incx, scl, sumsq )
            return
         end if
      end if
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scl = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scl = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range
!
      scl = one
      sumsq = amed
   end if
   DNRM2 = scl*sqrt( sumsq )
   return
contains
!
!  Only the exceptional top-range path uses the following helpers.  The
!  accumulator digits are small exact integers stored in the working real
!  kind; no wider real or integer kind, IEEE module, or rounding-mode change
!  is needed.  With B = radix**b and 3*b <= digits, each digit product and
!  carry is exact in every rounding mode.
!
!  Compute ceil(log2(N)) in the existing integer kind.  N is positive,
!  so N-1 and every quotient remain in range, including the largest N.
!  If an unusual model cannot scale normally, zero is conservative and
!  avoids introducing an underflow/inexact exception in the guard itself.
   function entry_threshold( nvec ) result(bound)
      integer, intent(in) :: nvec
      integer :: left, shift
      real(wp) :: bound

      left = nvec-1
      shift = 0
      do while (left > 0)
         left = left/2
         shift = shift+1
      end do
      bound = zero
      if (shift <= exponent(overflow_bound)-minexponent(zero)) &
         bound = scale( overflow_bound, -shift )
   end function entry_threshold

   function finite_nrm2( nvec, xvec, stride, old_scale, old_sumsq ) result(ans)
      integer, intent(in) :: nvec, stride
      real(wp), intent(in) :: xvec(*), old_scale, old_sumsq
      real(wp) :: ans, ax, value, root, boundary, delta, ulp
      real(wp) :: accum(bin_lo:(2*maxexponent(zero)+bit_size(nvec))/bin_bits+2)
      real(wp) :: bound(bin_lo:(2*maxexponent(zero)+bit_size(nvec))/bin_bits+2)
      integer :: j, ix2, high, first, power, half, shift, relation

      accum = zero
      ix2 = 1
      if (stride < 0) ix2 = 1 - (nvec-1)*stride
      do j = 1, nvec
         ax = abs(xvec(ix2))
         if (ax > maxN .or. ax /= ax) then
!           Preserve the original propagation of infinity and NaN.
            ans = old_scale*sqrt(old_sumsq)
            return
         end if
         call add_product( ax, ax, accum )
         ix2 = ix2 + stride
      end do

      bound = zero
      call add_product( maxN, maxN, bound )
      relation = compare_bins( accum, bound )
      if (relation == 0) then
         ans = maxN
         return
      end if
      if (relation > 0) then
!        Classify the exact square against the nearest overflow midpoint:
!          4*S  versus  4*HUGE**2 + 4*HUGE*ulp + ulp**2.
!        This never forms an overflowing real product or an unrepresentable
!        midpoint.  Directed rounding does not overflow until the exact
!        norm reaches HUGE+ulp, so also compare against that square.  Choose
!        an increment from the matching interval; using more than one ulp
!        below that second boundary would raise a spurious overflow flag.
!        Data dependence prevents folding in a different rounding mode.
         call times_four( accum )
         call times_four( bound )
         ulp = scale( one, maxexponent(zero)-digits(zero) )
         do j = 1, 4
            call add_product( maxN, ulp, bound )
         end do
         call add_product( ulp, ulp, bound )
         relation = compare_bins( accum, bound )
!        Extend the bound to 4*(HUGE+ulp)**2 without forming HUGE+ulp.
         do j = 1, 4
            call add_product( maxN, ulp, bound )
         end do
         do j = 1, 3
            call add_product( ulp, ulp, bound )
         end do
         do high = ubound(accum,1), bin_lo, -1
            if (accum(high) /= zero) exit
         end do
         if (compare_bins(accum,bound) >= 0) then
            delta = scale( one+fraction(accum(high)), &
                           maxexponent(zero)-digits(zero) )
         else if (relation < 0 .or. &
                  (relation == 0 .and. mod(radix(zero),2) /= 0)) then
            delta = scale( fraction(accum(high)), &
                           maxexponent(zero)-digits(zero)-2 )
         else
!           Half an ulp < delta < one ulp: nearest/upward overflow, but
!           downward/toward-zero return HUGE without signaling overflow.
            delta = (ulp/2)*(one+fraction(accum(high)))
         end if
         ans = maxN + delta
         return
      end if

!     The exact norm fits.  Form a working-precision approximation from its
!     leading bins.  Two guard bins bound the omitted tail far below one
!     working ulp; a sticky upper endpoint prevents dropping it outright.
      do high = ubound(accum,1), bin_lo, -1
         if (accum(high) /= zero) exit
      end do
      if (high < bin_lo) then
         ans = zero
         return
      end if
      first = max( bin_lo, high-bin_count-1 )
      value = zero
      do j = bin_lo, first-1
         if (accum(j) /= zero) value = bin_base
      end do
      do j = first, high
         value = value/bin_base + accum(j)
      end do
      power = high*bin_bits
      half = power/2
      shift = mod( power, 2 )
      if (shift < 0) then
         half = half-1
         shift = shift+2
      end if
      root = sqrt( scale(value,shift) )
      boundary = scale( maxN, -half )
!     Avoid raising overflow before correcting a finite boundary result.
      ans = scale( min(root,boundary), half )
   end function finite_nrm2

   subroutine add_product( a, b, bins )
      real(wp), intent(in) :: a, b
      real(wp), intent(inout) :: bins(bin_lo:)
      real(wp) :: da(bin_count), db(bin_count), fa, fb, carry, digit, q, total
      integer :: j, k, index, power, shift, start

      if (a == zero .or. b == zero) return
      fa = fraction(a)
      fb = fraction(b)
      do j = 1, bin_count
         fa = fa*bin_base
         fb = fb*bin_base
         da(j) = aint(fa)
         db(j) = aint(fb)
         fa = fa-da(j)
         fb = fb-db(j)
      end do
      power = exponent(a)+exponent(b)-2*bin_count*bin_bits
      start = power/bin_bits
      shift = mod( power, bin_bits )
      if (shift < 0) then
         start = start-1
         shift = shift+bin_bits
      end if
      do j = 1, bin_count
         do k = 1, bin_count
            carry = scale( da(j)*db(k), shift )
            index = start+2*bin_count-j-k
            do while (carry > zero)
               q = aint(carry/bin_base)
               digit = carry-q*bin_base
               total = bins(index)+digit
               if (total >= bin_base) then
                  total = total-bin_base
                  q = q+one
               end if
               bins(index) = total
               carry = q
               index = index+1
            end do
         end do
      end do
   end subroutine add_product

   subroutine times_four( bins )
      real(wp), intent(inout) :: bins(bin_lo:)
      real(wp) :: total, carry
      integer :: j

      carry = zero
      do j = bin_lo, ubound(bins,1)
         total = 4*bins(j)+carry
         carry = aint(total/bin_base)
         bins(j) = total-carry*bin_base
      end do
   end subroutine times_four

   integer function compare_bins( lhs, rhs ) result(order)
      real(wp), intent(in) :: lhs(bin_lo:), rhs(bin_lo:)
      integer :: j

      order = 0
      do j = ubound(lhs,1), bin_lo, -1
         if (lhs(j) > rhs(j)) then
            order = 1
            return
         else if (lhs(j) < rhs(j)) then
            order = -1
            return
         end if
      end do
   end function compare_bins

end function
