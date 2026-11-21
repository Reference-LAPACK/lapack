! Regression for finite top-range norms (issue #1436).
! Check values and IEEE exceptions without a wider working real kind.
program test_snrm2_overflow
   use, intrinsic :: ieee_arithmetic
   use, intrinsic :: iso_fortran_env, only: int32
   implicit none
   integer, parameter :: wp = kind(1.e0)
   real(wp), external :: snrm2
   real(wp) :: x(129), spaced(129), answer, small, h, u
   integer :: i, checks

   checks = 0
   x = 0.0_wp
   x(1:2) = (/3.0_wp,4.0_wp/)
   call require( snrm2(2,x,1) == 5.0_wp, 'ordinary norm' )
   call require( snrm2(0,x,1) == 0.0_wp, 'zero length' )
   call require( snrm2(-1,x,1) == 0.0_wp, 'negative length' )
   h = huge(0.0_wp)
   u = spacing(h)
   x(1) = h
   call require( snrm2(1,x,1) == h, 'scalar HUGE' )

   if (radix(h) == 2 .and. digits(h) == 24 .and. &
       maxexponent(h) == 128) then
!     The first square is just below the overflow boundary.  Each of the
!     seven smaller squares is just above half an ulp of the scaled sum,
!     so separate working-precision additions round upwards repeatedly.
!     The exact norm of this vector is nevertheless below HUGE.
      small = scale( 11863284.0_wp, maxexponent(h)-36 )
      x(1) = h-2*u
      x(2:8) = small
      call finite_boundary( snrm2(8,x,1), 'constructed forward' )
      call finite_boundary( snrm2(8,x,-1), 'constructed reverse' )
      spaced = 0.0_wp
      do i = 1, 8
         spaced(2*i-1) = x(i)
      end do
      call finite_boundary( snrm2(8,spaced,2), 'constructed stride two' )
      x(1:64) = h/8
      call require( snrm2(64,x,1) == h, 'exact HUGE with 64 entries' )
      call ieee_boundary_tests()
   end if

   x(1:2) = h
   answer = snrm2(2,x,1)
   call require( answer > h, 'true overflow must not be clamped' )
   write(*,'(a,i0)') 'NRM2 top-range checks passed: ', checks

contains
   subroutine require( condition, label )
      logical, intent(in) :: condition
      character(*), intent(in) :: label
      if (.not.condition) then
         write(*,'(a,a)') 'NRM2 top-range failure: ', label
         stop 1
      end if
      checks = checks+1
   end subroutine require

   subroutine finite_boundary( result, label )
      real(wp), intent(in) :: result
      character(*), intent(in) :: label
      call require( result <= h .and. result >= h-4*u, label )
   end subroutine finite_boundary
   subroutine ieee_boundary_tests()
      type(ieee_round_type) :: original, selected
      real(wp) :: values(129)
      integer :: mode, k, nvec, j

      if (.not.ieee_support_datatype(h) .or. &
          .not.ieee_support_flag(ieee_overflow,h) .or. &
          .not.ieee_support_flag(ieee_inexact,h) .or. &
          .not.ieee_support_flag(ieee_invalid,h) .or. &
          .not.ieee_support_flag(ieee_underflow,h) .or. &
          .not.ieee_support_flag(ieee_divide_by_zero,h)) return
      do mode = 1, 4
         selected = ieee_nearest
         if (mode == 2) selected = ieee_up
         if (mode == 3) selected = ieee_down
         if (mode == 4) selected = ieee_to_zero
         if (.not.ieee_support_rounding(selected,h)) return
      end do
      call ieee_get_rounding_mode(original)
      do mode = 1, 4
         selected = ieee_nearest
         if (mode == 2) selected = ieee_up
         if (mode == 3) selected = ieee_down
         if (mode == 4) selected = ieee_to_zero
         call ieee_set_rounding_mode(selected)
!        A positive increment below the nearest midpoint overflows only up.
         values(1) = h
         values(2) = u
         call flag_case( 2, values, mode == 2, mode == 2, 'low half', mode )
!        Exact norm HUGE+ulp/2: ties go up, but down/zero must not overflow.
         values(1:4) = h/2
         k = (digits(h)+2)/2
         nvec = 5
         values(5) = (scale(1.0_wp,k)-1.0_wp)*(u/2)
         if (mod(digits(h)+2,2) /= 0) then
            values(5) = scale(u,k-1)
            values(6) = (scale(1.0_wp,k)-1.0_wp)*(u/2)
            nvec = 6
         end if
         do j = 0, (k-1)/2-1
            values(nvec+1:nvec+3) = scale(u,j)
            nvec = nvec+3
         end do
         call flag_case( nvec, values, mode <= 2, mode <= 2, 'midpoint', mode )
!        Exact norm HUGE+ulp: every direction must signal genuine overflow.
         values(1:4) = scale(1.0_wp,maxexponent(h)-1)
         call flag_case( 4, values, mode <= 2, .true., 'next power', mode )
         values(1:64) = h/8
         call flag_case( 64, values, .false., .false., 'exact HUGE', mode )
      end do
      call ieee_set_rounding_mode(original)
   end subroutine ieee_boundary_tests

   subroutine flag_case( nvec, values, inf_expected, overflow_expected, label, mode )
      integer, intent(in) :: nvec, mode
      real(wp), intent(in) :: values(*)
      logical, intent(in) :: inf_expected, overflow_expected
      character(*), intent(in) :: label
      logical(kind=int32) :: flags(5), flag_expected
      logical :: good
      real(wp) :: result

!     GNU's IEEE module uses logical kind4 even with default index64.
      call ieee_set_flag(ieee_all,.false._int32)
      result = snrm2(nvec,values,1)
      call ieee_get_flag(ieee_overflow,flags(1))
      call ieee_get_flag(ieee_invalid,flags(2))
      call ieee_get_flag(ieee_inexact,flags(3))
      call ieee_get_flag(ieee_underflow,flags(4))
      call ieee_get_flag(ieee_divide_by_zero,flags(5))
      flag_expected = overflow_expected
      good = (flags(1) .eqv. flag_expected) .and. .not.flags(2) .and. &
             flags(3) .and. .not.flags(4) .and. .not.flags(5)
      if (inf_expected) then
         if (.not.(result > h)) good = .false.
      else
         if (.not.(result == h)) good = .false.
      end if
      if (.not.good) write(*,'(a,i0,5l2)') 'IEEE mode/flags: ', mode, flags
      call require( good, 'IEEE '//label )
   end subroutine flag_case
end program test_snrm2_overflow
