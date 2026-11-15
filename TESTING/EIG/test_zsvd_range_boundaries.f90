! Check RANGE='V' endpoint selection for the complex*16 dense SVDs
! (issue #1492).  The documented interval is [VL,VU).
program test_zsvd_range_boundaries
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer, parameter :: lwork = 1000
  real(dp), parameter :: canary = 12345.0_dp
  complex(dp) :: a(2,2), u(2,2), vt(2,2), work(lwork), original
  real(dp) :: s_guard(4), rwork(lwork), vl, vu, tol
  integer :: iwork(lwork), iroutine, icase, n, ns, info, expected_ns
  integer :: failures
  character(10) :: routine
  character :: jobz
  external :: zgesvdx, zgesvdmr3

  failures = 0
  tol = 100.0_dp*epsilon(1.0_dp)
  do iroutine = 1, 2
    do icase = 1, 4
      a = cmplx(0.0_dp,0.0_dp,dp)
      u = cmplx(0.0_dp,0.0_dp,dp)
      vt = cmplx(0.0_dp,0.0_dp,dp)
      s_guard = canary
      n = 1
      expected_ns = 1
      vl = 2.0_dp
      vu = 3.0_dp
      jobz = 'N'
      original = cmplx(0.0_dp,2.0_dp,dp)
      a(1,1) = original

      select case (icase)
      case (2)
        ! The upper endpoint is excluded.
        vl = 1.0_dp
        vu = 2.0_dp
        expected_ns = 0
      case (3)
        ! Zero is included when VL is zero.
        n = 2
        vl = 0.0_dp
        vu = 1.0_dp
        expected_ns = 1
        a(1,1) = cmplx(0.0_dp,0.0_dp,dp)
        a(2,2) = cmplx(2.0_dp,0.0_dp,dp)
      case (4)
        ! The 1x1 complex singular vectors reconstruct the input.
        jobz = 'V'
      end select

      if (iroutine == 1) then
        routine = 'ZGESVDX'
        call zgesvdx(jobz,jobz,'V',n,n,a,2,vl,vu,1,n,ns, &
                     s_guard(2),u,2,vt,2,work,lwork,rwork,iwork,info)
      else
        routine = 'ZGESVDMR3'
        call zgesvdmr3(jobz,jobz,'V',n,n,a,2,vl,vu,1,n,ns, &
                       s_guard(2),u,2,vt,2,work,lwork,rwork,iwork,info)
      end if

      if (info /= 0 .or. ns /= expected_ns) then
        print *, trim(routine), ' case=', icase, ' INFO=', info, &
                 ' NS=', ns, ' expected NS=', expected_ns
        failures = failures + 1
      end if
      if (s_guard(1) /= canary .or. any(s_guard(n+2:4) /= canary)) then
        print *, trim(routine), ' wrote beyond S, case=', icase
        failures = failures + 1
      end if
      if (info /= 0 .or. ns /= expected_ns) cycle
      if (icase == 3) then
        if (abs(s_guard(2)) > tol) then
          print *, trim(routine), ' wrong zero-range values'
          failures = failures + 1
        end if
      else if (ns == 1) then
        if (abs(s_guard(2)-2.0_dp) > tol) then
          print *, trim(routine), ' wrong endpoint value, case=', icase
          failures = failures + 1
        end if
      end if
      if (icase == 4 .and. ns == 1) then
        if (abs(u(1,1)*cmplx(s_guard(2),0.0_dp,dp)*vt(1,1) - &
                original) > tol .or. abs(abs(u(1,1))-1.0_dp) > tol .or. &
            abs(abs(vt(1,1))-1.0_dp) > tol) then
          print *, trim(routine), ' wrong 1x1 vectors'
          failures = failures + 1
        end if
      end if
    end do
  end do

  print '(a,i0)', 'ZSVD range boundary failures = ', failures
  if (failures /= 0) stop 1
end program
