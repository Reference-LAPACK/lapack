! Check RANGE='V' endpoint selection for the bidiagonal and dense SVDs
! (issue #1492).  The documented interval is [VL,VU).  Dense drivers
! also reach bidiagonal order 1 for one-row and one-column matrices.
program test_dsvd_range_boundaries
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer, parameter :: lwork = 1000
  real(dp) :: diagonal(2), d(2), e(2), s(2), z(4,3)
  real(dp) :: a(3,3), u(3,3), vt(3,3), work(lwork)
  real(dp) :: vl, vu, expected_s, tol
  integer :: iwork(lwork), icase, ishape, iroutine, m, n, k, ns, info
  integer :: expected_ns, failures, first_routine
  character(10) :: routine
  external :: dbdsvdx, dbdsvdmr3
  external :: dgesvdx, dgesvdmr3

  failures = 0
  do icase = 1, 3
    select case (icase)
    case (1)
      ! The lower endpoint is included.
      diagonal = [2.0_dp, 4.0_dp]
      vl = 2.0_dp
      vu = 3.0_dp
      expected_ns = 1
      expected_s = 2.0_dp
    case (2)
      ! The upper endpoint is excluded.
      diagonal = [2.0_dp, 4.0_dp]
      vl = 1.0_dp
      vu = 2.0_dp
      expected_ns = 0
      expected_s = 0.0_dp
    case (3)
      ! Zero is included when VL is zero.
      diagonal = [0.0_dp, 2.0_dp]
      vl = 0.0_dp
      vu = 1.0_dp
      expected_ns = 1
      expected_s = 0.0_dp
    end select

    do ishape = 1, 4
      select case (ishape)
      case (1)
        m = 1; n = 1
      case (2)
        m = 2; n = 2
      case (3)
        m = 3; n = 1
      case (4)
        m = 1; n = 3
      end select
      k = min(m,n)
      first_routine = 1
      if (ishape >= 3) first_routine = 3
      do iroutine = first_routine, 4
        d = diagonal
        e = 0.0_dp
        a = 0.0_dp
        a(1,1) = diagonal(1)
        if (ishape == 2) a(2,2) = diagonal(2)
        s = -1.0_dp
        ns = -1
        info = -1

        select case (iroutine)
        case (1)
          routine = 'DBDSVDX'
          call dbdsvdx('U','N','V',k,d,e,vl,vu,1,k,ns,s,z,4, &
                       work,iwork,info)
        case (2)
          routine = 'DBDSVDMR3'
          call dbdsvdmr3('U','N','V',k,d,e,vl,vu,1,k,ns,s,z,4, &
                         work,lwork,iwork,lwork,info)
        case (3)
          routine = 'DGESVDX'
          call dgesvdx('N','N','V',m,n,a,3,vl,vu,1,k,ns,s,u,3, &
                       vt,3,work,lwork,iwork,info)
        case (4)
          routine = 'DGESVDMR3'
          call dgesvdmr3('N','N','V',m,n,a,3,vl,vu,1,k,ns,s,u,3, &
                         vt,3,work,lwork,iwork,info)
        end select

        if (info /= 0 .or. ns /= expected_ns) then
          print *, trim(routine), ' case=', icase, ' m=', m, ' n=', n, &
                   ' INFO=', info, ' NS=', ns, ' expected NS=', expected_ns
          failures = failures + 1
        else if (ns > 0) then
          tol = 10.0_dp*epsilon(1.0_dp)* &
                max(1.0_dp, expected_s)
          if (.not.(abs(s(1)-expected_s) <= tol)) then
            print *, trim(routine), ' case=', icase, ' m=', m, ' n=', n, &
                     ' S(1)=', s(1), ' expected=', expected_s
            failures = failures + 1
          end if
        end if
      end do
    end do
  end do

  call check_singleton_vectors(failures)
  print '(a,i0)', 'DSVD range boundary failures = ', failures
  if (failures /= 0) stop 1

contains

  subroutine check_singleton_vectors(failures)
    integer, intent(inout) :: failures
    real(dp) :: d(1), e(1), s(1), z(2,2)
    real(dp) :: w(lwork), tol
    integer :: iw(lwork), iroutine, ns, info

    tol = 100.0_dp*epsilon(1.0_dp)
    do iroutine = 1, 2
      d = 2.0_dp
      e = 0.0_dp
      z = 0.0_dp
      select case (iroutine)
      case (1)
        call dbdsvdx('U','V','V',1,d,e,2.0_dp,3.0_dp,1,1, &
                      ns,s,z,2,w,iw,info)
      case (2)
        call dbdsvdmr3('U','V','V',1,d,e,2.0_dp,3.0_dp,1,1, &
                        ns,s,z,2,w,lwork,iw,lwork,info)
      end select
      if (info /= 0 .or. ns /= 1) then
        print *, 'Singleton vectors routine=', iroutine, &
                 ' INFO=', info, ' NS=', ns
        failures = failures + 1
        cycle
      end if
      if (abs(s(1)-2.0_dp) > tol) then
        print *, 'Singleton singular value routine=', iroutine
        failures = failures + 1
      end if
      if (abs(z(1,1)*z(2,1)-1.0_dp) > tol) failures = failures + 1
    end do
  end subroutine
end program
