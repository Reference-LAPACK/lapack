! Check RANGE='V' endpoint selection for the bidiagonal and dense SVDs
! (issue #1492).  The documented interval is [VL,VU).  Dense drivers
! also reach bidiagonal order 1 for one-row and one-column matrices.
program test_ssvd_range_boundaries
  use, intrinsic :: iso_fortran_env, only: sp => real32
  implicit none
  integer, parameter :: lwork = 1000
  real(sp) :: diagonal(2), d(2), e(2), s(2), z(4,3)
  real(sp) :: a(3,3), u(3,3), vt(3,3), work(lwork)
  real(sp) :: vl, vu, expected_s, tol
  integer :: iwork(lwork), icase, ishape, iroutine, m, n, k, ns, info
  integer :: expected_ns, failures, first_routine
  character(10) :: routine
  external :: sbdsvdx, sbdsvdmr3
  external :: sgesvdx, sgesvdmr3

  failures = 0
  do icase = 1, 3
    select case (icase)
    case (1)
      ! The lower endpoint is included.
      diagonal = [2.0_sp, 4.0_sp]
      vl = 2.0_sp
      vu = 3.0_sp
      expected_ns = 1
      expected_s = 2.0_sp
    case (2)
      ! The upper endpoint is excluded.
      diagonal = [2.0_sp, 4.0_sp]
      vl = 1.0_sp
      vu = 2.0_sp
      expected_ns = 0
      expected_s = 0.0_sp
    case (3)
      ! Zero is included when VL is zero.
      diagonal = [0.0_sp, 2.0_sp]
      vl = 0.0_sp
      vu = 1.0_sp
      expected_ns = 1
      expected_s = 0.0_sp
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
        e = 0.0_sp
        a = 0.0_sp
        a(1,1) = diagonal(1)
        if (ishape == 2) a(2,2) = diagonal(2)
        s = -1.0_sp
        ns = -1
        info = -1

        select case (iroutine)
        case (1)
          routine = 'SBDSVDX'
          call sbdsvdx('U','N','V',k,d,e,vl,vu,1,k,ns,s,z,4, &
                       work,iwork,info)
        case (2)
          routine = 'SBDSVDMR3'
          call sbdsvdmr3('U','N','V',k,d,e,vl,vu,1,k,ns,s,z,4, &
                         work,lwork,iwork,lwork,info)
        case (3)
          routine = 'SGESVDX'
          call sgesvdx('N','N','V',m,n,a,3,vl,vu,1,k,ns,s,u,3, &
                       vt,3,work,lwork,iwork,info)
        case (4)
          routine = 'SGESVDMR3'
          call sgesvdmr3('N','N','V',m,n,a,3,vl,vu,1,k,ns,s,u,3, &
                         vt,3,work,lwork,iwork,info)
        end select

        if (info /= 0 .or. ns /= expected_ns) then
          print *, trim(routine), ' case=', icase, ' m=', m, ' n=', n, &
                   ' INFO=', info, ' NS=', ns, ' expected NS=', expected_ns
          failures = failures + 1
        else if (ns > 0) then
          tol = 10.0_sp*epsilon(1.0_sp)* &
                max(1.0_sp, expected_s)
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
  print '(a,i0)', 'SSVD range boundary failures = ', failures
  if (failures /= 0) stop 1

contains

  subroutine check_singleton_vectors(failures)
    integer, intent(inout) :: failures
    real(sp) :: d(1), e(1), s(1), z(2,2)
    real(sp) :: w(lwork), tol
    integer :: iw(lwork), iroutine, ns, info

    tol = 100.0_sp*epsilon(1.0_sp)
    do iroutine = 1, 2
      d = 2.0_sp
      e = 0.0_sp
      z = 0.0_sp
      select case (iroutine)
      case (1)
        call sbdsvdx('U','V','V',1,d,e,2.0_sp,3.0_sp,1,1, &
                      ns,s,z,2,w,iw,info)
      case (2)
        call sbdsvdmr3('U','V','V',1,d,e,2.0_sp,3.0_sp,1,1, &
                        ns,s,z,2,w,lwork,iw,lwork,info)
      end select
      if (info /= 0 .or. ns /= 1) then
        print *, 'Singleton vectors routine=', iroutine, &
                 ' INFO=', info, ' NS=', ns
        failures = failures + 1
        cycle
      end if
      if (abs(s(1)-2.0_sp) > tol) then
        print *, 'Singleton singular value routine=', iroutine
        failures = failures + 1
      end if
      if (abs(z(1,1)*z(2,1)-1.0_sp) > tol) failures = failures + 1
    end do
  end subroutine
end program
