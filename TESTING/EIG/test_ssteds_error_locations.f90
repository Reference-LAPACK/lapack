! Exercise diagnostic encoding only.  SSTEQR below is a deliberate mock:
! it supplies an orthogonal leaf basis, then forces a selected leaf to fail.
! SSTEDS and SLAES0 are compiled into this program by the build files, so
! SLAES0 calls this mock even when LAPACK is a shared library.
! No eigenvalue accuracy or performance is measured by this program.
program test_ssteds_error_locations
  implicit none
  integer :: n, i, mode, info, expected, nerr, ntest, split, lw, liw
  integer :: rows_done, fail_mode, fail_first, fail_last, wanted
  integer :: iq(1), isel(1)
  common /ssteds_failure_test/ rows_done, fail_mode, fail_first, &
                              fail_last, wanted
  real :: wq(1), dummy(1)
  real, allocatable :: d(:),e(:),w(:)
  integer, allocatable :: iw(:)
  external :: ssteds,slasum
  n=50000
  allocate(d(n),e(n))
  isel=0
  dummy=0e0
  call ssteds(n,d,e,0,isel,dummy,n,wq,-1,iq,-1,info)
  if(info/=0) then
    deallocate(d,e)
    stop 1
  endif
  lw=int(wq(1))
  liw=iq(1)
  allocate(w(lw),iw(liw))
  nerr=0
  ntest=0
  do mode=1,3
    rows_done=0
    fail_mode=mode
    fail_first=0
    fail_last=0
    wanted=n
    d=2e0
    e=0.1e0
    ! Case 3 fails late in the second block: local encoding fits INT32,
    ! while the final global row/column encoding does not.
    if(mode==3) e(30000)=0e0
    call ssteds(n,d,e,0,isel,dummy,n,w,lw,iw,liw,info)
    ntest=ntest+1
    if(fail_first==0 .or. fail_last<fail_first) then
      nerr=nerr+1
      print *, 'Mock did not reach the requested failure: ',mode
    endif
    if(fail_first>(huge(0)-fail_last)/(n+1)) then
      expected=1
    else
      expected=fail_first*(n+1)+fail_last
    endif
    ntest=ntest+1
    if(info/=expected) then
      nerr=nerr+1
      print *, 'Wrong failure encoding: ',mode,info,expected
    endif
  enddo
  call slasum('SST',6,nerr,ntest)
  deallocate(d,e,w,iw)
  if(nerr/=0) stop 1
end program test_ssteds_error_locations

subroutine ssteqr(compz,n,d,e,z,ldz,work,info)
  implicit none
  character :: compz
  integer :: n,ldz,info,i,j
  real :: d(*),e(*),z(ldz,*),work(*)
  integer :: rows_done, fail_mode, fail_first, fail_last, wanted
  common /ssteds_failure_test/ rows_done, fail_mode, fail_first, &
                              fail_last, wanted
  rows_done=rows_done+n
  if((fail_mode==1 .and. rows_done==n) .or. &
      (fail_mode/=1 .and. rows_done==wanted)) then
    fail_first=rows_done-n+1
    fail_last=rows_done
    info=1
    return
  endif
  info=0
  ! Sorted poles and an identity basis keep subsequent rank-one merges
  ! algebraically valid.  They are not the eigenpairs of the input leaf.
  do i=1,n
    d(i)=real(i)/real(n+1)
    do j=1,n
      z(i,j)=0e0
    enddo
    z(i,i)=1e0
  enddo
end subroutine ssteqr
