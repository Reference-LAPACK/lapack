! Exercise the INFO encoding of a failing merge (secular equation) in
! SSTEDS.  SSTEQR and SLAED4 below are deliberate mocks: SSTEQR returns
! an orthogonal leaf basis, so that the first merge has two
! nondeflated poles, and SLAED4 then fails.  SSTEDS, SLAES0, SLAES1
! and SLAED9_BR are compiled into this program by the build files, so
! that both mocks are reached even when LAPACK is a shared library.
! No eigenvalue accuracy or performance is measured by this program.
program test_ssteds_merge_failures
  implicit none
  integer :: n, i, mode, info, expected, nerr, ntest, lw, liw
  integer :: rows_done, block_start, leaf1, leaf2, nsolve
  integer :: iq(1), isel(1)
  common /ssteds_merge_test/ rows_done, block_start, leaf1, leaf2, &
                             nsolve
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
  do mode=1,2
    rows_done=0
    leaf1=0
    leaf2=0
    nsolve=0
    d=2e0
    e=0.1e0
    ! Mode 1: the first merge of the matrix fails.
    ! Mode 2: rows 1:45000 split into single-leaf blocks, so the first
    ! merge, rows 45001 onward, fails where its global row/column
    ! encoding does not fit INT32.
    block_start=1
    if(mode==2) then
      block_start=45001
      do i=25,45000,25
        e(i)=0e0
      enddo
    endif
    call ssteds(n,d,e,0,isel,dummy,n,w,lw,iw,liw,info)
    ntest=ntest+1
    if(nsolve/=1 .or. leaf1==0 .or. leaf2==0) then
      nerr=nerr+1
      print *, 'Mock did not reach the requested merge: ',mode,nsolve
    endif
    if(block_start>(huge(0)-(block_start+leaf1+leaf2-1))/(n+1)) then
      expected=1
    else
      expected=block_start*(n+1)+block_start+leaf1+leaf2-1
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
end program test_ssteds_merge_failures

subroutine ssteqr(compz,n,d,e,z,ldz,work,info)
  implicit none
  character :: compz
  integer :: n,ldz,info,i,j
  real :: d(*),e(*),z(ldz,*),work(*)
  integer :: rows_done, block_start, leaf1, leaf2, nsolve
  common /ssteds_merge_test/ rows_done, block_start, leaf1, leaf2, &
                             nsolve
  ! Record the sizes of the first two leaves of the failing block.
  if(rows_done>=block_start-1) then
    if(leaf1==0) then
      leaf1=n
    else if(leaf2==0) then
      leaf2=n
    endif
  endif
  rows_done=rows_done+n
  info=0
  ! Sorted poles and an identity basis: the rank-one update of the
  ! first merge then couples only one pole of each child.
  do i=1,n
    d(i)=real(i)/real(n+1)
    do j=1,n
      z(i,j)=0e0
    enddo
    z(i,i)=1e0
  enddo
end subroutine ssteqr

subroutine slaed4(n,i,d,z,delta,rho,dlam,info)
  implicit none
  integer :: n,i,info
  real :: d(*),z(*),delta(*),rho,dlam
  integer :: rows_done, block_start, leaf1, leaf2, nsolve
  common /ssteds_merge_test/ rows_done, block_start, leaf1, leaf2, &
                             nsolve
  nsolve=nsolve+1
  info=1
end subroutine slaed4
