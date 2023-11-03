! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! Solve Riemann problems for the 1D advection equation q_t + u(x)*q_x = 0
! with variable u(x) in non-conservative form (the color equation)

! waves:     1
! equations: 1

! Conserved quantities:
!       1 q

! Auxiliary variables:
!       1 velocity

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q

! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routine step1, rp is called with ql = qr = q.


    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, maux, mbc, mx
    real(kind=8), intent(in) :: ql(meqn,1-mbc:maxmx+mbc)
    real(kind=8), intent(in) :: qr(meqn,1-mbc:maxmx+mbc)
    real(kind=8), intent(in) :: auxl(maux,1-mbc:maxmx+mbc)
    real(kind=8), intent(in) :: auxr(maux,1-mbc:maxmx+mbc)

    real(kind=8), intent(out) :: s(mwaves,1-mbc:maxmx+mbc)
    real(kind=8), intent(out) :: wave(meqn, mwaves,1-mbc:maxmx+mbc)
    real(kind=8), intent(out) :: amdq(meqn,1-mbc:maxmx+mbc)
    real(kind=8), intent(out) :: apdq(meqn,1-mbc:maxmx+mbc)

    real(kind=8) :: u
    integer :: i
    common /comrp/ u



    do 30 i=2-mbc,mx+mbc
    
    !        # Compute the wave and speed
    
        u = dmin1(auxr(1,i-1),0.d0) +  dmax1(auxl(1,i),0.d0)
    
    !        # the formula below seems to work as well...
    !        u = 0.5d0 * (auxr(1,i-1) +  auxl(1,i))

        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        s(1,i) = u
        amdq(1,i) = dmin1(auxr(1,i-1), 0.d0) * wave(1,1,i)
        apdq(1,i) = dmax1(auxl(1,i),   0.d0) * wave(1,1,i)
    30 END DO

    return
    end subroutine rp1
