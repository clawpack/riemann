! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! solve Riemann problems for the 1D advection equation q_t + u*q_x = 0.
! For constant advection velocity u, passed in common block.

! waves:     1
! equations: 1

! Conserved quantities:
!       1 q

! The advection speed u is passed in the common block cparam
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
    !Input
    integer, intent(in) :: maxmx, meqn, mwaves, maux, mbc, mx
    double precision, intent(in) ::   ql( meqn,1-mbc:maxmx+mbc)
    double precision, intent(in) ::   qr( meqn,1-mbc:maxmx+mbc)
    double precision, intent(in) :: auxl(maux,1-mbc:maxmx+mbc)
    double precision, intent(in) :: auxr(maux,1-mbc:maxmx+mbc)
    
    double precision, intent(out) ::   amdq( meqn,1-mbc:maxmx+mbc)
    double precision, intent(out) ::   apdq( meqn,1-mbc:maxmx+mbc)
    double precision, intent(out) ::    s(mwaves,1-mbc:maxmx+mbc)
    double precision, intent(out) :: wave(meqn, mwaves,1-mbc:maxmx+mbc)

    !Local
    double precision :: u
    integer :: i
    common /cparam/ u

    do i=2-mbc,mx+mbc
    
    !        # Compute the wave and speed
!       # (speed negative because we are going backward in time)
    
        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        s(1,i) = -u
        amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
        apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)
    END DO

    return
    end subroutine rp1
