! =============================================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =============================================================================
!
! Riemann problems for the 1D cubic equation q_t + (q^3)_x = 0.
!
! waves:     1
! equations: 1
!
! Conserved quantities:
!       1 q
!
! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
!
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q
!            (where the flux difference is f(qr(i-1)) - f(ql(i)))
!
! Note that the i'th Riemann problem has left state qr(:, i-1)
!                                    and right state ql(:, i)
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

    integer :: i


    do i = 2-mbc, mx+mbc

        wave(1, 1, i) = ql(1, i) - qr(1, i-1)
        ! (ul^3 - ur^3) / (ul - ur) == ul^2 + ul*ur + ur^2
        s(1, i) = ql(1, i)*ql(1, i) + ql(1, i)*qr(1, i-1) + qr(1, i-1)*qr(1, i-1)

        amdq(1,i) = 0.d0
        apdq(1,i) = ql(1, i)*ql(1, i)*ql(1, i) - qr(1, i-1)*qr(1, i-1)*qr(1, i-1)

    enddo

    return
end
