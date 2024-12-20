! =============================================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =============================================================================
!
! Riemann problems for the 1D Burgers' equation with entropy fix for 
! transonic rarefaction. See "Finite Volume Method for Hyperbolic Problems",
! R. J. LeVeque.

! waves: 1
! equations: 1

! Conserved quantities:
!       1 q
    
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
    logical :: efix
 
    efix = .true.   !# Compute correct flux for transonic rarefactions
 

    do i=2-mbc,mx+mbc
        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        s(1,i) = 0.5d0 * (qr(1,i-1) + ql(1,i))

        amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
        apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)

        if (efix) then
            if (ql(1,i).gt.0.d0 .and. qr(1,i-1).lt.0.d0) then
                amdq(1,i) = - 1.d0/2.d0 * qr(1,i-1)**2
                apdq(1,i) =   1.d0/2.d0 * ql(1,i)**2
            endif
        endif
    enddo

    return
    end subroutine
