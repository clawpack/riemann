! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================

!     # Riemann solver in the transverse direction for the advection equation.

    implicit none
    !Input
    integer, intent(in)  :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1,aux2,aux3
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: asdq

    !Output
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: bmasdq,bpasdq

    !Local
    integer :: i,i1,kv

    kv = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1
    do 10 i=2-mbc,mx+mbc
        i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
        bmasdq(1,i) = dmin1(aux2(kv,i1), 0.d0) * asdq(1,i)
        bpasdq(1,i) = dmax1(aux3(kv,i1), 0.d0) * asdq(1,i)
    10 END DO

    return
    end subroutine rpt2
