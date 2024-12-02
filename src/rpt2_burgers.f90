!   =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
!   =====================================================

!     # Riemann solver in the transverse direction for 2D Burgers' equation
!
!     # Split asdq into eigenvectors of Roe matrix B.
!     # For the scalar equation, this simply amounts to computing the
!     # transverse wave speed from the opposite Riemann problem.
        implicit none
        !Input
        integer, intent(in)  :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx
        double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
        double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1,aux2,aux3
        double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: asdq

        !Output
        double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: bmasdq,bpasdq

        !Local
        integer :: i
        double precision :: sb
        
!     # x- and y- Riemann problems are identical, so it doesn't matter if
!     # ixy=1 or 2.

        do 10 i = 2-mbc, mx+mbc
            sb = 0.5d0*(qr(1,i-1) + ql(1,i))
            bmasdq(1,i) = dmin1(sb, 0.d0) * asdq(1,i)
            bpasdq(1,i) = dmax1(sb, 0.d0) * asdq(1,i)
        10 continue
        
    return
    end subroutine rpt2