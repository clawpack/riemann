!   =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
!   =====================================================

!     # Riemann solver in the transverse direction for the acoustics equations.

!     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.
    implicit none
    !Input
    integer, intent(in)  :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1,aux2,aux3
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: asdq

    !Output
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: bmasdq,bpasdq

    !Local
    integer :: i,mu,mv
    double precision :: a1,a2
    double precision :: rho,bulk,cc,zz
    

!     # density, bulk modulus, and sound speed, and impedence of medium:
!     # (should be set in setprob.f)
    common /cparam/ rho,bulk,cc,zz


    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do 20 i = 2-mbc, mx+mbc
        a1 = (asdq(1,i)*zz - asdq(mv,i)) / (2.d0*zz)
        a2 = (zz*asdq(1,i) + asdq(mv,i)) / (2.d0*zz)
    
    !        # The down-going flux difference bmasdq is the product  -c * wave
    
        bmasdq(1,i) = -cc * a2
        bmasdq(mu,i) = 0.d0
        bmasdq(mv,i) = -cc * a2 * zz
    
    !        # The up-going flux difference bpasdq is the product  c * wave
    
        bpasdq(1,i) = cc * a1
        bpasdq(mu,i) = 0.d0
        bpasdq(mv,i) = -cc * a1 * zz
    
    20 END DO

    return
    end subroutine rpt2
