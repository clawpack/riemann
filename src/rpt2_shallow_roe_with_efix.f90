! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================

!     # Riemann solver in the transverse direction for the shallow water
!     equations .
!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)

    implicit none
    !Input
    integer, intent(in)  :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1,aux2,aux3
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: asdq

    !Output
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: bmasdq,bpasdq

    !Local
    integer :: i,m,mw,mu,mv
    double precision :: hsqrtl,hsqrtr,hsq2,a1,a2,a3
    integer, parameter :: maxm2 = 1800
    !   # Roe averages quantities of each interface
    double precision, dimension(-6:maxm2+7) :: u,v,a,h
    double precision :: waveb(3,3),sb(3)
    double precision :: grav

!   # grav must be set elsewhere
    common /cparam/ grav


    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif
            
    do 10 i = 2-mbc, mx+mbc
        h(i) = (qr(1,i-1)+ql(1,i))*0.50d0
        hsqrtl = dsqrt(qr(1,i-1))
        hsqrtr = dsqrt(ql(1,i))
        hsq2 = hsqrtl + hsqrtr
        u(i) = (qr(mu,i-1)/hsqrtl + ql(mu,i)/hsqrtr) / hsq2
        v(i) = (qr(mv,i-1)/hsqrtl + ql(mv,i)/hsqrtr) / hsq2
        a(i) = dsqrt(grav*h(i))
    10 END DO


    do 20 i = 2-mbc, mx+mbc
        a1 = (0.50d0/a(i))*((v(i)+a(i))*asdq(1,i)-asdq(mv,i))
        a2 = asdq(mu,i) - u(i)*asdq(1,i)
        a3 = (0.50d0/a(i))*(-(v(i)-a(i))*asdq(1,i)+asdq(mv,i))
    
        waveb(1,1) = a1
        waveb(mu,1) = a1*u(i)
        waveb(mv,1) = a1*(v(i)-a(i))
        sb(1) = v(i) - a(i)
    
        waveb(1,2) = 0.0d0
        waveb(mu,2) = a2
        waveb(mv,2) = 0.0d0
        sb(2) = v(i)
    
        waveb(1,3) = a3
        waveb(mu,3) = a3*u(i)
        waveb(mv,3) = a3*(v(i)+a(i))
        sb(3) = v(i) + a(i)
    
!       # compute the flux differences bmasdq and bpasdq
    
        do m=1,meqn
            bmasdq(m,i) = 0.d0
            bpasdq(m,i) = 0.d0
            do mw=1,mwaves
                bmasdq(m,i) = bmasdq(m,i) &
                + dmin1(sb(mw), 0.d0) * waveb(m,mw)
                bpasdq(m,i) = bpasdq(m,i) &
                + dmax1(sb(mw), 0.d0) * waveb(m,mw)
            end do
        end do
    
    20 END DO

    return
    end subroutine rpt2
