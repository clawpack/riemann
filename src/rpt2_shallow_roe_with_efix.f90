

!     =====================================================
    subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx, &
    ql,qr,aux1,aux2,aux3, &
    ilr,asdq,bmasdq,bpasdq,num_aux)
!     =====================================================
    implicit double precision (a-h,o-z)

!     # Riemann solver in the transverse direction for the shallow water
!     equations .
!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)

!     # Uses Roe averages and other quantities which were
!     # computed in rpn2sh and stored in the common block comroe.

    dimension     ql(meqn, 1-mbc:maxm+mbc)
    dimension     qr(meqn, 1-mbc:maxm+mbc)
    dimension   asdq(meqn, 1-mbc:maxm+mbc)
    dimension bmasdq(meqn, 1-mbc:maxm+mbc)
    dimension bpasdq(meqn, 1-mbc:maxm+mbc)

    dimension waveb(3,3),sb(3)
    parameter (maxm2 = 1800)
    common /comroe/ u(-6:maxm2+7),v(-6:maxm2+7),a(-6:maxm2+7), &
    h(-6:maxm2+7)


    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif
            

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
    
    !           # compute the flux differences bmasdq and bpasdq
    
        do 10 m=1,meqn
            bmasdq(m,i) = 0.d0
            bpasdq(m,i) = 0.d0
            do 10 mw=1,mwaves
                bmasdq(m,i) = bmasdq(m,i) &
                + dmin1(sb(mw), 0.d0) * waveb(m,mw)
                bpasdq(m,i) = bpasdq(m,i) &
                + dmax1(sb(mw), 0.d0) * waveb(m,mw)
        10 END DO
    
    20 END DO

    return
    end subroutine rpt2
