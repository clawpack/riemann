

!     =====================================================
    subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx, &
                ql,qr,aux1,aux2,aux3,ilr,asdq,bmasdq,bpasdq,num_aux)
!     =====================================================
    implicit double precision (a-h,o-z)

!     # Riemann solver in the transverse direction for the Euler equations
!     #  with a tracer variable.
!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)

!     # Uses Roe averages and other quantities which were
!     # computed in rpn2eu and stored in the common block comroe.

    dimension     ql(meqn, 1-mbc:maxm+mbc)
    dimension     qr(meqn, 1-mbc:maxm+mbc)
    dimension   asdq(meqn, 1-mbc:maxm+mbc)
    dimension bmasdq(meqn, 1-mbc:maxm+mbc)
    dimension bpasdq(meqn, 1-mbc:maxm+mbc)

    common /cparam/  gamma,gamma1
    dimension waveb(5,4),sb(4)
    parameter (maxm2 = 1800)
!     # assumes at most maxm2 * maxm2 grid with mbc<=7
    common /comroe/ u2v2(-6:maxm2+7), &
    u(-6:maxm2+7),v(-6:maxm2+7), &
    enth(-6:maxm2+7),a(-6:maxm2+7), &
    g1a2(-6:maxm2+7),euv(-6:maxm2+7)

    if (mbc > 7 .OR. maxm2 < maxm) then
        write(6,*) 'need to increase maxm2 or 7 in rpt'
        stop
    endif

    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do 20 i = 2-mbc, mx+mbc
        a3 = g1a2(i) * (euv(i)*asdq(1,i) &
        + u(i)*asdq(mu,i) + v(i)*asdq(mv,i) - asdq(4,i))
        a2 = asdq(mu,i) - u(i)*asdq(1,i)
        a4 = (asdq(mv,i) + (a(i)-v(i))*asdq(1,i) - a(i)*a3) &
        / (2.d0*a(i))
        a1 = asdq(1,i) - a3 - a4
    
        waveb(1,1) = a1
        waveb(mu,1) = a1*u(i)
        waveb(mv,1) = a1*(v(i)-a(i))
        waveb(4,1) = a1*(enth(i) - v(i)*a(i))
        waveb(5,1) = 0.d0
        sb(1) = v(i) - a(i)
    
        waveb(1,2) = a3
        waveb(mu,2) = a3*u(i) + a2
        waveb(mv,2) = a3*v(i)
        waveb(4,2) = a3*0.5d0*u2v2(i) + a2*u(i)
        waveb(5,2) = 0.d0
        sb(2) = v(i)
    
        waveb(1,3) = a4
        waveb(mu,3) = a4*u(i)
        waveb(mv,3) = a4*(v(i)+a(i))
        waveb(4,3) = a4*(enth(i)+v(i)*a(i))
        waveb(5,3) = 0.d0
        sb(3) = v(i) + a(i)
    
        waveb(1,4) = 0.d0
        waveb(mu,4) = 0.d0
        waveb(mv,4) = 0.d0
        waveb(4,4) = 0.d0
        waveb(5,4) = asdq(5,i)
        sb(4) = v(i)
    
    !           # compute the flux differences bmasdq and bpasdq
    
        do 10 m=1,meqn
            bmasdq(m,i) = 0.d0
            bpasdq(m,i) = 0.d0
            do 10 mw=1,4
                bmasdq(m,i) = bmasdq(m,i) &
                + dmin1(sb(mw), 0.d0) * waveb(m,mw)
                bpasdq(m,i) = bpasdq(m,i) &
                + dmax1(sb(mw), 0.d0) * waveb(m,mw)
        10 END DO
    
    20 END DO

    return
    end subroutine rpt2
