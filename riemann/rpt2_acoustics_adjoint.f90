!   =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
!   =====================================================
    implicit double precision (a-h,o-z)

!     # Riemann solver in the transverse direction for the acoustics equations.

!     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.

    dimension     ql(meqn, 1-mbc:maxm+mbc)
    dimension     qr(meqn, 1-mbc:maxm+mbc)
    dimension   asdq(meqn, 1-mbc:maxm+mbc)
    dimension bmasdq(meqn, 1-mbc:maxm+mbc)
    dimension bpasdq(meqn, 1-mbc:maxm+mbc)

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
