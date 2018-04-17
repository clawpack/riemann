!   =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
!   =====================================================
    implicit double precision (a-h,o-z)

!     # Riemann solver in the transverse direction for the acoustics equations.

!     # auxN(1,i) holds rho
!     # auxN(2,i) holds c
!     #  N = 1 for row below
!     #      2 for this row
!     #      3 for row above

!     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.

    dimension     ql(meqn, 1-mbc:maxm+mbc)
    dimension     qr(meqn, 1-mbc:maxm+mbc)
    dimension   asdq(meqn, 1-mbc:maxm+mbc)
    dimension bmasdq(meqn, 1-mbc:maxm+mbc)
    dimension bpasdq(meqn, 1-mbc:maxm+mbc)
    dimension   aux1(maux, 1-mbc:maxm+mbc)
    dimension   aux2(maux, 1-mbc:maxm+mbc)
    dimension   aux3(maux, 1-mbc:maxm+mbc)


    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do 20 i = 2-mbc, mx+mbc

!        # imp is used to flag whether wave is going to left or right,
!        # since material properties are different on the two sides

        if (imp == 1) then
!            # asdq = amdq, moving to left
            i1 = i-1
        else
!            # asdq = apdq, moving to right
            i1 = i
        endif

!        # sound speed in each row of cells:
        cm = aux1(2,i1)
        c = aux2(2,i1)
        cp = aux3(2,i1)

!        # impedances:
        zm = aux1(1,i1)*aux1(2,i1)
        zz = aux2(1,i1)*aux2(2,i1)
        zp = aux3(1,i1)*aux3(2,i1)

        a1 = (asdq(1,i)*zz - asdq(mv,i)) / (zp + zz)
        a2 = (zz*asdq(1,i) + asdq(mv,i)) / (zz + zm)
    
    !        # The down-going flux difference bmasdq is the product  -c * wave
    
        bmasdq(1,i) = -cc * a2
        bmasdq(mu,i) = 0.d0
        bmasdq(mv,i) = -cc * a2 * zm
    
    !        # The up-going flux difference bpasdq is the product  c * wave
    
        bpasdq(1,i) = cc * a1
        bpasdq(mu,i) = 0.d0
        bpasdq(mv,i) = -cc * a1 * zp
    
    20 END DO

    return
    end subroutine rpt2
