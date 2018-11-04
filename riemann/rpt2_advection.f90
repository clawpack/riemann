! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
    implicit double precision (a-h,o-z)

!     # Riemann solver in the transverse direction for the scalar equation

!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)


    dimension     ql(meqn,1-mbc:maxm+mbc)
    dimension     qr(meqn,1-mbc:maxm+mbc)
    dimension   asdq(meqn,1-mbc:maxm+mbc)
    dimension bmasdq(meqn,1-mbc:maxm+mbc)
    dimension bpasdq(meqn,1-mbc:maxm+mbc)
    common /cparam/ u,v

!     # transverse wave speeds have been computed in rpn2
!     # maux=0 and aux arrays are unused in this example.

    if (ixy == 1) then
        stran = v
    else
        stran = u
    endif

    stranm = dmin1(stran, 0.d0)
    stranp = dmax1(stran, 0.d0)
          
    do 10 i = 2-mbc, mx+mbc
        bmasdq(1,i) = stranm * asdq(1,i)
        bpasdq(1,i) = stranp * asdq(1,i)
    10 END DO

    return
    end subroutine rpt2
