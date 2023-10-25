! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================

!     # Riemann solver in the transverse direction.
!     # This is a dummy routine that returns zeros and is only intended
!     # to illustrate the format of this routine.  See various example
!     # directories for better examples.

!     # This dummy routine can be used if transverse solves are not being
!     # used, i.e. if method(3)=0.

!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)


    implicit double precision (a-h,o-z)
    dimension     ql(meqn,1-mbc:maxm+mbc)
    dimension     qr(meqn,1-mbc:maxm+mbc)
    dimension   asdq(meqn,1-mbc:maxm+mbc)
    dimension bmasdq(meqn,1-mbc:maxm+mbc)
    dimension bpasdq(meqn,1-mbc:maxm+mbc)
    dimension   aux1(maux,1-mbc:maxm+mbc)
    dimension   aux2(maux,1-mbc:maxm+mbc)
    dimension   aux3(maux,1-mbc:maxm+mbc)
    dimension s (2,1-mbc:maxm+mbc)

    do 10 i = 2-mbc, mx+mbc
        if (imp == 1) then !left going fluctuation
            ix = i-1
        else !right going fluctuation
            ix = i
        endif

    ! material properties
        pjm=aux1(1,ix) !density at (ix, j-1) where ix=i or i-1
        pj=aux2(1,ix)
        pjp=aux3(1,ix)
        Ejm=aux1(2,ix)
        Ej=aux2(2,ix)
        Ejp=aux3(2,ix)

    ! linearity of the material
        linearity_matjm = nint(aux1(3,ix))
        linearity_matj  = nint(aux2(3,ix))
        linearity_matjp = nint(aux3(3,ix))

    ! eps (strain) at different rows
        epsjm=aux1(4,ix)
        epsj=aux2(4,ix)
        epsjp=aux3(4,ix)

    ! sigmap
        sigmapjm=sigmap(epsjm,Ejm,linearity_matjm)
        sigmapj=sigmap(epsj,Ej,linearity_matj)
        sigmapjp=sigmap(epsjp,Ejp,linearity_matjp)

    ! computation of components of eigenvectors
        r11=1/dsqrt(sigmapjm*pjm)
        r13=-1/dsqrt(sigmapjp*pjp)

    ! shock speeds (eigenvalues of A and B)
        s(1,i)=-dsqrt(sigmapjm/pjm)  !lambda_1, lambda_2=0
        s(2,i)=dsqrt(sigmapjp/pjp)     !lambda_3

        if(ixy == 1) then !x direction, use eigenvectors of matrix B
        ! or the decomposition of the fluctuations
        ! computation of coefficients of the decomposition of the fluctuations.
        ! I called them gamma since beta is used in the decomposition
        ! of the flux difference (notation used in paper for the f-wave)
            gamma1=(asdq(1,i)-r13*asdq(3,i))/(r11-r13)
            gamma3=(-asdq(1,i)+r11*asdq(3,i))/(r11-r13)
        ! computation of the fluctuations
        ! own going part of the normal fluctuation
            bmasdq(1,i)=s(1,i)*gamma1*r11
            bmasdq(2,i)=s(1,i)*gamma1*0
            bmasdq(3,i)=s(1,i)*gamma1*1
        ! p going part of the normal fluctuation
            bpasdq(1,i)=s(2,i)*gamma3*r13
            bpasdq(2,i)=s(2,i)*gamma3*0
            bpasdq(3,i)=s(2,i)*gamma3*1

        else !y direction, use eigenvectors of matrix A
        ! or the decomposition of the fluctuations
            gamma1=(asdq(1,i)-r13*asdq(2,i))/(r11-r13)
            gamma3=(-asdq(1,i)+r11*asdq(2,i))/(r11-r13)
        ! computation of fluctuations
        ! down" going part
            bmasdq(1,i)=s(1,i)*gamma1*r11
            bmasdq(2,i)=s(1,i)*gamma1*1
            bmasdq(3,i)=s(1,i)*gamma1*0
        ! up" going part
            bpasdq(1,i)=s(2,i)*gamma3*r13
            bpasdq(2,i)=s(2,i)*gamma3*1
            bpasdq(3,i)=s(2,i)*gamma3*0
        endif


    10 END DO

    return
    end subroutine rpt2
