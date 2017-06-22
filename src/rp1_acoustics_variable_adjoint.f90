! =====================================================
subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================

!     # Adjoint Riemann solver for the acoustics equations in 1d, with
!     #  variable coefficients (heterogeneous media)

!     # auxl(1,i) should contain the impedance Z in cell i
!     # auxl(2,i) should contain the sound speed c in cell i

    implicit double precision (a-h,o-z)

    dimension auxl(maux, 1-mbc:maxm+mbc)
    dimension auxr(maux, 1-mbc:maxm+mbc)
    dimension fwave(meqn, mwaves, 1-mbc:maxm+mbc)
    dimension    s(mwaves, 1-mbc:maxm+mbc)
    dimension   ql(meqn, 1-mbc:maxm+mbc)
    dimension   qr(meqn, 1-mbc:maxm+mbc)
    dimension apdq(meqn, 1-mbc:maxm+mbc)
    dimension amdq(meqn, 1-mbc:maxm+mbc)

!     local arrays
!     ------------
    dimension delta(2)


!   # split the jump in f(q) at each interface into f-waves

!   # find a1 and a2, the coefficients of the 2 eigenvectors:
    do 20 i = 2-mbc, mx+mbc
!       # material properties
        zi = auxl(1,i)
        zim = auxr(1,i-1)
        ci = auxl(2,i)
        cim = auxr(2,i-1)
        rhoi = zi/ci
        rhoim = zim/cim
        bulki = rhoi*ci**2
        bulkim = rhoim*cim**2

!       # f-wave splitting
        delta(1) = -(ql(2,i)/rhoi - qr(2,i-1)/rhoim)
        delta(2) = -(bulki*ql(1,i) - bulkim*qr(1,i-1))

        beta1 = (delta(2) + zi*delta(1)) / (zim + zi)
        beta2 = (-delta(2) + zim*delta(1)) / (zim + zi)
    
!       # Compute the waves.  Eigenvectors switched for adjoint:
    
        fwave(1,1,i) = beta1
        fwave(2,1,i) = beta1*zim
        s(1,i) = -cim
    
        fwave(1,2,i) = beta2
        fwave(2,2,i) = -beta2*zi
        s(2,i) = ci
    
    20 END DO

!     # compute the leftgoing and rightgoing fluctuations:
!     # Note s(1,i) < 0   and   s(2,i) > 0.

!   # f-wave spitting, so do not multiply by wave speeds:
    do 220 m=1,meqn
        do 220 i = 2-mbc, mx+mbc
            amdq(m,i) = fwave(m,1,i)
            apdq(m,i) = fwave(m,2,i)
    220 END DO

    return
    end subroutine rp1
