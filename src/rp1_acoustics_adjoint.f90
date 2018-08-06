! =====================================================
subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================
!     # Adjoint Riemann solver for the acoustics equations in 1d.

! waves:     2
! equations: 2

! Conserved quantities:
!       1 pressure
!       2 velocity

!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell

!     # On output, wave contains the waves,
!     #            s the speeds,
!     #            amdq the  left-going flux difference  A^- \Delta q
!     #            apdq the right-going flux difference  A^+ \Delta q


!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     # From the basic clawpack routines, this routine is called with ql = qr

!     # aux arrays not used in this solver.

    implicit double precision (a-h,o-z)

    dimension fwave(meqn, mwaves, 1-mbc:maxm+mbc)
    dimension    s(mwaves, 1-mbc:maxm+mbc)
    dimension   ql(meqn, 1-mbc:maxm+mbc)
    dimension   qr(meqn, 1-mbc:maxm+mbc)
    dimension apdq(meqn, 1-mbc:maxm+mbc)
    dimension amdq(meqn, 1-mbc:maxm+mbc)
    dimension auxl(maux, 1-mbc:maxm+mbc)
    dimension auxr(maux, 1-mbc:maxm+mbc)

!     local arrays
!     ------------
    dimension delta(2)

!     # density, bulk modulus, and sound speed, and impedence of medium:
!     # (should be set in setprob.f)
    common /cparam/ rho,bulk,cc,zz

!     # split the jump in f(q) at each interface into waves

!     # find b1 and b2, the coefficients of the 2 eigenvectors:
    do 20 i = 2-mbc, mx+mbc
!       # material properties
        rhoi = zz/cc
        bulki = rhoi*cc**2

!       # f-wave splitting
        delta(1) = -(ql(2,i)/rhoi - qr(2,i-1)/rhoi)
        delta(2) = -(bulki*ql(1,i) - bulki*qr(1,i-1))

        beta1 = (zz*delta(1) + delta(2)) / (2.d0*zz)
        beta2 = (zz*delta(1) - delta(2)) / (2.d0*zz)
    
!       # Compute the waves. Eigenvectors switched for adjoint:
!       # (speeds negative of eigenvalues because we are going backward in time)

        fwave(1,1,i) = beta1
        fwave(2,1,i) = beta1*zz
        s(1,i) = -cc
    
        fwave(1,2,i) = beta2
        fwave(2,2,i) = -beta2*zz
        s(2,i) = cc
    
    20 END DO


!     # compute the leftgoing and rightgoing flux differences:
!     # Note s(i,1) < 0   and   s(i,2) > 0.

!   # f-wave splitting, so do not multiply by wave speeds:
    forall (m=1:meqn,  i=2-mbc: mx+mbc)
    amdq(m,i) = fwave(m,1,i)
    apdq(m,i) = fwave(m,2,i)
    end forall

    return
    end subroutine rp1
