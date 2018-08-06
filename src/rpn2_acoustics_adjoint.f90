! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================
!     # Adjoint Riemann solver for the acoustics equations in 2d,

!     # Note that although there are 3 eigenvectors, the second eigenvalue
!     # is always zero and so we only need to compute 2 waves.
!     #
!     # solve Riemann problems along one slice of data.

!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell

!     # This data is along a slice in the x-direction if ixy=1
!     #                            or the y-direction if ixy=2.
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
    dimension delta(3)

!     # density, bulk modulus, and sound speed, and impedence of medium:
!     # (should be set in setprob.f)
    common /cparam/ rho,bulk,cc,zz



!     # set mu to point to  the component of the system that corresponds
!     # to velocity in the direction of this slice, mv to the orthogonal
!     # velocity:

    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

!     # note that notation for u and v reflects assumption that the
!     # Riemann problems are in the x-direction with u in the normal
!     # direciton and v in the orthogonal direcion, but with the above
!     # definitions of mu and mv the routine also works with ixy=2
!     # in which case waves come from the
!     # Riemann problems u_t + g(u)_y = 0 in the y-direction.


!     # split the jump in f(q) at each interface into waves

!     # find b1 and b2, the coefficients of the 2 eigenvectors:
    do 20 i = 2-mbc, mx+mbc
!       # material properties
        rhoi = zz/cc
        bulki = rhoi*cc**2

!       # f-wave splitting
        delta(1) = -(ql(mu,i)/rhoi - qr(mu,i-1)/rhoi)
        delta(2) = -(bulki*ql(1,i) - bulki*qr(1,i-1))

        beta1 = (zz*delta(1) + delta(2)) / (2.d0*zz)
        beta2 = (zz*delta(1) - delta(2)) / (2.d0*zz)
    
!       # Compute the waves. Eigenvectors switched for adjoint:
!       # (speeds negative of eigenvalues because we are going backward in time)

        fwave(1,1,i) = beta1
        fwave(mu,1,i) = beta1*zz
        fwave(mv,1,i) = 0.d0
        s(1,i) = -cc
    
        fwave(1,2,i) = beta2
        fwave(mu,2,i) = -beta2*zz
        fwave(mv,2,i) = 0.d0
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
    end subroutine rpn2
