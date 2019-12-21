! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================
!     # Adjoint Riemann solver for the acoustics equations in 2d, with varying
!     # material properties rho and kappa

!     waves: 2
!     equations: 3
!     aux fields: 2

!     Conserved quantities:
!       1 pressure
!       2 x_momentum
!       3 y_momentum

!     Auxiliary variables:
!         1  density
!         2  sound_speed

!     # Note that although there are 3 eigenvectors, the second eigenvalue
!     # is always zero and so we only need to compute 2 waves.
!     #
!     # solve Riemann problems along one slice of data.

!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell

!     # Here it is assumed that auxl=auxr gives the cell values.

!     # This data is along a slice in the x-direction if ixy=1
!     #                            or the y-direction if ixy=2.
!     # On output, wave contains the waves,
!     #            s the speeds,
!     #            amdq the  left-going flux difference  A^- \Delta q
!     #            apdq the right-going flux difference  A^+ \Delta q


!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     # From the basic clawpack routines, this routine is called with ql = qr

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
!     # The jump is split into a leftgoing wave traveling at speed -c
!     # relative to the material properties to the left of the interface,
!     # and a rightgoing wave traveling at speed +c
!     # relative to the material properties to the right of the interface,

!     # find b1 and b2, the coefficients of the 2 eigenvectors:
    do 20 i = 2-mbc, mx+mbc
!       # material properties

!        # impedances:
        zi = auxl(1,i)*auxl(2,i)
        zim = auxr(1,i-1)*auxr(2,i-1)

!       # sound speed
        ci = auxl(2,i)
        cim = auxr(2,i-1)

!       # density
        rhoi = auxl(1,i)
        rhoim = auxr(1,i-1)

!       # bulk modulus
        bulki = rhoi*ci**2
        bulkim = rhoim*cim**2

!       # f-wave splitting
        delta(1) = -(ql(mu,i)/rhoi - qr(mu,i-1)/rhoim)
        delta(2) = -(bulki*ql(1,i) - bulkim*qr(1,i-1))

        beta1 = (zi*delta(1) + delta(2)) / (zi+zim)
        beta2 = (zim*delta(1) - delta(2)) / (zi+zim)
    
!       # Compute the waves. Eigenvectors switched for adjoint:

        fwave(1,1,i) = beta1
        fwave(mu,1,i) = beta1*zim
        fwave(mv,1,i) = 0.d0
        s(1,i) = -cim
    
        fwave(1,2,i) = beta2
        fwave(mu,2,i) = -beta2*zi
        fwave(mv,2,i) = 0.d0
        s(2,i) = ci
    
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
