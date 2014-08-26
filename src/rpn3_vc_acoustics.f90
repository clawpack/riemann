! ==================================================================
subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! ==================================================================

! Riemann solver for the acoustics equations in 3d, with varying
! material properties.
!
! waves: 2
! equations: 4
! aux fields: 2
!
! Conserved quantities:
!       1 pressure
!       2 x_velocity
!       3 y_velocity
!       4 z_velocity
!
! Auxiliary variables:
!       1 impedance
!       2 sound_speed

! Note that although there are 4 eigenvectors, two eigenvalues are
! always zero and so we only need to compute 2 waves.

! Solve Riemann problems along one slice of data.
! This data is along a slice in the x-direction if ixyz=1
!                               the y-direction if ixyz=2.
!                               the z-direction if ixyz=3.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr


    implicit none
    double precision :: wave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision ::    s(mwaves,1-mbc:maxm+mbc)
    double precision ::   ql(meqn,1-mbc:maxm+mbc)
    double precision ::   qr(meqn,1-mbc:maxm+mbc)
    double precision :: amdq(meqn,1-mbc:maxm+mbc)
    double precision :: apdq(meqn,1-mbc:maxm+mbc)
    double precision :: auxl(maux,1-mbc:maxm+mbc)
    double precision :: auxr(maux,1-mbc:maxm+mbc)
    double precision :: delta(3), zi, zim, a1, a2 
    integer :: ixyz, mu, mv, mw, maxm,meqn,mwaves,mbc,mx, maux
    integer :: i,m



!     # set mu to point to  the component of the system that corresponds
!     # to velocity in the direction of this slice, mv to the orthogonal
!     # velocity.


    if (ixyz == 1) then
        mu = 2
        mv = 3
        mw = 4
    else if (ixyz == 2) then
        mu = 3
        mv = 4
        mw = 2
    else if (ixyz == 3) then
        mu = 4
        mv = 2
        mw = 3
    endif


!     # split the jump in q at each interface into waves
!     # The jump is split into a leftgoing wave traveling at speed -c
!     # relative to the material properties to the left of the interface,
!     # and a rightgoing wave traveling at speed +c
!     # relative to the material properties to the right of the interface,

!     # find a1 and a2, the coefficients of the 2 eigenvectors:
    do 20 i = 2-mbc, mx+mbc
        delta(1) = ql(1,i) - qr(1,i-1)
        delta(2) = ql(mu,i) - qr(mu,i-1)
    !        # impedances:
        zi = auxl(1,i)
        zim = auxr(1,i-1)

        a1 = (-delta(1) + zi*delta(2)) / (zim + zi)
        a2 =  (delta(1) + zim*delta(2)) / (zim + zi)


    
    !        # Compute the waves.
    
        wave(1,1,i) = -a1*zim
        wave(mu,1,i) = a1
        wave(mv,1,i) = 0.d0
        wave(mw,1,i) = 0.d0
        s(1,i) = -auxr(2,i-1)
    
        wave(1,2,i) = a2*zi
        wave(mu,2,i) = a2
        wave(mv,2,i) = 0.d0
        wave(mw,2,i) = 0.d0
        s(2,i) = auxl(2,i)
    
    20 END DO


!     # compute the leftgoing and rightgoing flux differences:
!     # Note s(i,1) < 0   and   s(i,2) > 0.

    do 220 m=1,meqn
        do 220 i = 2-mbc, mx+mbc
            amdq(m,i) = s(1,i)*wave(m,1,i)
            apdq(m,i) = s(2,i)*wave(m,2,i)
    220 END DO


    return
    end subroutine rpn3

