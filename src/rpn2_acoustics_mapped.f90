! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================
! Riemann solver for the acoustics equations in 2d 
! on general quadrilateral grid, with variable coefficients

! waves: 2
! equations: 3

! Conserved quantities:
!       1 pressure
!       2 x_velocity
!       3 y_velocity

! Auxiliary quantities:
!       1 a_x
!       2 a_y   where (a_x,a_y) is unit normal to left face
!       3 length_ratio_left ratio of length of left face to dyc
!       4 b_x
!       5 b_y   where (b_x,b_y) is unit normal to bottom face
!       6 length_ratio_bottom   ratio of length of bottom face to dxc
!       7 cell_area   ratio of cell area to dxc*dyc
!         (approximately Jacobian of mapping function)
!       8 Z (impedance)
!       9 c (sound speed)
!
! Rotate velocity and then call standard Riemann solver.
! The resulting waves and flux differences are then rotated
! back to x-y.

! solve Riemann problems along one slice of data.
!
! On input, ql contains the state vector at the left edge of each cell
! qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1 
! or the y-direction if ixy=2.
! On output, wave contains the waves, s the speeds, 
! and amdq, apdq the decomposition of the flux difference
! f(qr(i-1)) - f(ql(i))  
! into leftgoing and rightgoing parts respectively.
! With the Roe solver we have   
!     amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated
! into the flux differences.
!
! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                   and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr.

    implicit none
!
    integer, intent(in) :: maxm, meqn, mwaves, mbc, mx
    double precision, intent(out) :: wave(meqn, mwaves, 1-mbc:maxm+mbc)
    double precision, intent(out) :: s(mwaves,  1-mbc:maxm+mbc)
    double precision, intent(in)  ::   ql(meqn, 1-mbc:maxm+mbc)
    double precision, intent(in)  ::   qr(meqn, 1-mbc:maxm+mbc)
    double precision, intent(out) :: apdq(meqn, 1-mbc:maxm+mbc)
    double precision, intent(out) :: amdq(meqn, 1-mbc:maxm+mbc)
    double precision, intent(in)  :: auxl(9,    1-mbc:maxm+mbc)
    double precision, intent(in)  :: auxr(9,    1-mbc:maxm+mbc)
!
!   ------------
    double precision :: delta(3)
    double precision :: unorl, unorr, a1, a2
    double precision :: alpha, beta, zi, zim, ci, cim
    integer :: ixy, inx, iny, ilenrat, i, m
!
! Rotate the velocity vector (q(2), q(3)) so that it is aligned with the face
! normal.  The normal vector for the face at the i'th Riemann problem
! is stored in the aux array
! in locations (1,2) if ixy=1 or (4,5) if ixy=2.  The ratio of the
! length of the cell side to the length of the computational cell
! is stored in aux(3) or aux(6), respectively.

    if (ixy.eq.1) then
        inx = 1
        iny = 2
        ilenrat = 3
    else
        inx = 4
        iny = 5
        ilenrat = 6
    endif

! Determine rotation matrix:
!               [ alpha  beta ]
!               [-beta  alpha ]
! Note that this reduces to the identity on a standard cartesian grid.

! Determine normal velocity components at this edge:
    do i=2-mbc,mx+mbc
        alpha = auxl(i,inx)
        beta  = auxl(i,iny)
        unorl = alpha*ql(i,2) + beta*ql(i,3)
        unorr = alpha*qr(i-1,2) + beta*qr(i-1,3)
    enddo


    do i = 2-mbc, mx+mbc
        delta(1) = ql(i,1) - qr(i-1,1)
        delta(2) = unorl - unorr

        zi  = auxl(i,8)
        zim = auxl(i-1,8)
        ci  = auxl(i,9)
        cim = auxl(i-1,9)

        a1 = (-delta(1) + zi*delta(2)) / (zim + zi)
        a2 =  (delta(1) + zim*delta(2)) / (zim + zi)

!       Scale the velocities by the length ratios on each side.

        wave(i,1,1) = -a1*zim
        wave(i,2,1) = a1 * alpha
        wave(i,3,1) = a1 * beta
        s(i,1) = -cim * auxl(i,ilenrat)

        wave(i,1,2) = a2*zi
        wave(i,2,2) = a2 * alpha
        wave(i,3,2) = a2 * beta
        s(i,2) = ci * auxl(i,ilenrat)
    enddo

!   Compute the left-going and right-going fluctuations.
!   Note that s(i,1) < 0   and   s(i,2) > 0.

    do m=1,meqn
        do i = 2-mbc, mx+mbc
            amdq(i,m) = s(i,1)*wave(i,m,1)
            apdq(i,m) = s(i,2)*wave(i,m,2)
        enddo
    enddo

    return
    end subroutine rpn2
