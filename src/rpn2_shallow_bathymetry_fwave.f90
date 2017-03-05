subroutine rpn2(ixy, maxm, num_eqn, num_waves, num_aux, num_ghost, &
                num_cells, ql, qr, auxl, auxr, fwave, s, amdq, apdq)

! Riemann solver for the 2D shallow water equations

! (h)_t + (h*u)_x = 0
! (hu)_t + (h*u^2 + 1/2*grav*h^2)_x + (h*u*v)_y = -grav*h*(b)_x
! (hv)_t + (h*u*v)_x + (h*v^2 + 1/2*grav*h^2)_y = -grav*h*(b)_y

! using the f-wave algorithm and Roe's approximate Riemann solver.  

! waves:     3
! equations: 3

! Conserved quantities:
!       1 depth
!       2 x_momentum
!       3 y_momentum

! Auxiliary fields:
!       1 bathymetry

! The gravitational constant grav should be in the common block cparam.
!
! See http://www.clawpack.org/riemann.html for a detailed explanation
! of the Riemann solver API.

    implicit none 
    
    integer, intent(in) :: ixy, maxm, num_eqn, num_waves, num_ghost, num_aux, num_cells
    real(kind=8), intent(in) :: ql(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: qr(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: auxl(num_aux, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: auxr(num_aux, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: s(num_waves, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: fwave(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: amdq(num_eqn,1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: apdq(num_eqn,1-num_ghost:maxm+num_ghost)
    
    real(kind=8) :: grav
    common /cparam/ grav
    
    real(kind=8) :: hl, ul, vl, bl, hr, ur, vr, br, hbar, uhat, vhat, chat, n_1, n_2
    real(kind=8) :: R(3,3), L(3,3)
    real(kind=8) :: fluxdiff(3), beta(3)
    
    integer :: i, j, k
    
    ! Define normal and tangential directions of the current slice. This is needed
    ! to use the right components of the system that correspond to the momentum in
    ! the normal and tangential directions. Since this Riemann solver is used for
    ! both x- and y-slices we do not use names "n" for the normal component and "t"
    ! for the tangential component, but we simply define a vector n = (n_1,n_2)^T
    ! whose component depends on the current slice direction.
    !
    ! See for instance: D. Ambrosi
    !                    Approximation of shallow water equations by Roe's Riemann solver
    !                    International Journal of numerical methods in fluids
    !                    VOL. 20, 157-168 (1995)

    if (ixy .eq. 1) then
        n_1 = 1.d0  ! # This means vertical interface,
        n_2 = 0.d0  ! # e.g. normal vector pointing in the x direction
    else
        n_1 = 0.d0 ! # This means horizontal interface,
        n_2 = 1.d0 ! # e.g. normal vector pointing in the y direction
    endif

    !Now, it is possible to loop over the cells in this slice
    do i=2-num_ghost,num_cells+num_ghost
        ! # Left states
        hl = qr(1,i-1)
        ul = qr(2,i-1)/hl
        vl = qr(3,i-1)/hl
        bl = auxr(1,i-1)
    
        ! # Right states
        hr = ql(1,i)
        ur = ql(2,i)/hr
        vr = ql(3,i)/hr
        br = auxl(1,i)
    
        ! # Roe average states (Roe's linearization)
        hbar = 1.d0/2.d0*(hr + hl)
        uhat = (dsqrt(hr)*ur + dsqrt(hl)*ul)/(dsqrt(hr) + dsqrt(hl)) 
        vhat = (dsqrt(hr)*vr + dsqrt(hl)*vl)/(dsqrt(hr) + dsqrt(hl))
        chat = dsqrt(grav*hbar)
    
        ! # Flux differences
        ! ##################
        ! # In order to compute these differences in an efficient way (avoid if statements)
        ! # we compute the scalar product between the velocity vector at the interface and the
        ! # normal vector "{n} = {n_1,n_2}^T" which has been defined above.
        fluxdiff(1) = hr*(ur*n_1 + vr*n_2) - hl*(ul*n_1 + vl*n_2) ! Easy to understand
    
        ! # If the slice is in the y-direction we have:
        ! # fluxdiff(2) = (hr*ur^2 + 1/2*grav*hr^2) - (hl*ul^2 + 1/2*grav*hl^2)
        ! # fluxdiff(3) = (hr*ur*vr) - (hl*ul*vl)
        ! #
        ! # If the slice is in the x-direction we have:
        ! # fluxdiff(2) = (hr*ur*vr) - (hl*ul*vl)
        ! # fluxdiff(3) = (hr*vr^2 + 1/2*grav*hr^2) - (hl*vl^2 + 1/2*grav*hl^2)
        ! #
        ! # Using the vector component n_1 and n_2 defined above,
        ! # this two possibilities can be achieved in the following way:
        fluxdiff(2) = (hr*ur*(ur*n_1 + vr*n_2)+0.5*grav*hr*hr*n_1)-(hl*ul*(ul*n_1 + vl*n_2)+0.5*grav*hl*hl*n_1) &
        & + grav*hbar*(br-bl)*n_1
        
        fluxdiff(3) = (hr*vr*(ur*n_1 + vr*n_2)+0.5*grav*hr*hr*n_2)-(hl*vl*(ul*n_1 + vl*n_2)+0.5*grav*hl*hl*n_2) &
        & + grav*hbar*(br-bl)*n_2
    
        ! # Wave speeds
        s(1,i) = (uhat*n_1 + vhat*n_2) - chat
        s(2,i) = (uhat*n_1 + vhat*n_2)
        s(3,i) = (uhat*n_1 + vhat*n_2) + chat
        
        !write(*,*) i,s(i,1),s(i,2),s(i,3)
        
        ! # Right eigenvectors (columns)
        R(1,1) = 1.d0
        R(2,1) = uhat - chat*n_1
        R(3,1) = vhat - chat*n_2
        
        R(1,2) = 0.d0
        R(2,2) = -1.d0*n_2
        R(3,2) = n_1
        
        R(1,3) = 1.d0
        R(2,3) = uhat + chat*n_1
        R(3,3) = vhat + chat*n_2
        
        
        ! # Left eigenvectors (rows)
        L(1,1) = ((uhat*n_1 + vhat*n_2) + chat)/(2.d0*chat)
        L(1,2) = -1.d0*n_1/(2.d0*chat)
        L(1,3) = -1.d0*n_2/(2.d0*chat)
        
        L(2,1) = uhat*n_2 - vhat*n_1 
        L(2,2) = -1.d0*n_2
        L(2,3) = n_1
        
        L(3,1) = (chat - (uhat*n_1 + vhat*n_2))/(2.d0*chat)
        L(3,2) = n_1/(2.d0*chat)
        L(3,3) = n_2/(2.d0*chat)
        
        
        ! # Coefficients beta which multiply the right eigenvectors (see below)
        do j=1,num_eqn
            beta(j) = 0.d0
            do k=1,num_eqn
                beta(j) = beta(j) + L(j,k)*fluxdiff(k)
            enddo
        enddo

        ! # Flux waves
        do j=1,num_eqn
            do k=1,num_eqn
                fwave(j,k,i) = beta(k)*R(j,k)
            enddo
        enddo
    
    
        ! # Fluctuations
        do j=1,num_eqn
            amdq(j,i) = 0.d0
            apdq(j,i) = 0.d0
            do k=1,num_eqn
                if (s(k,i) .lt. 1.0e-14) then
                    amdq(j,i) = amdq(j,i) + fwave(j,k,i)
                elseif (s(k,i) .gt. 1.0e-14) then
                    apdq(j,i) = apdq(j,i) + fwave(j,k,i)
                else
                    amdq(j,i) = amdq(j,i) + 1.d0/2.d0*fwave(j,k,i)
                    apdq(j,i) = apdq(j,i) + 1.d0/2.d0*fwave(j,k,i)
                endif
            enddo
        enddo
    
    enddo
    
    return
end subroutine rpn2
