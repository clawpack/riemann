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
    
    real(kind=8) :: grav, dry_tolerance, sea_level
    common /cparam/ grav, dry_tolerance, sea_level
    
    real(kind=8) :: hl, ul, vl, bl, hr, ur, vr, br, hbar, uhat, vhat, chat
    real(kind=8) :: phil, phir, n_1, n_2, dry_state_l, dry_state_r
    real(kind=8) :: R(3,3)
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
        dry_state_l = merge(1.d0, 0.d0, qr(1, i - 1) < dry_tolerance)
        dry_state_r = merge(1.d0, 0.d0, ql(1, i) < dry_tolerance)

        ! Note that for the states below u is always the normal velocity and
        ! v is always the tangential velocity

        ! # Left states
        hl = qr(1, i - 1)
        ul = ((qr(2, i - 1) / hl) * n_1 + (qr(3, i - 1) / hl) * n_2) * dry_state_l
        vl = ((qr(2, i - 1) / hl) * n_2 + (qr(3, i - 1) / hl) * n_1) * dry_state_l
        phil = (0.5d0 * grav * hl**2 + hl * ul**2) * dry_state_l

        bl = auxr(1,i-1)
    
        ! # Right states
        hr = ql(1, i)
        ur = ((ql(2, i) / hr) * n_1 + (ql(3, i) / hr) * n_2) * dry_state_r
        vr = ((ql(2, i) / hr) * n_2 + (ql(3, i) / hr) * n_1) * dry_state_r
        phir = (0.5d0 * grav * hr**2 + hr * ur**2) * dry_state_r

        br = auxl(1, i)
    
        ! # Roe average states (Roe's linearization)
        hbar = 0.5d0 * (hr + hl)
        uhat = (sqrt(hr) * ur + sqrt(hl) * ul) / (sqrt(hr) + sqrt(hl)) 
        vhat = (sqrt(hr) * vr + sqrt(hl) * vl) / (sqrt(hr) + sqrt(hl))
        chat = sqrt(grav * hbar)
    
        ! # Flux differences
        ! ##################
        ! # In order to compute these differences in an efficient way (avoid if statements)
        ! # we compute the scalar product between the velocity vector at the interface and the
        ! # normal vector "{n} = {n_1,n_2}^T" which has been defined above.
        ! fluxdiff(1) = hr*(ur*n_1 + vr*n_2) - hl*(ul*n_1 + vl*n_2) ! Easy to understand
        fluxdiff(1) = hr * ur - hl * ul

        fluxdiff(2) = phir - phil + grav * hbar * (br - bl) ! Always normal

        fluxdiff(3) = hr * ur * vr - hl * ul * vl ! Always tangential
    
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
        ! fluxdiff(2) = (hr*ur*(ur*n_1 + vr*n_2)+0.5*grav*hr*hr*n_1)-(hl*ul*(ul*n_1 + vl*n_2)+0.5*grav*hl*hl*n_1) &
        ! & + grav*hbar*(br-bl)*n_1
        
        ! fluxdiff(3) = (hr*vr*(ur*n_1 + vr*n_2)+0.5*grav*hr*hr*n_2)-(hl*vl*(ul*n_1 + vl*n_2)+0.5*grav*hl*hl*n_2) &
        ! & + grav*hbar*(br-bl)*n_2
    
        ! # Wave speeds
        s(1, i) = min(uhat - chat, ul - sqrt(grav * hl))
        s(3, i) = max(uhat + chat, ur + sqrt(grav * hr))
        s(2, i) = 0.5d0 * (s(1, i) + s(3, i))

        ! s(1,i) = (uhat*n_1 + vhat*n_2) - chat
        ! s(2,i) = (uhat*n_1 + vhat*n_2)
        ! s(3,i) = (uhat*n_1 + vhat*n_2) + chat
        
        !write(*,*) i,s(i,1),s(i,2),s(i,3)
        
        ! # Right eigenvectors (columns)
        R(1, 1) = 1.d0
        R(2, 1) = s(1, i) * n_1 + vhat * n_2
        R(3, 1) = s(1, i) * n_1 + vhat * n_2
        ! R(2,1) = uhat - chat*n_1
        ! R(3,1) = vhat - chat*n_2
        
        R(1, 2) = 0.d0
        R(2, 2) = -1.d0*n_2
        R(3, 2) = n_1
        
        R(1, 3) = 1.d0
        R(2, 3) = s(3, i) * n_1 + vhat * n_2
        R(3, 3) = s(3, i) * n_1 + vhat * n_2
        ! R(2,3) = uhat + chat*n_1
        ! R(3,3) = vhat + chat*n_2
        
        ! # Left eigenvectors (rows)
        ! L(1,1) = ((uhat*n_1 + vhat*n_2) + chat)/(2.d0*chat)
        ! L(1,2) = -1.d0*n_1/(2.d0*chat)
        ! L(1,3) = -1.d0*n_2/(2.d0*chat)
        
        ! L(2,1) = uhat*n_2 - vhat*n_1 
        ! L(2,2) = -1.d0*n_2
        ! L(2,3) = n_1
        
        ! L(3,1) = (chat - (uhat*n_1 + vhat*n_2))/(2.d0*chat)
        ! L(3,2) = n_1/(2.d0*chat)
        ! L(3,3) = n_2/(2.d0*chat)
        
        
        ! # Coefficients beta which multiply the right eigenvectors (see below)
        ! do j=1,num_eqn
        !     beta(j) = 0.d0
        !     do k=1,num_eqn
        !         beta(j) = beta(j) + L(j,k)*fluxdiff(k)
        !     enddo
        ! enddo
        beta(1) = (s(3, i) * fluxdiff(1) - fluxdiff(2)) / (s(3, i) - s(1, i))
        beta(3) = (fluxdiff(2) - s(1, i) * fluxdiff(1)) / (s(3, i) - s(1, i))
        beta(2) = hr * ur * vr - hl * ul * vl - beta(1) * vl - beta(3) * vr

        ! # Flux waves
        do k=1,num_eqn
            fwave(:, k, i) = beta(k) * R(:, k)
        enddo
    
        ! # Fluctuations
        do k=1, num_waves
            if (s(k, i) < 1.0e-14) then
                amdq(:, i) = amdq(:, i) + fwave(:, k, i)
            elseif (s(k, i) > 1.0e-14) then
                apdq(:, i) = apdq(:, i) + fwave(:, k, i)
            else
                amdq(:, i) = amdq(:, i) + 0.5d0 * fwave(:, k, i)
                apdq(:, i) = apdq(:, i) + 0.5d0 * fwave(:, k, i)
            endif
        enddo
    
    enddo
    
    return
end subroutine rpn2
