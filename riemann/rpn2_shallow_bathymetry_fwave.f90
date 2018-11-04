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
    
    real(kind=8) :: hl, ul, vl, bl, hr, ur, vr, br, hbar, uhat, chat
    real(kind=8) :: phil, phir, dry_state_l, dry_state_r
    real(kind=8) :: R(3,3)
    real(kind=8) :: fluxdiff(3), beta(3)
    
    integer :: i, k, normal_index, transverse_index
    
    ! Determine normal and tangential directions
    if (ixy == 1) then
        normal_index = 2
        transverse_index = 3
    else
        normal_index = 3
        transverse_index = 2
    end if

    amdq = 0.d0
    apdq = 0.d0

    ! Primary loop over each cell
    do i = 2 - num_ghost, num_cells + num_ghost
        
        ! Check for dry states - need merge here to convert to float
        dry_state_l = merge(0.d0, 1.d0, qr(1, i - 1) < dry_tolerance)
        dry_state_r = merge(0.d0, 1.d0, ql(1, i) < dry_tolerance)

        ! Note that for the states below u is always the normal velocity and
        ! v is always the tangential velocity

        ! Left states
        hl = qr(1, i - 1) * dry_state_l
        ul = qr(normal_index, i - 1) / qr(1, i - 1) * dry_state_l
        vl = qr(transverse_index, i - 1) / qr(1, i - 1) * dry_state_l
        phil = (0.5d0 * grav * hl**2 + hl * ul**2) * dry_state_l

        bl = auxr(1, i - 1)
    
        ! Right states
        hr = ql(1, i) * dry_state_r
        ur = ql(normal_index, i) / ql(1, i) * dry_state_r
        vr = ql(transverse_index, i) / ql(1, i) * dry_state_r
        phir = (0.5d0 * grav * hr**2 + hr * ur**2) * dry_state_r

        br = auxl(1, i)
    
        ! Roe average states (Roe's linearization)
        hbar = 0.5d0 * (hr + hl)
        uhat = (sqrt(hr) * ur + sqrt(hl) * ul) / (sqrt(hr) + sqrt(hl))
        chat = sqrt(grav * hbar)
    
        ! Flux differences
        fluxdiff(1) = hr * ur - hl * ul
        fluxdiff(2) = phir - phil + grav * hbar * (br - bl)
        fluxdiff(3) = hr * ur * vr - hl * ul * vl
    
        ! Wave speeds
        s(1, i) = min(uhat - chat, ul - sqrt(grav * hl))
        s(3, i) = max(uhat + chat, ur + sqrt(grav * hr))
        s(2, i) = 0.5d0 * (s(1, i) + s(3, i))
        
        ! Right eigenvectors (columns)
        ! could possibly use vhat instead of vl and vr
        R(1, 1) = 1.d0
        R(normal_index, 1) = s(1, i)
        R(transverse_index, 1) = vl
        
        R(1, 2) = 0.d0
        R(normal_index, 2) = 0.0
        R(transverse_index, 2) = 1.0
        
        R(1, 3) = 1.d0
        R(normal_index, 3) = s(3, i)
        R(transverse_index, 3) = vr
        
        ! Wave strengths
        beta(1) = (s(3, i) * fluxdiff(1) - fluxdiff(2)) / (s(3, i) - s(1, i))
        beta(3) = (fluxdiff(2) - s(1, i) * fluxdiff(1)) / (s(3, i) - s(1, i))
        beta(2) = fluxdiff(3) - beta(1) * vl - beta(3) * vr

        ! f-waves
        do k = 1, num_waves
            fwave(:, k, i) = beta(k) * R(:, k)
        enddo
    
        ! Fluctuations - could probably rewrite this to be a masking operation
        ! instead...
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

    enddo ! End of main loop

end subroutine rpn2
