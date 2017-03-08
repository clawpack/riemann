subroutine rp1(maxm,num_eqn,num_waves,num_aux,num_ghost,num_cells,ql,qr,auxl,auxr,fwave,s,amdq,apdq)

! Riemann solver for the 1D shallow water equations:

! (h)_t + (h*u)_x = 0
! (hu)_t + (h*u^2 + 1/2*grav*h^2)_x = -grav*h*(b)_x
!
! using the f-wave algorithm and Roe's approximate Riemann solver.

! waves:     2
! equations: 2

! Conserved quantities:
!       1 depth
!       2 momentum

! Auxiliary fields:
!       1 bathymetry

! The gravitational constant grav should be in the common block cparam.
!
! See http://www.clawpack.org/riemann.html for a detailed explanation
! of the Riemann solver API.

    implicit none 
    
    integer, intent(in) :: maxm, num_eqn, num_waves, num_ghost, num_aux, num_cells
    real(kind=8), intent(in) :: ql(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: qr(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: auxl(num_aux, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: auxr(num_aux, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: s(num_waves, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: fwave(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: amdq(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(out) :: apdq(num_eqn, 1-num_ghost:maxm+num_ghost)
    
    real(kind=8) :: grav, dry_tolerance, sea_level
    common /cparam/ grav, dry_tolerance, sea_level
    
    real(kind=8) :: hl, ul, hr, ur, hbar, uhat, chat, bl, br, phil, phir
    real(kind=8) :: R(2,2), fluxdiff(2), beta(2)
    
    integer :: i, k

    amdq = 0.d0
    apdq = 0.d0
    
    do i=2-num_ghost,num_cells+num_ghost
        ! # Left states
        hl = qr(1, i - 1)
        if (hl > dry_tolerance) then
            ul = qr(2, i - 1) / hl
            phil = qr(1, i - 1) * ul**2 + 0.5d0 * grav * qr(1, i - 1)**2
        else
            ul = 0.d0
            phil = 0.d0
        end if
        bl = auxr(1, i - 1)
        
        ! # Right states
        hr = ql(1, i)
        if (hr > dry_tolerance) then
            ur = ql(2, i)/hr
            phir = ql(1, i) * ur**2 + 0.5d0 * grav * ql(1, i)**2
        else
            ur = 0.d0
            phir = 0.d0
        end if
        br = auxl(1, i)
        
        ! # Roe average states
        hbar = 0.5 * (hr + hl)
        uhat = (sqrt(hl) * ul + sqrt(hr) * ur) / (sqrt(hl) + sqrt(hr))
        chat = sqrt(grav * hbar)
        
        ! # Flux differences
        fluxdiff(1) = (hr * ur) - (hl * ul)
        fluxdiff(2) = phir - phil + grav * hbar * (br - bl)
        
        ! # Wave speeds
        s(1,i) = min(ul - sqrt(grav * hl), uhat - chat)
        s(2,i) = max(ur + sqrt(grav * hr), uhat + chat)
        
        ! # Right eigenvectors (column)
        R(1,1) = 1.d0
        R(2,1) = s(1, i)
        
        R(1,2) = 1.d0
        R(2,2) = s(2, i)
    
        ! Wave strengths
        beta(1) = (s(2, i) * fluxdiff(1) - fluxdiff(2)) / (s(2, i) - s(1, i))
        beta(2) = (fluxdiff(2) - s(1, i) * fluxdiff(1)) / (s(2, i) - s(1, i))

        ! # Flux waves
        do k=1,num_waves
            fwave(:,k,i) = beta(k) * R(:,k)
        enddo
        
        ! # Fluctuations
        do k=1, num_waves
            if (s(k, i) < 0.d0) then
                amdq(:, i) = amdq(:, i) + fwave(:, k, i)
            else if (s(k, i) > 0.d0) then
                apdq(:, i) = apdq(:, i) + fwave(:, k, i)
            else
                amdq(:, i) = amdq(:, i) + 0.5d0 * fwave(:, k, i)
                apdq(:, i) = apdq(:, i) + 0.5d0 * fwave(:, k, i)
            end if
        enddo
    enddo
    
end subroutine rp1
