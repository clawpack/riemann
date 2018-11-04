subroutine rpt2(ixy, imp, maxm, num_eqn, num_waves, num_aux, num_ghost,   &
                num_cells, ql, qr, aux1, aux2, aux3, asdq, bmasdq, bpasdq)
! Transverse solver for shallow water with bathymetry using an fwave solver
! Note that this solver does NOT handle dry-states
!
! (h)_t + (h*u)_x = 0
! (hu)_t + (h*u^2 + 1/2*grav*h^2)_x + (h*u*v)_y = -grav*h*(b)_x
! (hv)_t + (h*u*v)_x + (h*v^2 + 1/2*grav*h^2)_y = -grav*h*(b)_y
!
! using the f-wave algorithm and Roe's approximate Riemann solver.  
!
! waves:     3
! equations: 3
!
! Conserved quantities:
!       1 depth
!       2 x_momentum
!       3 y_momentum
!
! Auxiliary fields:
!       1 bathymetry
!
! The gravitational constant grav should be in the common block cparam.
!
! See http://www.clawpack.org/riemann.html for a detailed explanation
! of the Riemann solver API.

    implicit none

    ! Input
    integer, intent(in) :: ixy, imp, maxm, num_eqn, num_waves, num_aux
    integer, intent(in) :: num_ghost, num_cells
    real(kind=8), intent(in) :: ql(num_eqn, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in) :: qr(num_eqn, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in) :: aux1(num_aux, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in) :: aux2(num_aux, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in) :: aux3(num_aux, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in) :: asdq(num_eqn, 1 - num_ghost:maxm + num_ghost)

    ! Output
    real(kind=8), intent(out) :: bmasdq(num_eqn, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(out) :: bpasdq(num_eqn, 1 - num_ghost:maxm + num_ghost)

    ! Local
    integer :: i, k, normal_index, transverse_index
    real(kind=8) :: h_l, h_r, u_l, u_r, v_l, v_r, h_bar, h_root, u_hat, v_hat
    real(kind=8) :: dry_state_l, dry_state_r, s(3), beta(3), R(3, 3)

    ! Parameters
    real(kind=8) :: grav, dry_tolerance, sea_level
    common /cparam/ grav, dry_tolerance, sea_level

    ! Deterimine direction of solve
    if (ixy == 1) then
        normal_index = 2
        transverse_index = 3
    else
        normal_index = 3
        transverse_index = 2
    end if

    ! Primary loop over cell edges
    do i = 2 - num_ghost, num_cells + num_ghost
        
        ! Check for dry states - need merge here to convert to float
        dry_state_l = merge(0.d0, 1.d0, qr(1, i - 1) < dry_tolerance)
        dry_state_r = merge(0.d0, 1.d0, ql(1, i) < dry_tolerance)

        ! Extract state variables
        h_l = qr(1, i - 1) * dry_state_l
        u_l = qr(normal_index, i - 1) / qr(1, i - 1) * dry_state_l
        v_l = qr(transverse_index, i - 1) / qr(1, i - 1) * dry_state_l

        h_r = ql(1, i) * dry_state_r
        u_r = ql(normal_index, i) / ql(1, i) * dry_state_r
        v_r = ql(transverse_index, i) / ql(1, i) * dry_state_r

        ! Determine speeds for the Jacobian
        h_bar = 0.5d0 * (h_r + h_l)
        h_root = sqrt(h_r) + sqrt(h_l)
        v_hat = (v_r * sqrt(h_r)) / h_root + (v_l * sqrt(h_l)) / h_root
        u_hat = (u_r * sqrt(h_r)) / h_root + (u_l * sqrt(h_l)) / h_root

        s(1) = min(v_hat - sqrt(grav * h_bar), v_l - sqrt(grav * h_l))
        s(3) = max(v_hat + sqrt(grav * h_bar), v_r + sqrt(grav * h_r))
        s(2) = 0.5d0 * (s(1) + s(2))

        ! Determine asdq decompisition
        beta(1) = s(3) * asdq(1, i) / (s(3) - s(1))                 &
                                - asdq(transverse_index, i) / (s(3) - s(1))
        beta(2) = -s(2) * asdq(1, i) + asdq(normal_index, i)
        beta(3) = asdq(transverse_index, i) / (s(3) - s(1))         &
                                - (s(1) * asdq(1, i) / (s(3) - s(1)))

        ! Eigenvectors
        R(1, :) = [1.d0, 0.d0, 1.d0]
        R(2, :) = [s(2), 1.d0, s(2)]
        R(3, :) = [s(1), 0.d0, s(3)]

        ! Compute fluctuations
        do k = 1, num_waves
            if (s(k) < 0.d0) then
                bmasdq(1, i) = bmasdq(1, i) + s(k) * beta(k) * R(1, k)
                bmasdq(normal_index, i) = bmasdq(normal_index, i)              &
                                                + s(k) * beta(k) * R(2, k)
                bmasdq(transverse_index, i) = bmasdq(transverse_index, i)      &
                                                + s(k) * beta(k) * R(3, k)
            else
                bpasdq(1, i) = bpasdq(1, i) + s(k) * beta(k) * R(1, k)
                bpasdq(normal_index, i) = bpasdq(normal_index, i)              &
                                                + s(k) * beta(k) * R(2, k)
                bpasdq(transverse_index, i) = bpasdq(transverse_index, i)      &
                                                + s(k) * beta(k) * R(3, k)
            end if
        end do

    end do  ! End of primary loop

end subroutine rpt2