! ==============================================================================
! Riemann solver for the acoustics equations in 2d.
!
! waves: 3
! equations: 3
!
! Conserved quantities:
!       1 pressure
!       2 x_velocity
!       3 y_velocity
!
! Note that although there are 3 eigenvectors, the second eigenvalue
! is always zero and so we only need to compute 2 waves.
!
! solve Riemann problems along one slice of data.

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q
!
! ==============================================================================

subroutine rpn2_ptwise(ixy, num_eqn, num_aux, num_waves, q_l, q_r,            &
                       aux_l, aux_r, wave, s, amdq, apdq)

    implicit none

    ! Input
    integer, intent(in) :: num_eqn, num_aux, num_waves, ixy
    real(kind=8), intent(in) :: q_l(num_eqn), q_r(num_eqn)
    real(kind=8), intent(in) :: aux_l(num_aux), aux_r(num_aux)

    ! Output
    real(kind=8), intent(in out) :: wave(num_eqn, num_waves)
    real(kind=8), intent(in out) :: s(num_waves)
    real(kind=8), intent(in out) :: apdq(num_eqn), amdq(num_eqn)

    ! Locals
    integer :: normal_index, transverse_index
    real(kind=8) :: delta(3), a(2)

    ! Common block
    real(kind=8) :: rho, bulk, cc, zz
    common /cparam/ rho, bulk, cc, zz

    if (ixy == 1) then
        normal_index = 2
        transverse_index = 3
    else
        normal_index = 3
        transverse_index = 2
    end if

    ! note that notation for u and v reflects assumption that the
    ! Riemann problems are in the x-direction with u in the normal
    ! direciton and v in the orthogonal direcion, but with the above
    ! definitions of mu and mv the routine also works with ixy=2
    ! in which case waves come from the
    ! Riemann problems u_t + g(u)_y = 0 in the y-direction.


    ! Split the jump in q at each interface into waves
    delta = q_r - q_l
    ! Find coefficients for the 2 eigenvectors:
    a(1) = (-delta(1) + zz * delta(normal_index)) / (2.d0 * zz)
    a(2) = ( delta(1) + zz * delta(normal_index)) / (2.d0 * zz)

    ! Compute waves
    wave(1, 1) = -a(1) * zz
    wave(normal_index, 1) = a(1)
    wave(transverse_index, 1) = 0.d0
    s(1) = -cc

    wave(1, 2) = a(2) * zz
    wave(normal_index, 2) = a(2)
    wave(transverse_index, 2) = 0.d0
    s(2) = cc

    ! Compute the leftgoing and rightgoing flux differences:
    ! Note s(1) < 0   and   s(2) > 0.
    amdq = s(1) * wave(:, 1)
    apdq = s(2) * wave(:, 2)

end subroutine rpn2_ptwise
