! ==============================================================================
!  Point-wise version of the acoustics problem
! ==============================================================================

subroutine rp1_ptwise(num_eqn, num_aux, num_waves, q_l, q_r, aux_l, aux_r,    &
                      wave, s, amdq, apdq)

    implicit none

    ! Input Arguments
    integer, intent(in) :: num_eqn, num_aux, num_waves
    real(kind=8), intent(in out) :: q_l(num_eqn), q_r(num_eqn)
    real(kind=8), intent(in out) :: aux_l(num_aux), aux_r(num_aux)

    ! Output arguments
    real(kind=8), intent(in out) :: wave(num_eqn, num_waves)
    real(kind=8), intent(in out) :: s(num_waves)
    real(kind=8), intent(in out) :: apdq(num_eqn), amdq(num_eqn)

    ! Locals
    real(kind=8) :: delta(2), a(2)

    ! Common block
    real(kind=8) :: rho, bulk, cc, zz

    common /cparam/ rho, bulk, cc, zz

    delta = q_r - q_l

    a(1) = (-delta(1) + zz * delta(2)) / (2.d0 * zz)
    a(2) = ( delta(1) + zz * delta(2)) / (2.d0 * zz)

    wave(1, 1) = -a(1) * zz
    wave(2, 1) =  a(1)
    s(1) = -cc

    wave(1, 2) = a(2) * zz
    wave(2, 2) = a(2)
    s(2) = cc

    amdq = s(1) * wave(:, 1)
    apdq = s(2) * wave(:, 2)

end subroutine rp1_ptwise