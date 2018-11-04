! ==============================================================================
!  Point-wise version of a simple scalar advection problem
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

    ! Common block
    real(kind=8) :: u

    common /cparam/ u

    s = u
    wave(1, 1) = q_r(1) - q_l(1)

    amdq = min(u, 0.0) * wave(1, 1)
    apdq = max(u, 0.0) * wave(1, 1)

end subroutine rp1_ptwise