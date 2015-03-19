! ==============================================================================
! Riemann solver in the transverse direction for the acoustics equations.
! Split asdq into down-going flux bmasdq and up-going flux bpasdq.
! ==============================================================================

subroutine rpt2_ptwise(ixy, imp, num_eqn, num_aux, num_waves, q_l, q_r,      &
                       aux_l, aux_r, asdq, bmasdq, bpasdq)

    implicit none

    ! Input
    integer, intent(in) :: num_eqn, num_aux, num_waves, ixy, imp
    real(kind=8), intent(in) :: q_l(num_eqn), q_r(num_eqn)
    real(kind=8), intent(in) :: aux_l(num_aux, 3), aux_r(num_aux, 3)
    real(kind=8), intent(in out) :: asdq(num_eqn)

    ! Output
    real(kind=8), intent(in out) :: bmasdq(num_eqn), bpasdq(num_eqn)

    ! Locals
    integer :: normal_index, transverse_index
    real(kind=8) :: a(2)

    ! Common block
    ! density, bulk modulus, and sound speed, and impedence of medium:
    real(kind=8) :: rho, bulk, cc, zz
    common /cparam/ rho, bulk, cc, zz
    
    if (ixy == 1) then
        normal_index = 2
        transverse_index = 3
    else
        normal_index = 3
        transverse_index = 2
    endif

    a(1) = (-asdq(1) + zz * asdq(transverse_index)) / (2.d0 * zz)
    a(2) = ( asdq(1) + zz * asdq(transverse_index)) / (2.d0 * zz)

    ! Down-going flux difference bmasdq is the product of -c * wave
    bmasdq(1) = cc * a(1) * zz
    bmasdq(normal_index) = 0.d0
    bmasdq(transverse_index) = -cc * a(1)

    ! Up-going flux difference bpasdq is the product of c * wave
    bpasdq(1) = cc * a(2) * zz
    bpasdq(normal_index) = 0.d0
    bpasdq(transverse_index) = cc * a(2)
    
end subroutine rpt2_ptwise
