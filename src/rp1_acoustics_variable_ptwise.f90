! =====================================================
subroutine rp1_ptwise(num_eqn, num_aux, num_waves, q_l, q_r, aux_l, aux_r,    &
                       wave, s, amdq, apdq)
! =====================================================

!     # Pointwise Riemann solver for the acoustics equations in 1d, with
!     #  variable coefficients (heterogeneous media)

! waves:     2
! equations: 2
! aux fields: 2

! Conserved quantities:
!       1 pressure
!       2 velocity

! Auxiliary variables:
!       1 impedance
!       2 sound_speed

!     # auxl(1,i) should contain the impedance Z in cell i
!     # auxl(2,i) should contain the sound speed c in cell i

!     # On input, q_l contains the state vector on the left
!     #           q_r contains the state vector on the right

!     # On output, wave contains the waves,
!     #            s the speeds,
!     #
!     #            amdq = A^- Delta q,
!     #            apdq = A^+ Delta q,
!     #                   the decomposition of the flux difference
!     #                       f(q_l) - f(q_r)
!     #                   into leftgoing and rightgoing parts respectively.
!     #


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
    real(kind=8) :: delta(2), a(2), zi, zim

!   # Split the jump in q into waves:

!   # Find a(1) and a(2), the coefficients of the 2 eigenvectors
    delta = q_r - q_l

!   # Impedances
    zi = aux_r(1)
    zim = aux_l(1)

    a(1) = (-delta(1) + zi*delta(2)) / (zim + zi)
    a(2) =  (delta(1) + zim*delta(2)) / (zim + zi)
    
!   # Compute the waves:
    
    wave(1,1) = -a(1)*zim
    wave(2,1) = a(1)
    s(1) = -aux_l(2)
    
    wave(1,2) = a(2)*zi
    wave(2,2) = a(2)
    s(2) = aux_r(2)

!   # Compute the leftgoing and rightgoing fluctuations:
!   # Note s(1,i) < 0 and s(2,i) > 0.

    amdq = s(1)*wave(:,1)
    apdq = s(2)*wave(:,2)

    end subroutine rp1_ptwise
