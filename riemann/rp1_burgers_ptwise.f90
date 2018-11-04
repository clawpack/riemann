! =============================================================================
subroutine rp1_ptwise(num_eqn, num_aux, num_waves, q_l, q_r, aux_l, aux_r,    &
wave, s, amdq, apdq)
! =============================================================================
!
! Riemann problems for the 1D Burgers' equation with entropy fix for 
! transonic rarefaction. See "Finite Volume Method for Hyperbolic Problems",
! R. J. LeVeque.

! waves: 1
! equations: 1

! Conserved quantities:
!       1 q
    
    implicit none

    ! Input Arguments
    integer, intent(in) :: num_eqn, num_aux, num_waves
    real(kind=8), intent(in out) :: q_l(num_eqn), q_r(num_eqn)
    real(kind=8), intent(in out) :: aux_l(num_aux), aux_r(num_aux)

    ! Output arguments
    real(kind=8), intent(in out) :: wave(num_eqn, num_waves)
    real(kind=8), intent(in out) :: s(num_waves)
    real(kind=8), intent(in out) :: apdq(num_eqn), amdq(num_eqn)

    ! Local
    logical :: efix
 
    efix = .true.   !# Compute correct flux for transonic rarefactions

    wave(1,1) = q_r(1) - q_l(1)
    s(1) = 0.5d0 * (q_l(1) + q_r(1))

    amdq(1) = dmin1(s(1), 0.d0) * wave(1,1)
    apdq(1) = dmax1(s(1), 0.d0) * wave(1,1)

    if (efix) then
        if (q_r(1).gt.0.d0 .and. q_l(1).lt.0.d0) then
            amdq(1) = - 1.d0/2.d0 * q_l(1)**2
            apdq(1) =   1.d0/2.d0 * q_r(1)**2
        endif
    endif

    end subroutine rp1_ptwise
