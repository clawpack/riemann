! =====================================================
subroutine rpn2_ptwise(ixy, num_eqn, num_aux, num_waves, q_l, q_r,            &
aux_l, aux_r, wave, s, amdq, apdq)
! =====================================================
! Riemann-solver for the advection equation
!    q_t  +  u*q_x + v*q_y = 0
! where u and v are a given velocity field.

! waves: 1
! equations: 1
! aux fields: 2

! Conserved quantities:
!       1 q

! Auxiliary variables:
!         1  x_velocity
!         2  y_velocity

! solve Riemann problem along one slice of data.
! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! On output, wave contains the waves, s the speeds,
! and amdq, apdq the left-going and right-going flux differences,
! respectively.  Note that in this advective form, the sum of
! amdq and apdq is not equal to a difference of fluxes except in the
! case of constant velocities.

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


!     # Set wave, speed, and flux differences:
!     ------------------------------------------

        wave(1,1) = q_r(1) - q_l(1)
        s(1) = aux_r(ixy)
    !        # The flux difference df = s*wave  all goes in the downwind direction:
        amdq(1) = dmin1(aux_r(ixy), 0.d0) * wave(1,1)
        apdq(1) = dmax1(aux_r(ixy), 0.d0) * wave(1,1)

    return
    end subroutine rpn2_ptwise
