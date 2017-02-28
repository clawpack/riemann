! =====================================================
subroutine rpn2_ptwise(ixy, num_eqn, num_aux, num_waves, q_l, q_r,            &
aux_l, aux_r, wave, s, amdq, apdq)
! =====================================================

! Riemann solver for the sample scalar equation
!  q_t + u*q_x + v*q_y = 0

! waves: 1
! equations: 1

! Conserved quantities:
!       1 q

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.
! On output, wave contains the waves,
!            s the speeds,
!
!            amdq = A^- Delta q,
!            apdq = A^+ Delta q,
!                   the decomposition of the flux difference
!                       f(qr(i-1)) - f(ql(i))
!                   into leftgoing and rightgoing parts respectively.
!

! maux=0 and aux arrays are unused in this example.

    implicit none

    ! Input
    integer, intent(in) :: num_eqn, num_aux, num_waves, ixy
    real(kind=8), intent(in) :: q_l(num_eqn), q_r(num_eqn)
    real(kind=8), intent(in) :: aux_l(num_aux), aux_r(num_aux)

    ! Output
    real(kind=8), intent(in out) :: wave(num_eqn, num_waves)
    real(kind=8), intent(in out) :: s(num_waves)
    real(kind=8), intent(in out) :: apdq(num_eqn), amdq(num_eqn)

    ! Common block
    real(kind=8) :: u, v
    common /cparam/ u, v

    wave(1,1) = q_r(1) - q_l(1)
    if (ixy == 1) then
        s(1) = u
    else
        s(1) = v
    endif
    
    ! Flux differences
    amdq(1) = dmin1(s(1), 0.d0) * wave(1,1)
    apdq(1) = dmax1(s(1), 0.d0) * wave(1,1)

    end subroutine rpn2_ptwise
