! =====================================================
subroutine rpt2_ptwise(ixy, imp, num_eqn, num_aux, num_waves, q_l, q_r,      &
aux_l, aux_r, asdq, bmasdq, bpasdq)
! =====================================================

!     # Riemann solver in the transverse direction for the scalar equation

!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)

    implicit none

    ! Input
    integer, intent(in) :: num_eqn, num_aux, num_waves, ixy, imp
    real(kind=8), intent(in) :: q_l(num_eqn), q_r(num_eqn)
    real(kind=8), intent(in) :: aux_l(num_aux), aux_r(num_aux)
    real(kind=8), intent(in out) :: asdq(num_eqn)

    ! Output
    real(kind=8), intent(in out) :: bmasdq(num_eqn), bpasdq(num_eqn)

    ! Local
    real(kind=8) :: stran

    ! Common block
    real(kind=8) :: u, v
    common /cparam/ u, v

!    # transverse wave speeds have been computed in rpn2
!    # maux=0 and aux arrays are unused in this example.

    if (ixy == 1) then
        stran = v
    else
        stran = u
    endif

    bmasdq(1) = dmin1(stran, 0.d0) * asdq(1)
    bpasdq(1) = dmax1(stran, 0.d0) * asdq(1)

    end subroutine rpt2_ptwise
