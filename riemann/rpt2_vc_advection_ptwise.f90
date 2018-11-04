!     # Riemann solver in the transverse direction for the advection equation.

! =====================================================
subroutine rpt2_ptwise(ixy, imp, num_eqn, num_aux, num_waves, q_l, q_r,      &
aux_l, aux_r, asdq, bmasdq, bpasdq)
! =====================================================

    implicit none

    ! Input
    integer, intent(in) :: num_eqn, num_aux, num_waves, ixy, imp
    real(kind=8), intent(in) :: q_l(num_eqn), q_r(num_eqn)
    real(kind=8), intent(in) :: aux_l(num_aux, 3), aux_r(num_aux, 3)
    real(kind=8), intent(in out) :: asdq(num_eqn)

    ! Output
    real(kind=8), intent(in out) :: bmasdq(num_eqn), bpasdq(num_eqn)

    ! Local
    real(kind=8) :: aux(num_aux, 3)
    integer :: kv

    kv = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1

    !#  =  i-1 for amdq,  i for apdq
    if (ixy == 1) then
        aux = aux_l
    else
        aux = aux_r
    endif

    bmasdq(1) = dmin1(aux(kv,2), 0.d0) * asdq(1)
    bpasdq(1) = dmax1(aux(kv,3), 0.d0) * asdq(1)

    return
    end subroutine rpt2_ptwise
