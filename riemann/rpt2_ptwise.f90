! ==============================================================================
!  This Riemann solver provides a harness for the point-wise evaluation of 
!  Riemann solutions and requires a point-wise definition of the Riemann solver
!  in a subroutine of the form "rpt2_ptwise".
! ==============================================================================

subroutine rpt2(ixy, imp, maxm, meqn, mwaves, maux, mbc, mx, ql, qr,          &
                aux1, aux2, aux3, asdq, bmasdq, bpasdq)

    implicit none

    ! Input
    integer, intent(in) :: ixy, imp, maxm, meqn, mwaves, maux, mbc, mx
    real(kind=8), intent(in out) :: ql(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: qr(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: aux1(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: aux2(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: aux3(maux, 1-mbc:maxm+mbc)

    ! Output
    real(kind=8), intent(in out) :: asdq(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: bmasdq(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: bpasdq(meqn, 1-mbc:maxm+mbc)

    ! Locals
    integer :: i
    real(kind=8) :: aux_r(maux, 3)
    real(kind=8) :: aux_l(maux, 3)

    do i = 2 - mbc, mx + mbc

        aux_r(:, 1) = aux1(:, i)
        aux_r(:, 2) = aux2(:, i)
        aux_r(:, 3) = aux3(:, i)
        aux_l(:, 1) = aux1(:, i - 1)
        aux_l(:, 2) = aux2(:, i - 1)
        aux_l(:, 3) = aux3(:, i - 1)
        call rpt2_ptwise(ixy, imp, meqn, maux, mwaves, qr(:, i - 1),    &
                                                       ql(:, i),        &
                                                       aux_l,           &
                                                       aux_r,           &
                                                       asdq(:, i),      &
                                                       bmasdq(:, i),    &
                                                       bpasdq(:, i))

    end do

end subroutine rpt2
