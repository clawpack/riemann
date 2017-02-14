! Riemann solver for the LWR traffic model:

! q_t + (u_max q(1-q))_x = 0

! waves: 1
! equations: 1

! Conserved quantities:
!       1 q

! The speed limit umax should be in the common block cparam.

! See http://www.clawpack.org/riemann.html for a detailed explanation
! of the Riemann solver API.

subroutine rp1(maxm,num_eqn,num_waves,num_aux,num_ghost,num_cells, &
               ql,qr,auxl,auxr,wave,s,amdq,apdq)

    implicit none

    ! Inputs
    integer, intent(in) :: maxm, num_eqn, num_waves, num_aux, num_ghost, num_cells
    double precision, intent(in), dimension(num_eqn, 1-num_ghost:maxm+num_ghost) :: ql,qr
    double precision, intent(in), dimension(num_aux, 1-num_ghost:maxm+num_ghost) :: auxl,auxr
    ! Outputs
    double precision, intent(out) :: s(num_waves,1-num_ghost:num_cells+num_ghost)
    double precision, intent(out) :: wave(num_eqn,num_waves,1-num_ghost:num_cells+num_ghost)
    double precision, intent(out), dimension(num_eqn, 1-num_ghost:maxm+num_ghost) :: amdq,apdq
    ! Locals
    double precision :: sim1, si, flux0, fluxim1, fluxi
    integer :: i

    common /cparam/ umax
    double precision :: umax

    do i=2-num_ghost,num_cells+num_ghost
        ! Compute the wave and speed

        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        s(1,i) = umax * (1.d0 - (qr(1,i-1) + ql(1,i)))

        ! compute left-going and right-going flux differences:
        !------------------------------------------------------
        amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
        apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)

        ! entropy fix for transonic rarefactions:
        sim1 = umax*(1.d0 - 2.d0*qr(1,i-1))
        si = umax*(1.d0 - 2.d0*ql(1,i))
        if (sim1.lt.0.d0 .and. si.gt.0.d0) then
            flux0 = 0.25d0*umax
            fluxim1 = qr(1,i-1)*umax*(1.d0 - qr(1,i-1))
            fluxi   = ql(1,i)*umax*(1.d0 - ql(1,i))
            amdq(1,i) = flux0 - fluxim1
            apdq(1,i) = fluxi - flux0
        endif
    enddo

    return
end subroutine rp1
