! f-wave Riemann solver for the variable speed limit LWR traffic model

! q_t + (u_max(x,t) q(1-q))_x = 0

! Here the speed limit u_max (stored in aux(1,i) may vary with x and t.

! This solver seems to be less accurate that rp1_traffic_vc.f90, which does
! not use f-waves.  This solver also exhibits a one-cell entropy glitch
! for transonic rarefaction.  The other solver is the preferred default;
! this f-wave solver is left here mainly for educational purposes.

! waves: 1
! equations: 1
! aux fields: 1

! Conserved quantities:
!       1 q

! Auxiliary variables:
!       1 u_max

! See http://www.clawpack.org/riemann.html for a detailed explanation
! of the Riemann solver API.

subroutine rp1(maxm,num_eqn,num_waves,num_aux,num_ghost,num_cells, &
               ql,qr,auxl,auxr,fwave,s,amdq,apdq)
    implicit none
    ! Inputs
    integer, intent(in) :: maxm, num_eqn, num_waves, num_aux, num_ghost, num_cells
    double precision, intent(in), dimension(num_eqn, 1-num_ghost:maxm+num_ghost) :: ql,qr
    double precision, intent(in), dimension(num_aux, 1-num_ghost:maxm+num_ghost) :: auxl,auxr
    ! Outputs
    double precision, intent(out) :: s(num_waves,1-num_ghost:num_cells+num_ghost)
    double precision, intent(out) :: fwave(num_eqn,num_waves,1-num_ghost:num_cells+num_ghost)
    double precision, intent(out), dimension(num_eqn, 1-num_ghost:maxm+num_ghost) :: amdq,apdq
    ! Locals
    double precision :: fim1, fi, sim1, si, f0
    integer :: i

    do i = 2-num_ghost, num_cells+num_ghost
        ! Compute the wave, speeds, and flux difference

        ! compute characteristic speed in each cell:
        sim1 = auxr(1,i-1)*(1.d0 - 2.d0*qr(1,i-1))
        si   = auxl(1,i  )*(1.d0 - 2.d0*ql(1,i  ))

        ! compute flux in each cell and flux difference:
        fim1 = auxr(1,i-1)*qr(1,i-1)*(1.d0 - qr(1,i-1))
        fi   = auxl(1,i  )*ql(1,i  )*(1.d0 - ql(1,i  ))

        fwave(1,1,i) = fi - fim1

        if (sim1 .lt. 0.d0 .and. si .le. 0.d0) then
            ! left-going
            s(1,i) = sim1
            amdq(1,i) = fwave(1,1,i)
            apdq(1,i) = 0.d0
        else if (sim1 .ge. 0.d0 .and. si .gt. 0.d0) then
            ! right-going
            s(1,i) = si
            amdq(1,i) = 0.d0
            apdq(1,i) = fwave(1,1,i)
        else if (sim1 .lt. 0.d0 .and. si .gt. 0.d0) then
            ! transonic rarefaction
            s(1,i) = 0.5d0*(sim1 + si)

            ! entropy fix:  (perhaps doesn't work for all cases!!!)
            ! This assumes the flux in the transonic case should
            ! correspond to q=0.5 on the side with the smaller umax value.
            f0 = dmin1(auxr(1,i-1),auxl(1,i))*0.25d0
            ! split fwave between amdq and apdq:
            amdq(1,i) = f0 - fim1
            apdq(1,i) = fi - f0

        else
            ! transonic shock
            s(1,i) = 0.5d0*(sim1 + si)
            if (fi-fim1 .lt. 0.d0) then
                amdq(1,i) = fwave(1,1,i)
                apdq(1,i) = 0.d0
            else if (fi-fim1 .gt. 0.d0) then
                amdq(1,i) = 0.d0
                apdq(1,i) = fwave(1,1,i)
            else
                amdq(1,i) = 0.5d0 * fwave(1,1,i)
                apdq(1,i) = 0.5d0 * fwave(1,1,i)
            endif
        endif

    enddo
    return
end subroutine rp1
