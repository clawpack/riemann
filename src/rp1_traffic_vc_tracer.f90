! Riemann solver for the variable speed limit LWR traffic model
! with a tracer:

! q_t + (u_max q(1-q))_x = 0
! p_t + u_max(1-q) p_x = 0

! Here the speed limit u_max (stored in aux(1,i) may vary with x and t.

! waves: 2
! equations: 2
! aux fields: 1

! Conserved quantities:
!       1 q
!       2 p

! Auxiliary variables:
!       1 u_max

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
    double precision :: fim1, fi, sim1, si, f0, fluxdiff
    integer :: i

    do i = 2-num_ghost, num_cells+num_ghost
        ! Compute the wave, speeds, and flux difference

        ! compute flux in each cell and flux difference:
        fim1 = auxr(1,i-1)*qr(1,i-1)*(1.d0 - qr(1,i-1))
        fi   = auxl(1,i  )*ql(1,i  )*(1.d0 - ql(1,i  ))
        fluxdiff = fi - fim1
        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        wave(2,1,i) = 0.d0
        wave(1,2,i) = 0.d0
        wave(2,2,i) = ql(2,i) - qr(2,i-1)

        ! compute characteristic speed in each cell:
        sim1 = auxr(1,i-1)*(1.d0 - 2.d0*qr(1,i-1))
        si   = auxl(1,i  )*(1.d0 - 2.d0*ql(1,i  ))

        s(2,i) = auxl(1,i) * (1.d0 - ql(1,i))
        apdq(2,i) = s(2,i) * wave(2,2,i)
        amdq(2,i) = 0.d0 ! Traffic always moves right

        if (sim1 .lt. 0.d0 .and. si .le. 0.d0) then
            ! left-going
            s(1,i) = sim1
            amdq(1,i) = fluxdiff
            apdq(1,i) = 0.d0
        else if (sim1 .ge. 0.d0 .and. si .gt. 0.d0) then
            ! right-going
            s(1,i) = si
            amdq(1,i) = 0.d0
            apdq(1,i) = fluxdiff
        else if (sim1 .lt. 0.d0 .and. si .gt. 0.d0) then
            ! transonic rarefaction
            s(1,i) = 0.5d0*(sim1 + si)

            ! entropy fix:  (perhaps doesn't work for all cases!!!)
            ! This assumes the flux in the transonic case should
            ! correspond to q=0.5 on the side with the smaller umax value.
            f0 = dmin1(auxr(1,i-1),auxl(1,i))*0.25d0
            ! split fluxdiff between amdq and apdq:
            amdq(1,i) = f0 - fim1
            apdq(1,i) = fi - f0

        else
            ! transonic shock
            s(1,i) = 0.5d0*(sim1 + si)
            if (s(1,i) .lt. 0.d0) then 
                amdq(1,i) = fluxdiff
                apdq(1,i) = 0.d0
            else if (s(1,i) .gt. 0.d0) then 
                amdq(1,i) = 0.d0
                apdq(1,i) = fluxdiff
            else
                amdq(1,i) = 0.5d0 * fluxdiff 
                apdq(1,i) = 0.5d0 * fluxdiff
            endif
        endif

    enddo
    return
    end
