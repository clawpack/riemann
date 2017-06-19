! Riemann solver for the variable speed limit LWR traffic model
! with a tracer:

! q_t + (u_max(x,t) q(1-q))_x = 0
! p_t + u_max(x,t) (1-q) p_x = 0

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
    double precision :: f_l, f_r, s_l, s_r, f0, q_l, q_r, v_l, v_r
    integer :: i

    do i = 2-num_ghost, num_cells+num_ghost
        q_r = ql(1,i)
        q_l = qr(1,i-1)
        v_r = auxl(1,i)
        v_l = auxr(1,i-1)

        ! compute characteristic speed in each cell:
        s_l = v_l*(1.d0 - 2.d0*q_l)
        s_r = v_r*(1.d0 - 2.d0*q_r)
        s(2,i) = v_r * (1.d0 - q_r)

        ! compute flux in each cell and flux difference:
        f_l = v_l*q_l*(1.d0 - q_l)
        f_r = v_r*q_r*(1.d0 - q_r)

        ! This seems to work well even though there can be
        ! a stationary jump in q.
        wave(1,1,i) = q_r - q_l
        wave(2,1,i) = 0.d0

        ! Tracer wave and fluctuations
        wave(1,2,i) = 0.d0
        wave(2,2,i) = s(2,i)*(ql(2,i) - qr(2,i-1))
        apdq(2,i) = wave(2,2,i)
        amdq(2,i) = 0.d0 ! Traffic always moves right

        s(1,i) = 0.5d0*(s_l + s_r)

        ! Find Godunov flux in order to determine fluctuations
        if ((f_l .ge. 0.25d0*v_r) .and. (s_r .gt. 0.d0)) then
            ! left-going shock, right-going rarefaction
            f0 = 0.25d0*v_r
        elseif ((f_r .ge. 0.25d0*v_l) .and. (s_l .lt. 0.d0)) then
            ! right-going shock, left-going rarefaction
            f0 = 0.25d0*v_l
        elseif ((s_r .le. 0.d0) .and. (f_l .gt. f_r)) then
            ! left-going shock
            f0 = f_r
        elseif ((s_l .ge. 0.d0) .and. (f_r .gt. f_l)) then
            ! right-going shock
            f0 = f_l
        elseif ((s_l .le. 0.d0) .and. (s_r .ge. 0.d0)) then
            ! Transonic rarefaction
            if (v_r .le. v_l) then
                f0 = 0.25d0*v_r
            else
                f0 = 0.25d0*v_l
            endif
        elseif (f_l .le. f_r) then
            ! left-going rarefaction
            f0 = f_r
        else
            ! right-going rarefaction
            f0 = f_l
        endif
        amdq(1,i) = f0 - f_l
        apdq(1,i) = f_r - f0


    enddo
    return
end subroutine rp1
