subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
    ! Riemann-solver for the advection equation
    !    q_t  +  u*q_x + v*q_y + w*q_z = 0
    ! where u and v are a given velocity field.
    !
    ! In advective form, with interface velocities specified in
    ! the auxiliary variable
    ! auxl(i,ma) contains auxiliary data for cells along this slice,
    !   ma=1:   u-velocity at left edge of cell 
    !   ma=2:   v-velocity at bottom edge of cell
    !   ma=3:   w-velocity at back edge of cell 
    !
    ! solve Riemann problems along one slice of data.
    ! This data is along a slice in the x-direction if ixyz=1
    !                               the y-direction if ixyz=2.
    !                               the z-direction if ixyz=3.
    !
    ! On input, ql contains the state vector at the left edge of each cell
    !           qr contains the state vector at the right edge of each cell
    !
    ! On output, wave contains the waves, s the speeds,
    ! and amdq, apdq the left-going and right-going flux differences,
    ! respectively.  Note that in this advective form, the sum of
    ! amdq and apdq is not equal to a difference of fluxes except in the
    ! case of constant velocities.
    !
    ! Note that the i'th Riemann problem has left state qr(i-1,:)
    !                                    and right state ql(i,:)
    ! From the basic clawpack routines, this routine is called with ql = qr

    implicit none

    integer, intent(in) :: ixyz, maxm, meqn, mwaves, mbc, mx, maux

    double precision, intent(in) :: ql, qr, auxl, auxr
    double precision, intent(out) :: wave, s, amdq, apdq

    dimension wave(meqn,mwaves,1-mbc:maxm+mbc)
    dimension    s(mwaves,1-mbc:maxm+mbc)
    dimension   ql(meqn,1-mbc:maxm+mbc)
    dimension   qr(meqn,1-mbc:maxm+mbc)
    dimension amdq(meqn,1-mbc:maxm+mbc)
    dimension apdq(meqn,1-mbc:maxm+mbc)
    dimension auxl(maux,1-mbc:maxm+mbc)
    dimension auxr(maux,1-mbc:maxm+mbc)

    integer :: i

    ! Set wave, speed, and flux differences:
    do i = 2-mbc, mx+mbc
        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        s(1,i) = auxl(ixyz,i)
        ! The flux difference df = s*wave all goes in the downwind direction:
        amdq(1,i) = dmin1(auxl(ixyz,i), 0.d0) * wave(1,1,i)
        apdq(1,i) = dmax1(auxl(ixyz,i), 0.d0) * wave(1,1,i)
    end do

    return
end subroutine rpn3
