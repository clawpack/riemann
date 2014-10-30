subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! Riemann-solver for the 3d Burgers' equation

!    q_t  +  u*(.5*q^2)_x + v*(.5*q^2)_y + w*(.5*q^2)_z = 0

! where u,v,w are a given scalars, stored in the common block cparam.
!
! waves: 1
! equations: 1
!
! Conserved quantities:
!       1 q
!
! Solve Riemann problems along one slice of data.
! This data is along a slice in the x-direction if ixyz=1
!                               the y-direction if ixyz=2.
!                               the z-direction if ixyz=3.
!
! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
!
! On output, wave contains the waves, s the speeds,
! and amdq, apdq the left-going and right-going flux differences,
! respectively.  
!
! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr

    implicit none

    integer, intent(in) :: ixyz, maxm, meqn, mwaves, mbc, mx, maux

    double precision, intent(in)  :: ql, qr, auxl, auxr
    double precision, intent(out) :: wave, s, amdq, apdq
    dimension wave(meqn,mwaves,1-mbc:maxm+mbc)
    dimension    s(mwaves,1-mbc:maxm+mbc)
    dimension   ql(meqn,1-mbc:maxm+mbc)
    dimension   qr(meqn,1-mbc:maxm+mbc)
    dimension amdq(meqn,1-mbc:maxm+mbc)
    dimension apdq(meqn,1-mbc:maxm+mbc)
    logical, parameter :: efix = .true.

    double precision :: coeff
    common /cparam/ coeff(3)

    integer :: i

    ! Set wave, speed, and flux differences:
    do i = 2-mbc, mx+mbc
        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        s(1,i) = coeff(ixyz) * 0.5d0*(qr(1,i-1) + ql(1,i))

        ! The flux difference df = s*wave all goes in the downwind direction:
        amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
        apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)

        if (efix) then
            ! entropy fix for transonic rarefactions:
            if (qr(1,i-1).lt.0.d0 .and. ql(1,i).gt.0.d0) then
                amdq(1,i) = - coeff(ixyz) * 0.5d0*qr(1,i-1)**2
                apdq(1,i) =   coeff(ixyz) * 0.5d0*ql(1,i)**2
            endif
        endif
    end do

    return

end subroutine rpn3

