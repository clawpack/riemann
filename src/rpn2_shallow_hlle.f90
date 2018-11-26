! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! HLLE solver for the 2D shallow water equations.

! waves: 2
! equations: 3

! Conserved quantities:
!       1 depth
!       2 x-momentum
!       3 y-momentum

    implicit none

    integer, intent(in) :: ixy, maxm, meqn, mwaves, maux, mbc, mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, mwaves, 1-mbc:maxm+mbc), intent(out) :: wave
    double precision, dimension(mwaves, 1-mbc:maxm+mbc), intent(out) :: s
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: amdq, apdq

    double precision :: u_l, u_r, c_l, c_r, u_hat, c_hat, v_l, v_r, v_hat, h_hat
    double precision :: h_l, h_r, hsqrt_l, hsqrt_r, hsq2
    double precision :: grav
    double precision :: h_m, hu_m, hv_m, s1, s2
    integer :: depth, mu, mv
    integer :: i, m, mw

    common /cparam/  grav

!     # set mu to point to  the component of the system that corresponds
!     # to momentum in the direction of this slice, mv to the orthogonal
!     # momentum:
!

    depth = 1
    if (ixy.eq.1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do i=2-mbc,mx+mbc
        h_l = qr(depth,i-1)
        h_r = ql(depth,i)
        ! Velocity
        u_l = qr(mu,i-1) / qr(depth,i-1)
        u_r = ql(mu,i  ) / ql(depth,i  )
        v_l = qr(mv,i-1) / qr(depth,i-1)
        v_r = ql(mv,i  ) / ql(depth,i  )
        ! Sound speed
        c_l = dsqrt(grav*h_l)
        c_r = dsqrt(grav*h_r)

        hsqrt_l = dsqrt(qr(depth,i-1))
        hsqrt_r = dsqrt(ql(depth,i))
        hsq2 = hsqrt_l + hsqrt_r
        h_hat = 0.5*(h_l + h_r)
        u_hat = (u_l*hsqrt_l + u_r*hsqrt_r) / hsq2
        v_hat = (v_l*hsqrt_l + v_r*hsqrt_r) / hsq2
        c_hat = dsqrt(grav*h_hat)

        ! Speeds of non-shear waves
        s1 = min(u_l - c_l, u_hat - c_hat)
        s2 = max(u_r + c_r, u_hat + c_hat)

        ! middle" state
        h_m = (ql(mu,i) - qr(mu,i-1) - s2*ql(depth,i) + s1*qr(depth,i-1))/(s1-s2)
        hu_m = (ql(mu,i)*(u_r-s2) - qr(mu,i-1)*(u_l-s1) + 0.5*grav*(h_r**2 - h_l**2) ) / (s1-s2)
        hv_m = (ql(mv,i)*u_r - qr(mv,i-1)*u_l - s2*ql(mv,i) + s1*qr(mv,i-1))/(s1-s2)

        wave(depth,1,i) = h_m - h_l
        wave(mu,1,i) = hu_m - qr(mu,i-1)
        wave(mv,1,i) = hv_m - qr(mv,i-1)
        s(1,i) = s1

        wave(depth,2,i) = h_r - h_m
        wave(mu,2,i) = ql(mu,i) - hu_m
        wave(mv,2,i) = ql(mv,i) - hv_m
        s(2,i) = s2
    end do


    do m=1, meqn
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1, mwaves
                if (s(mw,i) .lt. 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
            end do
        end do
    end do

end subroutine rpn2
