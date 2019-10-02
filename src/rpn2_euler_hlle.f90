! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! HLLE solver for the Euler equations.

! waves: 2
! equations: 4

! Conserved quantities:
!       1 density
!       2 x-momentum
!       3 y-momentum
!       4 energy

    implicit none

    integer, intent(in) :: ixy, maxm, meqn, mwaves, maux, mbc, mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, mwaves, 1-mbc:maxm+mbc), intent(out) :: wave
    double precision, dimension(mwaves, 1-mbc:maxm+mbc), intent(out) :: s
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: amdq, apdq

    double precision :: u_l, u_r, c_l, c_r, u_hat, c_hat, v_l, v_r, v_hat
    double precision :: p_l, p_r, H_l, H_r, rhsqrt_l, rhsqrt_r
    double precision :: rhsq2, H_hat, gamma, gamma1, E_l, E_r
    double precision :: rho_m, rhou_m, rhov_m, E_m, s1, s2
    integer :: rho, mu, mv, E
    integer :: i, m, mw

    common /cparam/  gamma

    gamma1 = gamma - 1.d0

!     # set mu to point to  the component of the system that corresponds
!     # to momentum in the direction of this slice, mv to the orthogonal
!     # momentum:
!

    rho = 1
    E = 4
    if (ixy.eq.1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do i=2-mbc,mx+mbc
        ! Velocity
        u_l = qr(mu,i-1) / qr(rho,i-1)
        u_r = ql(mu,i  ) / ql(rho,i  )
        v_l = qr(mv,i-1) / qr(rho,i-1)
        v_r = ql(mv,i  ) / ql(rho,i  )
        ! Pressure
        p_l = gamma1 * (qr(E,i-1) - qr(rho,i-1)*(u_l**2+v_l**2)/2.d0)
        p_r = gamma1 * (ql(E,i  ) - ql(rho,i  )*(u_r**2+v_r**2)/2.d0)
        ! Enthalpy
        H_l = (qr(E,i-1) + p_l) / qr(rho,i-1)
        H_r = (ql(E,i  ) + p_l) / ql(rho,i  )
        ! Sound speed
        c_l = dsqrt(gamma1*H_l-u_l**2/2.)
        c_r = dsqrt(gamma1*H_r-u_r**2/2.)

        E_l = qr(E,i-1)
        E_r = ql(E,i  )

        rhsqrt_l = dsqrt(qr(rho,i-1))
        rhsqrt_r = dsqrt(ql(rho,i))
        rhsq2 = rhsqrt_l + rhsqrt_r
        u_hat = (qr(mu,i-1)/rhsqrt_l + ql(mu,i)/rhsqrt_r) / rhsq2
        v_hat = (qr(mv,i-1)/rhsqrt_l + ql(mv,i)/rhsqrt_r) / rhsq2
        H_hat = (((qr(E,i-1)+p_l)/rhsqrt_l + (ql(E,i)+p_r)/rhsqrt_r)) / rhsq2
        c_hat = dsqrt(gamma1*(H_hat - 0.5d0*(u_hat**2+v_hat**2)))

        ! Speeds of non-shear waves
        s1 = min(u_l - c_l, u_hat - c_hat)
        s2 = max(u_r + c_r, u_hat + c_hat)

        ! "middle" state
        rho_m = (ql(mu,i) - qr(mu,i-1) - s2*ql(rho,i) + s1*qr(rho,i-1))/(s1-s2)
        rhou_m = (ql(rho,i)*u_r**2 - qr(rho,i-1)*u_l**2 + p_r - p_l - s2*ql(mu,i) + s1*qr(mu,i-1))/(s1-s2)
        rhov_m = (ql(mv,i)*u_r - qr(mv,i-1)*u_l - s2*ql(mv,i) + s1*qr(mv,i-1))/(s1-s2)
        E_m = ( u_r*(E_r+p_r) - u_l*(E_l+p_l) -s2*E_r + s1*E_l)/(s1-s2)

        wave(rho,1,i) = rho_m - qr(rho,i-1)
        wave(mu,1,i) = rhou_m - qr(mu,i-1)
        wave(mv,1,i) = rhov_m - qr(mv,i-1)
        wave(E,1,i) = E_m - qr(E,i-1)
        s(1,i) = s1
    
        wave(rho,2,i) = ql(rho,i) - rho_m
        wave(mu,2,i) = ql(mu,i) - rhou_m
        wave(mv,2,i) = ql(mv,i) - rhov_m
        wave(E,2,i) = ql(E,i) - E_m
        s(2,i) = s2
    end do


    do m=1,meqn
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1,mwaves
                if (s(mw,i) .lt. 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
            end do
        end do
    end do

end subroutine rpn2
