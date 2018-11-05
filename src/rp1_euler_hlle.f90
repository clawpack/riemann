! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! solve Riemann problems for the 1D Euler equations using the HLLE
! approximate Riemann solver.

! waves: 2
! equations: 3

! Conserved quantities:
!       1 density
!       2 momentum
!       3 energy

    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx, maux
    double precision, dimension(meqn,1-mbc:maxmx+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxmx+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, mwaves, 1-mbc:maxmx+mbc), intent(out) :: wave
    double precision, dimension(meqn, 1-mbc:maxmx+mbc), intent(out) :: amdq, apdq
    double precision, dimension(mwaves, 1-mbc:maxmx+mbc), intent(out) :: s

    double precision :: u_l, u_r, c_l, c_r, u_hat, c_hat
    double precision :: p_l, p_r, H_l, H_r, rhsqrt_l, rhsqrt_r
    double precision :: rhsq2, H_hat, gamma, gamma1, E_l, E_r
    double precision :: rho_m, rhou_m, E_m, s1, s2
    integer :: m, i, mw

    common /cparam/  gamma

    gamma1 = gamma - 1.d0

    do i=2-mbc,mx+mbc
        ! Velocity
        u_l = qr(2,i-1) / qr(1,i-1)
        u_r = ql(2,i  ) / ql(1,i  )
        ! Pressure
        p_l = gamma1 * (qr(3,i-1) - qr(1,i-1)*u_l**2/2.d0)
        p_r = gamma1 * (ql(3,i  ) - ql(1,i  )*u_r**2/2.d0)
        ! Enthalpy
        H_l = (qr(3,i-1) + p_l) / qr(1,i-1)
        H_r = (ql(3,i  ) + p_l) / ql(1,i  )
        ! Sound speed
        c_l = dsqrt(gamma1*H_l-u_l**2/2.)
        c_r = dsqrt(gamma1*H_r-u_r**2/2.)

        E_l = qr(3,i-1)
        E_r = ql(3,i  )


        rhsqrt_l = dsqrt(qr(1,i-1))
        rhsqrt_r = dsqrt(ql(1,i))
        rhsq2 = rhsqrt_l + rhsqrt_r
        u_hat = (qr(2,i-1)/rhsqrt_l + ql(2,i)/rhsqrt_r) / rhsq2
        H_hat = (((qr(3,i-1)+p_l)/rhsqrt_l + (ql(3,i)+p_r)/rhsqrt_r)) / rhsq2
        c_hat = dsqrt(gamma1*(H_hat - .5d0*u_hat**2))

        s1 = min(u_l - c_l, u_hat - c_hat)
        s2 = max(u_r + c_r, u_hat + c_hat)
        rho_m = (ql(2,i) - qr(2,i-1) - s2*ql(1,i) + s1*qr(1,i-1))/(s1-s2)
        rhou_m = (ql(1,i)*u_r**2 - qr(1,i-1)*u_l**2 + p_r - p_l - s2*ql(2,i) + s1*qr(2,i-1))/(s1-s2)
        E_m = ( u_r*(E_r+p_r) - u_l*(E_l+p_l) -s2*E_r + s1*E_l)/(s1-s2)

        wave(1,1,i) = rho_m - qr(1,i-1)
        wave(2,1,i) = rhou_m - qr(2,i-1)
        wave(3,1,i) = E_m - qr(3,i-1)
        s(1,i) = s1
    
        wave(1,2,i) = ql(1,i) - rho_m
        wave(2,2,i) = ql(2,i) - rhou_m
        wave(3,2,i) = ql(3,i) - E_m
        s(2,i) = s2

    end do 

    do m=1,3
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                endif
            end do
        end do
    end do

end subroutine rp1
