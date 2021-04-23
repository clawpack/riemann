! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! solve Riemann problem for the 1D shallow water equations using the HLLE
! approximate Riemann solver.

! waves: 2
! equations: 2

! Conserved quantities:
!       1 depth (h)
!       2 momentum (hu)

    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx, maux
    double precision, dimension(meqn,1-mbc:maxmx+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxmx+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, mwaves, 1-mbc:maxmx+mbc), intent(out) :: wave
    double precision, dimension(meqn, 1-mbc:maxmx+mbc), intent(out) :: amdq, apdq
    double precision, dimension(mwaves, 1-mbc:maxmx+mbc), intent(out) :: s

    double precision :: u_l, u_r, h_l, h_r, c_l, c_r
    double precision :: hsqrt_l, hsqrt_r, u_hat, h_hat, c_hat, grav
    double precision :: h_m, hu_m
    double precision :: rho_m, rhou_m, E_m, s1, s2
    integer :: m, i, mw

    common /cparam/  grav

    do i=2-mbc,mx+mbc
        h_l = qr(1,i-1)
        h_r = ql(1,i)
        u_l = qr(2,i-1) / qr(1,i-1)
        u_r = ql(2,i  ) / ql(1,i  )
        c_l = dsqrt(grav*h_l)
        c_r = dsqrt(grav*h_r)
            
        if (h_l<1.e-14 .AND. h_r<1.e-14) then
            wave(1,1,i) =0.
            wave(2,1,i) =0.
            s(1,i) = 0.
    
            wave(1,2,i) = 0.
            wave(2,2,i) = 0.
            s(2,i) = 0.
            
        else
            ! Roe averages
            hsqrt_l = dsqrt(qr(1,i-1))
            hsqrt_r = dsqrt(ql(1,i))
            h_hat = 0.5*(h_l + h_r)
            u_hat = (hsqrt_l*u_l + hsqrt_r*u_r) / (hsqrt_l + hsqrt_r)
            c_hat = dsqrt(grav*h_hat)

            s1 = min(u_l - c_l, u_hat - c_hat)
            s2 = max(u_r + c_r, u_hat + c_hat)

            ! Middle state
            h_m = (ql(2,i) - qr(2,i-1) - s2*h_r + s1*h_l) / (s1-s2)
            hu_m = (ql(2,i)*(u_r-s2) - qr(2,i-1)*(u_l-s1) + 0.5*grav*(h_r**2-h_l**2)) / (s1-s2)
               
            wave(1,1,i) = h_m - h_l
            wave(2,1,i) = hu_m - qr(2,i-1)
            s(1,i) = s1
    
            wave(1,2,i) = h_r - h_m
            wave(2,2,i) = ql(2,i) - hu_m
            s(2,i) = s2
        end if

    end do 

    do m=1, meqn
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1, mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                endif
            end do
        end do
    end do

end subroutine rp1
