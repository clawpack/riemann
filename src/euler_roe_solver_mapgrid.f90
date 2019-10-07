    subroutine roe_solver(ixy,uv,enth,delta,wave_local,s_local,info)
    implicit none

    double precision :: uv(2),enth, delta(4)
    double precision :: wave_local(3,4), s_local(3)
    double precision :: euv, c2, c, u2v2, u,v
    double precision :: a1, a2, a3, a4, a5
    double precision :: gamma, gamma1
    integer :: m, p, mu, mv, ixy, i, j, k, info

    common /cparam/ gamma, gamma1


    info = 0

    gamma1 = gamma - 1.d0

!     # These are used even in the mapped case, but ixy is
!     # always set to 1;  rotation matrix rot does switches
!     # coordinates for us.
    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    u = uv(mu-1)
    v = uv(mv-1)

    u2v2 = u*u + v*v
    c2 = gamma1*(enth - 0.5d0*u2v2)
    if (c2 < 1e-13) then
        info = 1
        write(6,*) 'Roe solver : '
        write(6,'(A,E16.8, E16.8)') 'c2 .lt. 0; ', c2
        return
    !         stop
    endif
    c = sqrt(c2) !! sound speed
    euv = enth - u2v2

    a2 = (gamma1/c2) * (euv*delta(1) + u*delta(mu) + v*delta(mv) &
    - delta(4))
    a3 = delta(mv) - v*delta(1)
    a4 = (delta(mu) + (c-u)*delta(1) - c*a2) / (2.d0*c)
    a1 = delta(1) - a2 - a4

!     # 1-wave (rarefaction or shock)
    wave_local(1,1)  = a1
    wave_local(1,mu) = a1*(u-c)
    wave_local(1,mv) = a1*v
    wave_local(1,4)  = a1*(enth - u*c)
    s_local(1) = u - c

!     # contact discontinuity
    wave_local(2,1)  = a2
    wave_local(2,mu) = a2*u
    wave_local(2,mv) = a2*v  + a3
    wave_local(2,4)  = a2*0.5d0*u2v2 + a3*v
    s_local(2) = u

!     # 3-wave (rarefaction or shock)
    wave_local(3,1)  = a4
    wave_local(3,mu) = a4*(u+c)
    wave_local(3,mv) = a4*v
    wave_local(3,4)  = a4*(enth + u*c)
    s_local(3) = u + c

    end subroutine roe_solver
