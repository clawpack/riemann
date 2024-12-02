! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! solve Riemann problems for the 1D Euler equations using Roe's
! approximate Riemann solver.

! waves: 3
! equations: 3

! Conserved quantities:
!       1 density
!       2 momentum
!       3 energy

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routine step1, rp is called with ql = qr = q.

    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, maux, mbc, mx
    real(kind=8), intent(in) :: ql(meqn,1-mbc:maxmx+mbc)
    real(kind=8), intent(in) :: qr(meqn,1-mbc:maxmx+mbc)
    real(kind=8), intent(in) :: auxl(maux,1-mbc:maxmx+mbc)
    real(kind=8), intent(in) :: auxr(maux,1-mbc:maxmx+mbc)

    real(kind=8), intent(out) :: s(mwaves,1-mbc:maxmx+mbc)
    real(kind=8), intent(out) :: wave(meqn, mwaves,1-mbc:maxmx+mbc)
    real(kind=8), intent(out) :: amdq(meqn,1-mbc:maxmx+mbc)
    real(kind=8), intent(out) :: apdq(meqn,1-mbc:maxmx+mbc)

    real(kind=8) delta(3)
    real(kind=8) u(1-mbc:maxmx+mbc),enth(1-mbc:maxmx+mbc)
    real(kind=8) a(1-mbc:maxmx+mbc)
    logical :: efix
    real(kind=8) :: gamma, gamma1, rhsqrtl, rhsqrtr, pl, pr, rhsq2, a1, a2, a3
    real(kind=8) :: rhoim1, pim1, cim1, s0, rho1, rhou1, en1, p1, c1, s1, sfract
    real(kind=8) :: rhoi, pi, ci, s3, rho2, rhou2, en2, p2, c2, s2, df
    integer :: i, m, mw

    common /cparam/  gamma

    data efix /.true./    !# use entropy fix for transonic rarefactions

    gamma1 = gamma - 1.d0

!     # Compute Roe-averaged quantities:

    do 20 i=2-mbc,mx+mbc
        rhsqrtl = dsqrt(qr(1,i-1))
        rhsqrtr = dsqrt(ql(1,i))
        pl = gamma1*(qr(3,i-1) - 0.5d0*(qr(2,i-1)**2)/qr(1,i-1))
        pr = gamma1*(ql(3,i) - 0.5d0*(ql(2,i)**2)/ql(1,i))
        rhsq2 = rhsqrtl + rhsqrtr
        u(i) = (qr(2,i-1)/rhsqrtl + ql(2,i)/rhsqrtr) / rhsq2
        enth(i) = (((qr(3,i-1)+pl)/rhsqrtl &
        + (ql(3,i)+pr)/rhsqrtr)) / rhsq2
        a2 = gamma1*(enth(i) - .5d0*u(i)**2)
        a(i) = dsqrt(a2)
    20 END DO


    do 30 i=2-mbc,mx+mbc
    
    !        # find a1 thru a3, the coefficients of the 3 eigenvectors:
    
        delta(1) = ql(1,i) - qr(1,i-1)
        delta(2) = ql(2,i) - qr(2,i-1)
        delta(3) = ql(3,i) - qr(3,i-1)
        a2 = gamma1/a(i)**2 * ((enth(i)-u(i)**2)*delta(1) &
        + u(i)*delta(2) - delta(3))
        a3 = (delta(2) + (a(i)-u(i))*delta(1) - a(i)*a2) / (2.d0*a(i))
        a1 = delta(1) - a2 - a3
    
    !        # Compute the waves.
    
        wave(1,1,i) = a1
        wave(2,1,i) = a1*(u(i)-a(i))
        wave(3,1,i) = a1*(enth(i) - u(i)*a(i))
        s(1,i) = u(i)-a(i)
    
        wave(1,2,i) = a2
        wave(2,2,i) = a2*u(i)
        wave(3,2,i) = a2*0.5d0*u(i)**2
        s(2,i) = u(i)
    
        wave(1,3,i) = a3
        wave(2,3,i) = a3*(u(i)+a(i))
        wave(3,3,i) = a3*(enth(i)+u(i)*a(i))
        s(3,i) = u(i)+a(i)
    30 END DO

!     # compute Godunov flux f0:
!     --------------------------


    if (efix) go to 110

!     # no entropy fix
!     ----------------

!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

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
    go to 900

!-----------------------------------------------------

    110 continue

!     # With entropy fix
!     ------------------

!    # compute flux differences amdq and apdq.
!    # First compute amdq as sum of s*wave for left going waves.
!    # Incorporate entropy fix by adding a modified fraction of wave
!    # if s should change sign.

    do 200 i=2-mbc,mx+mbc
    
    !        # check 1-wave:
    !        ---------------
    
        rhoim1 = qr(1,i-1)
        pim1 = gamma1*(qr(3,i-1) - 0.5d0*qr(2,i-1)**2 / rhoim1)
        cim1 = dsqrt(gamma*pim1/rhoim1)
        s0 = qr(2,i-1)/rhoim1 - cim1     !# u-c in left state (cell i-1)

    !        # check for fully supersonic case:
        if (s0 >= 0.d0 .AND. s(1,i) > 0.d0)  then
        !            # everything is right-going
            do 60 m=1,3
                amdq(m,i) = 0.d0
            60 END DO
            go to 200
        endif
    
        rho1 = qr(1,i-1) + wave(1,1,i)
        rhou1 = qr(2,i-1) + wave(2,1,i)
        en1 = qr(3,i-1) + wave(3,1,i)
        p1 = gamma1*(en1 - 0.5d0*rhou1**2/rho1)
        c1 = dsqrt(gamma*p1/rho1)
        s1 = rhou1/rho1 - c1  !# u-c to right of 1-wave
        if (s0 < 0.d0 .AND. s1 > 0.d0) then
        !            # transonic rarefaction in the 1-wave
            sfract = s0 * (s1-s(1,i)) / (s1-s0)
        else if (s(1,i) < 0.d0) then
        !            # 1-wave is leftgoing
            sfract = s(1,i)
        else
        !            # 1-wave is rightgoing
            sfract = 0.d0   !# this shouldn't happen since s0 < 0
        endif
        do 120 m=1,3
            amdq(m,i) = sfract*wave(m,1,i)
        120 END DO
    
    !        # check 2-wave:
    !        ---------------
    
        if (s(2,i) >= 0.d0) go to 200  !# 2-wave is rightgoing
        do 140 m=1,3
            amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
        140 END DO
    
    !        # check 3-wave:
    !        ---------------
    
        rhoi = ql(1,i)
        pi = gamma1*(ql(3,i) - 0.5d0*ql(2,i)**2 / rhoi)
        ci = dsqrt(gamma*pi/rhoi)
        s3 = ql(2,i)/rhoi + ci     !# u+c in right state  (cell i)
    
        rho2 = ql(1,i) - wave(1,3,i)
        rhou2 = ql(2,i) - wave(2,3,i)
        en2 = ql(3,i) - wave(3,3,i)
        p2 = gamma1*(en2 - 0.5d0*rhou2**2/rho2)
        c2 = dsqrt(gamma*p2/rho2)
        s2 = rhou2/rho2 + c2   !# u+c to left of 3-wave
        if (s2 < 0.d0 .AND. s3 > 0.d0) then
        !            # transonic rarefaction in the 3-wave
            sfract = s2 * (s3-s(3,i)) / (s3-s2)
        else if (s(3,i) < 0.d0) then
        !            # 3-wave is leftgoing
            sfract = s(3,i)
        else
        !            # 3-wave is rightgoing
            go to 200
        endif
    
        do 160 m=1,3
            amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
        160 END DO
    200 END DO

!     # compute the rightgoing flux differences:
!     # df = SUM s*wave   is the total flux difference and apdq = df - amdq

    do m=1,3
        do i = 2-mbc, mx+mbc
            df = 0.d0
            do mw=1,mwaves
                df = df + s(mw,i)*wave(m,mw,i)
            end do
            apdq(m,i) = df - amdq(m,i)
        end do
    end do


    900 continue
    return
    end subroutine rp1
