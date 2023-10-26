! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! Roe-solver for the 2D shallow water equations
! solve Riemann problems along one slice of data.

! waves: 3
! equations: 3

! Conserved quantities:
!       1 depth
!       2 x_momentum
!       3 y_momentum

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.

! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
!   f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.
! With the Roe solver we have
!    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated
! into the flux differences.

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr


! This Riemann solver differs from the original clawpack Riemann solver
! for the interleaved indices

    implicit double precision (a-h,o-z)

    dimension   ql(meqn,           1-mbc:maxm+mbc)
    dimension   qr(meqn,           1-mbc:maxm+mbc)
    dimension    s(mwaves,         1-mbc:maxm+mbc)
    dimension wave(meqn,   mwaves, 1-mbc:maxm+mbc)
    dimension  apdq(meqn,          1-mbc:maxm+mbc)
    dimension  amdq(meqn,          1-mbc:maxm+mbc)

!   # Gravity constant set in the shallow1D.py or setprob.f90 file
    common /cparam/ grav

!   # Roe averages quantities of each interface
    parameter (maxm2 = 1800)
    double precision u(-6:maxm2+7),v(-6:maxm2+7),a(-6:maxm2+7), &
                     h(-6:maxm2+7)


!   local arrays
!   ------------
    dimension delta(3)
    logical :: efix

    data efix /.true./    !# Use entropy fix for transonic rarefactions

!   # Set mu to point to  the component of the system that corresponds
!   # to momentum in the direction of this slice, mv to the orthogonal
!   # momentum:

    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif


!   # Note that notation for u and v reflects assumption that the
!   # Riemann problems are in the x-direction with u in the normal
!   # direciton and v in the orthogonal direcion, but with the above
!   # definitions of mu and mv the routine also works with ixy=2
!   # and returns, for example, f0 as the Godunov flux g0 for the
!   # Riemann problems u_t + g(u)_y = 0 in the y-direction.


!   # Compute the Roe-averaged variables needed in the Roe solver.

    do 10 i = 2-mbc, mx+mbc
        h(i) = (qr(1,i-1)+ql(1,i))*0.50d0
        hsqrtl = dsqrt(qr(1,i-1))
        hsqrtr = dsqrt(ql(1,i))
        hsq2 = hsqrtl + hsqrtr
        u(i) = (qr(mu,i-1)/hsqrtl + ql(mu,i)/hsqrtr) / hsq2
        v(i) = (qr(mv,i-1)/hsqrtl + ql(mv,i)/hsqrtr) / hsq2
        a(i) = dsqrt(grav*h(i))
    10 END DO


!   # now split the jump in q at each interface into waves

!   # find a1 thru a3, the coefficients of the 3 eigenvectors:
    do 20 i = 2-mbc, mx+mbc
        delta(1) = ql(1,i) - qr(1,i-1)
        delta(2) = ql(mu,i) - qr(mu,i-1)
        delta(3) = ql(mv,i) - qr(mv,i-1)
        a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
        a2 = -v(i)*delta(1) + delta(3)
        a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))
    
!      # Compute the waves.
    
        wave(1,1,i) = a1
        wave(mu,1,i) = a1*(u(i)-a(i))
        wave(mv,1,i) = a1*v(i)
        s(1,i) = u(i)-a(i)
    
        wave(1,2,i) = 0.0d0
        wave(mu,2,i) = 0.0d0
        wave(mv,2,i) = a2
        s(2,i) = u(i)
    
        wave(1,3,i) = a3
        wave(mu,3,i) = a3*(u(i)+a(i))
        wave(mv,3,i) = a3*v(i)
        s(3,i) = u(i)+a(i)
    20 END DO


!    # compute flux differences amdq and apdq.
!    ---------------------------------------

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
!       # check 1-wave
        him1 = qr(1,i-1)
        s0 =  qr(mu,i-1)/him1 - dsqrt(grav*him1)
!       # check for fully supersonic case :
        if (s0 > 0.0d0 .AND. s(1,i) > 0.0d0) then
            do 60 m=1,3
                amdq(m,i)=0.0d0
            60 END DO
            goto 200
        endif
    
        h1 = qr(1,i-1)+wave(1,1,i)
        hu1= qr(mu,i-1)+wave(mu,1,i)
        s1 = hu1/h1 - dsqrt(grav*h1) !speed just to right of 1-wave
        if (s0 < 0.0d0 .AND. s1 > 0.0d0) then
!           # transonic rarefaction in 1-wave
            sfract = s0*((s1-s(1,i))/(s1-s0))
        else if (s(1,i) < 0.0d0) then
!           # 1-wave is leftgoing
            sfract = s(1,i)
        else
!           # 1-wave is rightgoing
            sfract = 0.0d0
        endif
        do 120 m=1,3
            amdq(m,i) = sfract*wave(m,1,i)
        120 END DO

!       # check 2-wave
        if (s(2,i) > 0.0d0) then
!           # 2 and 3 waves are right-going
            go to 200
        endif

        do 140 m=1,3
            amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
        140 END DO
    
!       check 3-wave
    
        hi = ql(1,i)
        s03 = ql(mu,i)/hi + dsqrt(grav*hi)
        h3=ql(1,i)-wave(1,3,i)
        hu3=ql(mu,i)-wave(mu,3,i)
        s3=hu3/h3 + dsqrt(grav*h3)
        if (s3 < 0.0d0 .AND. s03 > 0.0d0) then
!           # transonic rarefaction in 3-wave
            sfract = s3*((s03-s(3,i))/(s03-s3))
        else if (s(3,i) < 0.0d0) then
!           # 3-wave is leftgoing
            sfract = s(3,i)
        else
!           # 3-wave is rightgoing
            goto 200
        endif
        do 160 m=1,3
            amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
        160 END DO
    200 END DO

!   compute rightgoing flux differences :

    do m=1,3
        do i = 2-mbc,mx+mbc
            df = 0.0d0
            do mw=1,mwaves
                df = df + s(mw,i)*wave(m,mw,i)
            end do
            apdq(m,i)=df-amdq(m,i)
        end do
    end do


    900 continue
    return
    end subroutine rpn2


