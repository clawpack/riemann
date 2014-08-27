subroutine rpn2(ixy,maxm, meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)

! Roe-solver for the Euler equations with mapped grids

! waves: 3
! equations: 4

! Conserved quantities:
!       1 density
!       2 x_momentum
!       3 y_momentum
!       4 energy

! Need to fill in the auxiliary variables

! Solve Riemann problems along one slice of data.

    implicit none

    integer :: ixy, maxm, meqn, mwaves, mbc, mx, my, maux
    double precision ::    ql(meqn,   1-mbc:maxm+mbc)
    double precision ::    qr(meqn,   1-mbc:maxm+mbc)
    double precision ::     s(mwaves, 1-mbc:maxm+mbc)
    double precision ::  wave(meqn,   mwaves, 1-mbc:maxm+mbc)
    double precision ::  amdq(meqn,   1-mbc:maxm+mbc)
    double precision ::  apdq(meqn,   1-mbc:maxm+mbc)
    double precision ::  auxl(maux,      1-mbc:maxm+mbc)
    double precision ::  auxr(maux,      1-mbc:maxm+mbc)

    double precision :: gamma, gamma1
    double precision :: ql_state(4), qr_state(4)
    double precision :: rhol, ul, el, cl, pl, vl
    double precision :: rhor, ur, er, cr, pr, vr
    double precision :: rhsqrtl, rhsqrtr, rhsq2
    double precision :: u, enth, delta(4), rho, v
    double precision :: wave_local(3,4), s_local(3),uv(2)
    double precision :: speeds(2,3), u2v2l, u2v2r
    integer :: i, m, mw, mu, mv, ixy1, info
    logical :: efix

    integer :: mcapa,locrot, locarea
    double precision :: rot(4), area

    common /param/  gamma

    data efix /.true./

    gamma1 = gamma - 1.d0

    call get_aux_locations_n(ixy,mcapa,locrot,locarea)

    do i = 2-mbc,mx+mbc

        rot(1) = auxl(locrot,i)
        rot(2) = auxl(locrot+1,i)
        call compute_tangent(rot)

        do m = 1,meqn
            ql_state(m) = qr(m,i-1)
            qr_state(m) = ql(m,i)
        enddo
        call rotate2(rot,ql_state(2))
        call rotate2(rot,qr_state(2))

        rhol = ql_state(1)
        rhor = qr_state(1)

        ul = ql_state(2)/rhol
        ur = qr_state(2)/rhor

        vl = ql_state(3)/rhol
        vr = qr_state(3)/rhor

        el = ql_state(4)
        er = qr_state(4)

        u2v2l = ul*ul + vl*vl
        u2v2r = ur*ur + vr*vr
        pl = gamma1*(el - 0.5d0*rhol*u2v2l)
        pr = gamma1*(er - 0.5d0*rhor*u2v2r)

    !        # Get Roe averaged values
        rhsqrtl = sqrt(rhol)
        rhsqrtr = sqrt(rhor)
        rhsq2 = rhsqrtl + rhsqrtr

        uv(1) = (ul*rhsqrtl + ur*rhsqrtr) / rhsq2
        uv(2) = (vl*rhsqrtl + vr*rhsqrtr) / rhsq2
        enth = ((el + pl)/rhsqrtl + (er + pr)/rhsqrtr)/rhsq2

        do m = 1,meqn
            delta(m) = qr_state(m) - ql_state(m)
        enddo

        ixy1 = 1
        call roe_solver(ixy1,uv,enth,delta,wave_local,s_local,info)

        if (info /= 0) then
            write(6,*) 'ixy = ', ixy
            write(6,*) 'enth = ', enth
            write(6,*) 'Called from rpn2 '
            write(6,*) ' '
            do m = 1,meqn
                write(6,'(2E24.16)') ql_state(m), qr_state(m)
            enddo
            write(6,*) ' '
            stop
        endif

        do mw = 1,mwaves
            speeds(1,mw) = min(s_local(mw),0.d0)
            speeds(2,mw) = max(s_local(mw),0.d0)
        enddo

        if (efix) then
        !           # This modifies the speeds, but we will still have
        !           # s(mw) = speeds(mw,1) + speeds(mw,2)
            cl = sqrt(gamma*pl/rhol)
            cr = sqrt(gamma*pr/rhor)
            call apply_entropy_fix(ql_state,qr_state,cl,cr, &
            wave_local, speeds)
        endif

        area = auxl(locarea,i)
        do mw = 1,mwaves
            call rotate2_tr(rot,wave_local(mw,2))
            speeds(1,mw) = area*speeds(1,mw)
            speeds(2,mw) = area*speeds(2,mw)
        enddo

        do m = 1,meqn
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw = 1,mwaves
                wave(m,mw,i) = wave_local(mw,m)
                s(mw,i) = speeds(1,mw) + speeds(2,mw)
                amdq(m,i) = amdq(m,i) + speeds(1,mw)*wave(m,mw,i)
                apdq(m,i) = apdq(m,i) + speeds(2,mw)*wave(m,mw,i)
            enddo
        enddo
    enddo


    return
    end subroutine rpn2

    subroutine apply_entropy_fix(ql,qr,cl, cr,wave_local,speeds)
    implicit none

    double precision :: ql(4), qr(4),wave_local(3,4)
    double precision :: speeds(2,3), s1, s2, s3, s4
    double precision :: sl, sml, smr, sr, qml(4),qmr(4)
    double precision :: ul, cl, ur, cr, pml, pmr, cml,cmr
    double precision :: sfract
    double precision :: gamma, gamma1
    logical :: trans1, trans3

    integer :: m

    common /param/ gamma

    gamma1 = gamma - 1.d0

    s1 = speeds(1,1) + speeds(2,1)
    s2 = speeds(1,2) + speeds(2,2)
    s3 = speeds(1,3) + speeds(2,3)

    do m = 1,4
        qml(m) = ql(m) + wave_local(1,m)
    enddo
    sl = ql(2)/ql(1) - cl
    pml = gamma1*(qml(4) - 0.5d0*(qml(2)**2/qml(1) + &
    qml(3)**2/qml(1)))
    sml = qml(2)/qml(1) - sqrt(gamma*pml/qml(1))

!     # Check the 1-wave
    trans1 = .false.
    if (sl < 0 .AND. sml > 0) then
    !        # apply transonic entropy fix
        trans1 = .true.
        sfract = (sml - s1)/(sml - sl)
        speeds(1,1) = sfract*sl
        speeds(2,1) = (1-sfract)*sml
    endif

!     # If the 1-wave is transonic,then we are done...
!     # Otherwise, we have to check the 3-wave
    if ( .NOT. trans1) then
        do m = 1,4
            qmr(m) = qr(m) - wave_local(3,m)
        enddo
        sr = qr(2)/qr(1) + cr
        pmr = gamma1*(qmr(4) - 0.5d0*(qmr(2)**2/qmr(1) + &
        qmr(3)**2/qmr(1)))
        smr = qmr(2)/qmr(1) + sqrt(gamma*pmr/qmr(1))
        trans3 = .false.
        if (smr < 0 .AND. sr > 0) then
        !           # apply transonic entropy fix
            trans3 = .true.
            sfract = (sr - s3)/(sr - smr)
            speeds(1,3) = sfract*smr
            speeds(2,3) = (1-sfract)*sr
        endif
    endif

    if (trans1) then
        write(6,*) '1-wave is transonic'
    elseif (trans3) then
        write(6,*) '3-wave is transonic'
    endif

    end subroutine apply_entropy_fix
