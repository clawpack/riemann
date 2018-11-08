! =========================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =========================================================

!     # solve Riemann problems for the 1D Euler equations using Roe's
!     # approximate Riemann solver.

!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell
!     # On output, wave contains the waves,
!     #      the speeds,
!     #            amdq the  left-going flux difference  A^- \Delta q
!     #            apdq the right-going flux difference  A^+ \Delta q

!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     # From the basic clawpack routine step1, rp is called with ql = qr = q.


    implicit none

    integer :: maxm, meqn, mwaves, mbc, mx, imp, ixy, maux
    double precision ::    ql(meqn,   1-mbc:maxm+mbc)
    double precision ::    qr(meqn,   1-mbc:maxm+mbc)
    double precision ::     s(mwaves, 1-mbc:maxm+mbc)
    double precision ::  wave(meqn,   mwaves, 1-mbc:maxm+mbc)
    double precision ::  asdq(meqn,   1-mbc:maxm+mbc)
    double precision ::  bmasdq(meqn, 1-mbc:maxm+mbc)
    double precision ::  bpasdq(meqn, 1-mbc:maxm+mbc)
    double precision ::  aux1(maux,      1-mbc:maxm+mbc)
    double precision ::  aux2(maux,      1-mbc:maxm+mbc)
    double precision ::  aux3(maux,      1-mbc:maxm+mbc)

!     # For Roe solver
    double precision :: rhol, ul, vl, el, cl, pl
    double precision :: rhor, ur, vr, er, cr, pr
    double precision :: u, v, enth, rho, p, e
    double precision :: rhsqrtl, rhsqrtr, rhsq2
    double precision :: uvl, uvr, uv(2)


!     # For mappings
    integer :: meqn2,mwaves2
    parameter(meqn2 = 4, mwaves2 = 3)
    double precision :: ql_state(meqn2), qr_state(meqn2)
    double precision :: q_state(meqn2)
    double precision :: wave_local(mwaves2,meqn2)
    double precision :: s_local(mwaves2), delta(meqn2)
    double precision :: speeds(2,mwaves2)
    double precision :: deltam(meqn2), deltap(meqn2)
    double precision :: area

!     # for mapping
    double precision :: rotm(4), rotp(4), uvm_rot(2),uvp_rot(2)
    integer :: locrot,mcapa,locarea

!     # Miscellaneous
    integer :: i, j, m, mw, i1, ixy1
    logical :: useroe

!     # Problem parameters
    double precision :: gamma, gamma1

    double precision :: dtcom, dxcom, dycom, tcom
    integer :: info

    logical :: in_rpt

    common /cparam/  gamma

    gamma1 = gamma - 1.d0

    useroe = .true.

    call get_aux_locations_t(ixy, mcapa, locrot,locarea)

    do i = 2-mbc,mx+mbc
        i1 = i + imp - 2

        do m = 1,meqn
            ql_state(m) = qr(m,i-1)
            qr_state(m) = ql(m,i)
        enddo

        if (useroe) then
            rhol = ql_state(1)
            rhor = qr_state(1)

            ul = ql_state(2)/rhol
            ur = qr_state(2)/rhor

            vl = ql_state(3)/rhol
            vr = qr_state(3)/rhor

            el = ql_state(4)
            er = qr_state(4)

            uvl = ul*ul + vl*vl
            uvr = ur*ur + vr*vr
            pl = gamma1*(el - 0.5d0*rhol*uvl)
            pr = gamma1*(er - 0.5d0*rhor*uvr)

        !           # Get Roe averaged values
            rhsqrtl = sqrt(rhol)
            rhsqrtr = sqrt(rhor)
            rhsq2 = rhsqrtl + rhsqrtr


            uv(1) = (ul*rhsqrtl + ur*rhsqrtr) / rhsq2
            uv(2) = (vl*rhsqrtl + vr*rhsqrtr) / rhsq2
            enth = ((el + pl)/rhsqrtl + (er + pr)/rhsqrtr)/rhsq2
        else
        !           # This takes the values needed for the Roe matrix from the
        !           # cell centers in either the  left (imp == 1) or the
        !           # right (imp == 2) cell

            do m = 1,meqn
                if (imp == 1) then
                !                 # Left (minus) cell
                    q_state(m) = ql_state(m)
                else
                !                 # Right (plus) cell
                    q_state(m) = qr_state(m)
                endif
            enddo

            rho = q_state(1)
            uv(1) = q_state(2)/rho
            uv(2) = q_state(3)/rho
            e = q_state(4)
            p = gamma1*(e - 0.5d0*rho*(uv(1)**2 + uv(2)**2))
            enth = (e + p)/rho
        endif

        do j = 1,2
            uvm_rot(j) = uv(j)
            uvp_rot(j) = uv(j)
            rotm(j) = aux2(locrot+j-1,i1)
            rotp(j) = aux3(locrot+j-1,i1)
        enddo
        call compute_tangent(rotm)
        call compute_tangent(rotp)

        call rotate2(rotm,uvm_rot)
        call rotate2(rotp,uvp_rot)

        do m = 1,meqn
            deltap(m) = asdq(m,i)
            deltam(m) = asdq(m,i)
        enddo
        call rotate2(rotm,deltam(2))
        call rotate2(rotp,deltap(2))

    !        # ------------------------------------------
    !        # Solve for minus side
    !        # ------------------------------------------
        ixy1 = 1
        call roe_solver(ixy1,uvm_rot,enth,deltam, &
        wave_local,s_local,info)

        if (info /= 0) then
            write(6,*) 'ixy = ', ixy
            write(6,*) 'imp = ', imp
            write(6,*) 'enth = ', enth
            write(6,*) 'Called from rpt2  (A-DQ)'
            write(6,*) ' '
            write(6,'(24A,24A)') '        Left State      ', &
            '       Right State      '
            write(6,'(24A,24A)') '------------------------', &
            '-----------------------'
            do m = 1,meqn
                write(6,'(2E24.16)') ql_state(m), qr_state(m)
            enddo
            write(6,*) ' '
            stop
        endif

        area = aux2(locarea,i1)
        do mw = 1,mwaves
            call rotate2_tr(rotm,wave_local(mw,2))
            speeds(1,mw) = area*min(s_local(mw),0.d0)
        enddo

        do m = 1,meqn
            bmasdq(m,i) = 0.d0
            do mw = 1,mwaves
                bmasdq(m,i) = bmasdq(m,i) + speeds(1,mw)*wave_local(mw,m)
            enddo
        enddo

    !        # ------------------------------------------
    !        # Solve for plus side
    !        # ------------------------------------------
        ixy1 = 1
        call roe_solver(ixy1,uvp_rot,enth,deltap, &
        wave_local,s_local,info)

        if (info /= 0) then
            write(6,*) 'ixy = ', ixy
            write(6,*) 'imp = ', imp
            write(6,*) 'enth = ', enth
            write(6,*) 'Called from rpt2  (A+DQ)'
            write(6,*) ' '
            write(6,'(24A,24A)') '        Left State      ', &
            '       Right State      '
            write(6,'(24A,24A)') '------------------------', &
            '-----------------------'
            do m = 1,meqn
                write(6,'(2E24.16)') ql_state(m), qr_state(m)
            enddo
            write(6,*) ' '
            stop
        endif

        area = aux3(locarea,i1)
        do mw = 1,mwaves
            call rotate2_tr(rotp,wave_local(mw,2))
            speeds(2,mw) = area*max(s_local(mw),0.d0)
        enddo

        do m = 1,meqn
            bpasdq(m,i) = 0.d0
            do mw = 1,mwaves
                bpasdq(m,i) = bpasdq(m,i) + speeds(2,mw)*wave_local(mw,m)
            enddo
        enddo

    enddo !! end of i loop

    return
    end subroutine rpt2
