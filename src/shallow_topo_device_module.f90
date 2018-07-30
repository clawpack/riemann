! #include "cudaclaw/arrayIndex.H"

#define GRAVITY 1.d0

module shallow_topo_device_module
      use geoclaw_module, only: g => grav, drytol => dry_tolerance, rho
      use geoclaw_module, only: earth_radius, deg2rad
      use amr_module, only: mcapa

      use storm_module, only: pressure_forcing, pressure_index

    contains

    subroutine riemann_shallow_topo(ixy, q_l, q_r, aux_l, aux_r,    &
            fwave, s)

        implicit none

        integer, intent(in) :: ixy
        real(CLAW_REAL), intent(in) :: q_l(NEQNS), q_r(NEQNS)
        real(CLAW_REAL), intent(in) :: aux_r(NCOEFFS), aux_l(NCOEFFS)

        ! Output arguments
        real(CLAW_REAL), intent(inout) :: fwave(NEQNS,NWAVES)
        real(CLAW_REAL), intent(inout) :: s(NWAVES)

        !local only
        integer m,mw,maxiter,mu,nv
        double precision wall(3)
        double precision fw(3,3)
        double precision sw(3)

        double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL,pL,pR
        double precision bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
        double precision s1m,s2m
        double precision hstar,hstartest,hstarHLL,sLtest,sRtest
        double precision tw,dxdc

        logical rare1,rare2

        ! In case there is no pressure forcing
        pL = 0.d0
        pR = 0.d0

        !-----------------------Initializing-----------------------------------

        !Initialize Riemann problem for grid interface
        do mw=1,NWAVES
            s(mw)=0.d0
            fwave(1,mw)=0.d0
            fwave(2,mw)=0.d0
            fwave(3,mw)=0.d0
        enddo

        !        !set normal direction
        if (ixy.eq.1) then
            mu=2
            nv=3
        else
            mu=3
            nv=2
        endif

        !skip problem if in a completely dry area
        if (q_l(1) > drytol .or. q_r(1) > drytol) then

            !Riemann problem variables
            hL = q_l(1) 
            hR = q_r(1) 
            huL = q_l(mu) 
            huR = q_r(mu) 
            bL = aux_l(1)
            bR = aux_r(1)
            if (pressure_forcing) then
                pL = aux_l(pressure_index)
                pR = aux_r(pressure_index)
            end if

            hvL=q_l(nv) 
            hvR=q_r(nv)

            !check for wet/dry boundary
            if (hR.gt.drytol) then
                uR=huR/hR
                vR=hvR/hR
                phiR = 0.5d0*g*hR**2 + huR**2/hR
            else
                hR = 0.d0
                huR = 0.d0
                hvR = 0.d0
                uR = 0.d0
                vR = 0.d0
                phiR = 0.d0
            endif

            if (hL.gt.drytol) then
                uL=huL/hL
                vL=hvL/hL
                phiL = 0.5d0*g*hL**2 + huL**2/hL
            else
                hL=0.d0
                huL=0.d0
                hvL=0.d0
                uL=0.d0
                vL=0.d0
                phiL = 0.d0
            endif

            wall(1) = 1.d0
            wall(2) = 1.d0
            wall(3) = 1.d0
            if (hR.le.drytol) then
                call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m, &
                         rare1,rare2,1,drytol,g)
                hstartest=max(hL,hstar)
                ! hL+bL < bR && hstar+bL < bR, so water can't overtop right cell (move into right cell)
                if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
                    !                bR=hstartest+bL
                    wall(2)=0.d0
                    wall(3)=0.d0
                    hR=hL
                    huR=-huL
                    bR=bL
                    phiR=phiL
                    uR=-uL
                    vR=vL
                    ! hL+bL < bR && hstar+bL > bR, so we set bR to the water level in the left cell
                    ! so that water can possibly overtop the right cell (move into the right cell)
                elseif (hL+bL.lt.bR) then 
                    bR=hL+bL
                    ! hL+bL > bR || hstar+bL > bR, don't have to do anything about bR 
                endif
            elseif (hL.le.drytol) then ! right surface is lower than left topo
                call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m, &
                         rare1,rare2,1,drytol,g)
                hstartest=max(hR,hstar)
                if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
                    !               bL=hstartest+bR
                    wall(1)=0.d0
                    wall(2)=0.d0
                    hL=hR
                    huL=-huR
                    bL=bR
                    phiL=phiR
                    uL=-uR
                    vL=vR
                elseif (hR+bR.lt.bL) then
                    bL=hR+bR
                endif
            endif

            !determine wave speeds
            sL=uL-sqrt(g*hL) ! 1 wave speed of left state
            sR=uR+sqrt(g*hR) ! 2 wave speed of right state

            uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
            chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
            sRoe1=uhat-chat ! Roe wave speed 1 wave
            sRoe2=uhat+chat ! Roe wave speed 2 wave

            sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
            sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

            !--------------------end initializing...finally----------
            !solve Riemann problem.

            maxiter = 1

            call riemann_aug_JCP(maxiter,3,3,hL,hR,huL, &
                     huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2, &
                                                 drytol,g,rho,sw,fw)

            !          call riemann_ssqfwave(maxiter,NEQNS,NWAVES,hL,hR,huL,huR,
            !      &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,
            !      &     rho,sw,fw)

            !          call riemann_fwave(NEQNS,NWAVES,hL,hR,huL,huR,hvL,hvR,
            !      &      bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,rho,sw,
            !      &      fw)

            !        !eliminate ghost fluxes for wall
            do mw=1,NWAVES
                sw(mw)=sw(mw)*wall(mw)

                fw(1,mw)=fw(1,mw)*wall(mw) 
                fw(2,mw)=fw(2,mw)*wall(mw)
                fw(3,mw)=fw(3,mw)*wall(mw)
            enddo

            do mw=1,NWAVES
                s(mw)=sw(mw)
                fwave(1,mw)=fw(1,mw)
                fwave(mu,mw)=fw(2,mw)
                fwave(nv,mw)=fw(3,mw)
            enddo
        endif


! ! mcapa > 0 is not supported for GPU version
! #ifndef CUDA
!         !==========Capacity for mapping from latitude longitude to physical space====
!         if (mcapa.gt.0) then
!             if (ixy.eq.1) then
!                 dxdc=(earth_radius*deg2rad)
!             else
!                 dxdc=earth_radius*cos(aux_r(3))*deg2rad
!             endif
! 
!             do mw=1,NWAVES
!                 s(mw)=dxdc*s(mw)
!                 fwave(1,mw)=dxdc*fwave(1,mw)
!                 fwave(2,mw)=dxdc*fwave(2,mw)
!                 fwave(3,mw)=dxdc*fwave(3,mw)
!             enddo
!         endif
! #endif

        !===============================================================================
        return
    end subroutine 
end module shallow_topo_device_module
