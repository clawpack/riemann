#include "cudaclaw/arrayIndex.H"

#define MAXITER 1

! This module should only be used when CUDA is enable
! pressure_forcing from storm_module is not supported.
module shallow_topo_device_module

#ifdef CUDA
    real(CLAW_REAL), parameter :: g = 9.81
    real(CLAW_REAL), parameter :: drytol = 0.001
    real(CLAW_REAL), parameter :: rho = 1025.0
    attributes(constant) :: g,drytol,rho ! in device constant memory
#ifdef USE_CAPA
    real(CLAW_REAL), parameter :: earth_radius = 6367.5E3
    real(CLAW_REAL), parameter :: deg2rad = 4.d0*datan(1.d0) / 180.d0
    attributes(constant) :: earth_radius,deg2rad
#endif
#else
    ! should not be used
    use geoclaw_module, only: g => grav, drytol => dry_tolerance, rho
    use geoclaw_module, only: earth_radius, deg2rad
    use storm_module, only: pressure_forcing, pressure_index
    use amr_module, only: mcapa
#endif


    contains

#ifdef CUDA
    attributes(device) &
#endif
    subroutine riemann_shallow_topo_x(q_l, q_r, aux_l, aux_r,    &
            fwave, s) bind (C,name='riemann_shallow_topo_x')

        implicit none

        real(CLAW_REAL), intent(inout) :: q_l(NEQNS), q_r(NEQNS)
        real(CLAW_REAL), intent(in) :: aux_r(NCOEFFS), aux_l(NCOEFFS)

        ! Output arguments
        real(CLAW_REAL), intent(inout) :: fwave(*)
        real(CLAW_REAL), intent(inout) :: s(*)

        !local only
        integer m,mw
        real(CLAW_REAL) wall(3)
        real(CLAW_REAL) fw(3,3)
        real(CLAW_REAL) sw(3)

        real(CLAW_REAL) hR,hL,uR,uL,vR,vL,phiR,phiL
        real(CLAW_REAL) bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
        real(CLAW_REAL) s1m,s2m
        real(CLAW_REAL) hstar,hstartest,hstarHLL
#ifdef USE_CAPA
        real(CLAW_REAL) dxdc
#endif

        logical rare1,rare2

        !-----------------------Initializing-----------------------------------

        !Initialize Riemann problem for grid interface
        do mw=1,NWAVES
            s(GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, mw, NWAVES, blockDim%y, blockDim%x)) = 0.d0 
            fwave(GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 1, NWAVES, NEQNS, blockDim%y, blockDim%x)) = 0.d0 
            fwave(GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 2, NWAVES, NEQNS, blockDim%y, blockDim%x)) = 0.d0 
            fwave(GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 3, NWAVES, NEQNS, blockDim%y, blockDim%x)) = 0.d0 
        enddo


        !skip problem if in a completely dry area
        if (q_l(1) > drytol .or. q_r(1) > drytol) then

            !Riemann problem variables
            hL = q_l(1) 
            hR = q_r(1) 
            bL = aux_l(1)
            bR = aux_r(1)

            !check for wet/dry boundary
            if (hR.gt.drytol) then
                uR=q_r(2)/hR
                vR=q_r(3)/hR
                phiR = 0.5d0*g*hR**2 + q_r(2)**2/hR
            else
                hR = 0.d0
                uR = 0.d0
                vR = 0.d0
                phiR = 0.d0
            endif

            if (hL.gt.drytol) then
                uL=q_l(2)/hL
                vL=q_l(3)/hL
                phiL = 0.5d0*g*hL**2 + q_l(2)**2/hL
            else
                hL=0.d0
                uL=0.d0
                vL=0.d0
                phiL = 0.d0
            endif

            wall(1) = 1.d0
            wall(2) = 1.d0
            wall(3) = 1.d0
            if (hR.le.drytol) then
                call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m, &
                         rare1,rare2)
                hstartest=max(hL,hstar)
                if (hstartest+bL.lt.bR) then 
                ! hL+bL < bR && hstar+bL < bR, so water can't overtop right cell (move into right cell)
                ! so right state should become ghost values that mirror left for wall problem
                    !                bR=hstartest+bL
                    wall(2)=0.d0
                    wall(3)=0.d0
                    hR=hL
                    bR=bL
                    phiR=phiL
                    uR=-uL
                    vR=vL
                    ! Here we already have huR=-huL
                elseif (hL+bL.lt.bR) then 
                ! hL+bL < bR && hstar+bL > bR, so we set bR to the water level in the left cell
                ! so that water can possibly overtop the right cell (move into the right cell)
                    bR=hL+bL
                ! otherwise
                ! hL+bL > bR || hstar+bL > bR, don't have to do anything about bR 
                endif
            elseif (hL.le.drytol) then ! right surface is lower than left topo
                call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m, &
                         rare1,rare2)
                hstartest=max(hR,hstar)
                if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
                    !               bL=hstartest+bR
                    wall(1)=0.d0
                    wall(2)=0.d0
                    hL=hR
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

            call riemann_aug_JCP(hL,hR, &
                     bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2, &
                                                 sw,fw)

            !        !eliminate ghost fluxes for wall
            do mw=1,NWAVES
                sw(mw)=sw(mw)*wall(mw)

                fw(1,mw)=fw(1,mw)*wall(mw) 
                fw(2,mw)=fw(2,mw)*wall(mw)
                fw(3,mw)=fw(3,mw)*wall(mw)
            enddo

#ifdef USE_CAPA
            dxdc=(earth_radius*deg2rad)
#endif
            do mw=1,NWAVES
                s(GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, mw, NWAVES, blockDim%y, blockDim%x)) = &
#ifdef USE_CAPA
                    sw(mw)*dxdc
#else
                    sw(mw)
#endif
                fwave( &
                    GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
                    ) = &
#ifdef USE_CAPA
                    fw(1,mw)*dxdc
#else
                    fw(1,mw)
#endif
                fwave( &
                    GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
                    ) = &
#ifdef USE_CAPA
                    fw(2,mw)*dxdc
#else
                    fw(2,mw)
#endif
                fwave( &
                    GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
                    ) = &
#ifdef USE_CAPA
                    fw(3,mw)*dxdc
#else
                    fw(3,mw)
#endif
            enddo
        endif

        !===============================================================================
        return
    end subroutine riemann_shallow_topo_x

#ifdef CUDA
    attributes(device) &
#endif
    subroutine riemann_shallow_topo_y(q_l, q_r, aux_l, aux_r,    &
            fwave, s) bind (C,name='riemann_shallow_topo_y')

        implicit none

        real(CLAW_REAL), intent(inout) :: q_l(NEQNS), q_r(NEQNS)
        real(CLAW_REAL), intent(in) :: aux_r(NCOEFFS), aux_l(NCOEFFS)

        ! Output arguments
        real(CLAW_REAL), intent(inout) :: fwave(*)
        real(CLAW_REAL), intent(inout) :: s(*)

        !local only
        integer m,mw
        real(CLAW_REAL) wall(3)
        real(CLAW_REAL) fw(3,3)
        real(CLAW_REAL) sw(3)

        real(CLAW_REAL) hR,hL,uR,uL,vR,vL,phiR,phiL
        real(CLAW_REAL) bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
        real(CLAW_REAL) s1m,s2m
        real(CLAW_REAL) hstar,hstartest,hstarHLL
#ifdef USE_CAPA
        real(CLAW_REAL) dydc
#endif

        logical rare1,rare2

        !-----------------------Initializing-----------------------------------

        !Initialize Riemann problem for grid interface
        do mw=1,NWAVES
            s(GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, mw, NWAVES, blockDim%y, blockDim%x)) = 0.d0 
            fwave(GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 1, NWAVES, NEQNS, blockDim%y, blockDim%x)) = 0.d0 
            fwave(GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 2, NWAVES, NEQNS, blockDim%y, blockDim%x)) = 0.d0 
            fwave(GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 3, NWAVES, NEQNS, blockDim%y, blockDim%x)) = 0.d0 
        enddo


        !skip problem if in a completely dry area
        if (q_l(1) > drytol .or. q_r(1) > drytol) then

            !Riemann problem variables
            hL = q_l(1) 
            hR = q_r(1) 
            bL = aux_l(1)
            bR = aux_r(1)

            !check for wet/dry boundary
            if (hR.gt.drytol) then
                uR=q_r(3)/hR
                vR=q_r(2)/hR
                phiR = 0.5d0*g*hR**2 + q_r(3)**2/hR
            else
                hR = 0.d0
                uR = 0.d0
                vR = 0.d0
                phiR = 0.d0
            endif

            if (hL.gt.drytol) then
                uL=q_l(3)/hL
                vL=q_l(2)/hL
                phiL = 0.5d0*g*hL**2 + q_l(3)**2/hL
            else
                hL=0.d0
                uL=0.d0
                vL=0.d0
                phiL = 0.d0
            endif

            wall(1) = 1.d0
            wall(2) = 1.d0
            wall(3) = 1.d0
            if (hR.le.drytol) then
                call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m, &
                         rare1,rare2)
                hstartest=max(hL,hstar)
                if (hstartest+bL.lt.bR) then 
                ! hL+bL < bR && hstar+bL < bR, so water can't overtop right cell (move into right cell)
                ! so right state should become ghost values that mirror left for wall problem
                    !                bR=hstartest+bL
                    wall(2)=0.d0
                    wall(3)=0.d0
                    hR=hL
                    bR=bL
                    phiR=phiL
                    uR=-uL
                    vR=vL
                    ! Here we already have huR=-huL
                elseif (hL+bL.lt.bR) then 
                ! hL+bL < bR && hstar+bL > bR, so we set bR to the water level in the left cell
                ! so that water can possibly overtop the right cell (move into the right cell)
                    bR=hL+bL
                ! otherwise
                ! hL+bL > bR || hstar+bL > bR, don't have to do anything about bR 
                endif
            elseif (hL.le.drytol) then ! right surface is lower than left topo
                call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m, &
                         rare1,rare2)
                hstartest=max(hR,hstar)
                if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
                    !               bL=hstartest+bR
                    wall(1)=0.d0
                    wall(2)=0.d0
                    hL=hR
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

            call riemann_aug_JCP(hL,hR, &
                     bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2, &
                                                 sw,fw)

            ! eliminate ghost fluxes for wall
            ! TODO: merge this loop with the one below
            do mw=1,NWAVES
                sw(mw)=sw(mw)*wall(mw)

                fw(1,mw)=fw(1,mw)*wall(mw) 
                fw(2,mw)=fw(2,mw)*wall(mw)
                fw(3,mw)=fw(3,mw)*wall(mw)
            enddo

#ifdef USE_CAPA
            dydc = earth_radius*cos(aux_r(3))*deg2rad
#endif
            do mw=1,NWAVES
                s(GET_INDEX_SHARED_SPEED_1INDEX(threadIdx%y, threadIdx%x, mw, NWAVES, blockDim%y, blockDim%x)) = &
#ifdef USE_CAPA
                    sw(mw)*dydc
#else
                    sw(mw)
#endif
                fwave( &
                    GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 1, NWAVES, NEQNS, blockDim%y, blockDim%x) &
                    ) = &
#ifdef USE_CAPA
                    fw(1,mw)*dydc
#else
                    fw(1,mw)
#endif
                fwave( &
                    GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 3, NWAVES, NEQNS, blockDim%y, blockDim%x) &
                    ) = &
#ifdef USE_CAPA
                    fw(2,mw)*dydc
#else
                    fw(2,mw)
#endif
                fwave( &
                    GET_INDEX_SHARED_WAVE_1INDEX(threadIdx%y, threadIdx%x, mw, 2, NWAVES, NEQNS, blockDim%y, blockDim%x) &
                    ) = &
#ifdef USE_CAPA
                    fw(3,mw)*dydc
#else
                    fw(3,mw)
#endif
            enddo
        endif

        !===============================================================================
        return
    end subroutine riemann_shallow_topo_y

#ifdef CUDA
attributes(device) &
#endif
    subroutine riemann_aug_JCP(hL,hR, &
              bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2, &
              sw,fw)

      ! solve shallow water equations given single left and right states
      ! This solver is described in J. Comput. Phys. (6): 3089-3113, March 2008
      ! Augmented Riemann Solvers for the Shallow Equations with Steady States and Inundation

      implicit none

      !input
      real(CLAW_REAL) fw(NEQNS,NWAVES)
      real(CLAW_REAL) sw(NWAVES)
      real(CLAW_REAL) hL,hR,bL,bR,uL,uR,phiL,phiR,sE1,sE2
      real(CLAW_REAL) vL,vR


      !local
      integer m,mw,k,iter
      real(CLAW_REAL) A(3,3)
      real(CLAW_REAL) r(3,3)
      real(CLAW_REAL) lambda(3)
      real(CLAW_REAL) del(3)
      real(CLAW_REAL) beta(3)

      real(CLAW_REAL) delh,delhu,delphi,delb,delnorm
      real(CLAW_REAL) rare1st,rare2st,sdelta,raremin,raremax
      real(CLAW_REAL) criticaltol,convergencetol,raretol
      real(CLAW_REAL) criticaltol_2, hustar_interface
      real(CLAW_REAL) s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
      real(CLAW_REAL) huRstar,huLstar,uRstar,uLstar,hstarHLL
      real(CLAW_REAL) deldelh,deldelphi
      real(CLAW_REAL) s1m,s2m,hm
      real(CLAW_REAL) det1,det2,det3,determinant

      logical rare1,rare2,rarecorrector,rarecorrectortest,sonic

      !determine del vectors
      delh = hR-hL
      delhu = hR*uR-hL*uL
      delphi = phiR-phiL
      delb = bR-bL
      delnorm = delh**2 + delphi**2

      call riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2)


      lambda(1)= min(sE1,s2m) !Modified Einfeldt speed
      lambda(3)= max(sE2,s1m) !Modified Eindfeldt speed
      sE1=lambda(1)
      sE2=lambda(3)
      lambda(2) = 0.d0  ! ### Fix to avoid uninitialized value in loop on mw -- Correct?? ###

      
      hstarHLL = max((hL*uL-hR*uR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve

!     !determine the middle entropy corrector wave------------------------
      rarecorrectortest=.false.
      rarecorrector=.false.
      if (rarecorrectortest) then
         sdelta=lambda(3)-lambda(1)
         raremin = 0.5d0
         raremax = 0.9d0
         if (rare1.and.sE1*s1m.lt.0.d0) raremin=0.2d0
         if (rare2.and.sE2*s2m.lt.0.d0) raremin=0.2d0
         if (rare1.or.rare2) then
            !see which rarefaction is larger
            rare1st=3.d0*(sqrt(g*hL)-sqrt(g*hm))
            rare2st=3.d0*(sqrt(g*hR)-sqrt(g*hm))
            if (max(rare1st,rare2st).gt.raremin*sdelta.and.  &
                max(rare1st,rare2st).lt.raremax*sdelta) then
                  rarecorrector=.true.
               if (rare1st.gt.rare2st) then
                  lambda(2)=s1m
               elseif (rare2st.gt.rare1st) then
                  lambda(2)=s2m
               else
                  lambda(2)=0.5d0*(s1m+s2m)
               endif
            endif
         endif
         if (hstarHLL.lt.min(hL,hR)/5.d0) rarecorrector=.false.
      endif

!     ## Is this correct 2-wave when rarecorrector == .true. ??
      do mw=1,NWAVES
         r(1,mw)=1.d0
         r(2,mw)=lambda(mw)
         r(3,mw)=(lambda(mw))**2
      enddo
      if (.not.rarecorrector) then
         lambda(2) = 0.5d0*(lambda(1)+lambda(3))
!         lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
         r(1,2)=0.d0
         r(2,2)=0.d0
         r(3,2)=1.d0
      endif
!     !---------------------------------------------------


!     !determine the steady state wave -------------------
      !criticaltol = 1.d-6
      ! MODIFIED:
      criticaltol = max(drytol*g, 1d-6)
      criticaltol_2 = sqrt(criticaltol)
      deldelh = -delb
      deldelphi = -0.5d0 * (hR + hL) * (g * delb)

!     !determine a few quanitites needed for steady state wave if iterated
      hLstar=hL
      hRstar=hR
      uLstar=uL
      uRstar=uR
      huLstar=uLstar*hLstar
      huRstar=uRstar*hRstar

      !iterate to better determine the steady state wave
      convergencetol=1.d-6
      do iter=1,MAXITER
         !determine steady state wave (this will be subtracted from the delta vectors)
         if (min(hLstar,hRstar).lt.drytol.and.rarecorrector) then
              rarecorrector=.false.
              hLstar=hL
              hRstar=hR
              uLstar=uL
              uRstar=uR
              huLstar=uLstar*hLstar
              huRstar=uRstar*hRstar
              lambda(2) = 0.5d0*(lambda(1)+lambda(3))
              !           lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
              r(1,2)=0.d0
              r(2,2)=0.d0
              r(3,2)=1.d0
          endif

          hbar =  max(0.5d0*(hLstar+hRstar),0.d0)
          s1s2bar = 0.25d0*(uLstar+uRstar)**2 - g*hbar
          s1s2tilde= max(0.d0,uLstar*uRstar) - g*hbar

          !        !find if sonic problem
          ! MODIFIED from 5.3.1 version
          sonic = .false.
          if (abs(s1s2bar) <= criticaltol) then
              sonic = .true.
          else if (s1s2bar*s1s2tilde <= criticaltol**2) then
              sonic = .true.
          else if (s1s2bar*sE1*sE2 <= criticaltol**2) then
              sonic = .true.
          else if (min(abs(sE1),abs(sE2)) < criticaltol_2) then
              sonic = .true.
          else if (sE1 <  criticaltol_2 .and. s1m > -criticaltol_2) then
              sonic = .true.
          else if (sE2 > -criticaltol_2 .and. s2m <  criticaltol_2) then
              sonic = .true.
          else if ((uL+dsqrt(g*hL))*(uR+dsqrt(g*hR)) < 0.d0) then
              sonic = .true.
          else if ((uL- dsqrt(g*hL))*(uR- dsqrt(g*hR)) < 0.d0) then
              sonic = .true.
          end if

          !        !find jump in h, deldelh
          if (sonic) then
              deldelh =  -delb
          else
              deldelh = delb*g*hbar/s1s2bar
          endif
          !        !find bounds in case of critical state resonance, or negative states
          if (sE1.lt.-criticaltol.and.sE2.gt.criticaltol) then
              deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
              deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
          elseif (sE1.ge.criticaltol) then
              deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
              deldelh = max(deldelh,-hL)
          elseif (sE2.le.-criticaltol) then
              deldelh = min(deldelh,hR)
              deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
          endif

          !        !find jump in phi, deldelphi
          if (sonic) then
              deldelphi = -g*hbar*delb
          else
              deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar
          endif
          !        !find bounds in case of critical state resonance, or negative states
          deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
          deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))
          deldelphi = deldelphi

          del(1)=delh-deldelh
          del(2)=delhu
          del(3)=delphi-deldelphi

          !        !Determine determinant of eigenvector matrix========
          det1=r(1,1)*(r(2,2)*r(3,3)-r(2,3)*r(3,2))
          det2=r(1,2)*(r(2,1)*r(3,3)-r(2,3)*r(3,1))
          det3=r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))
          determinant=det1-det2+det3

          !        !solve for beta(k) using Cramers Rule=================
          do k=1,3
              do mw=1,3
                  A(1,mw)=r(1,mw)
                  A(2,mw)=r(2,mw)
                  A(3,mw)=r(3,mw)
              enddo
              A(1,k)=del(1)
              A(2,k)=del(2)
              A(3,k)=del(3)
              det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
              det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
              det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
              beta(k)=(det1-det2+det3)/determinant
          enddo

          !exit if things aren't changing
          if (abs(del(1)**2+del(3)**2-delnorm).lt.convergencetol) exit 
          delnorm = del(1)**2+del(3)**2
          !find new states qLstar and qRstar on either side of interface
          hLstar=hL
          hRstar=hR
          uLstar=uL
          uRstar=uR
          huLstar=uLstar*hLstar
          huRstar=uRstar*hRstar
          do mw=1,NWAVES
              if (lambda(mw).lt.0.d0) then
                  hLstar= hLstar + beta(mw)*r(1,mw)
                  huLstar= huLstar + beta(mw)*r(2,mw)
              endif
          enddo
          do mw=NWAVES,1,-1
              if (lambda(mw).gt.0.d0) then
                  hRstar= hRstar - beta(mw)*r(1,mw)
                  huRstar= huRstar - beta(mw)*r(2,mw)
              endif
          enddo

          if (hLstar.gt.drytol) then
              uLstar=huLstar/hLstar
          else
              hLstar=max(hLstar,0.d0)
              uLstar=0.d0
          endif
          if (hRstar.gt.drytol) then
              uRstar=huRstar/hRstar
          else
              hRstar=max(hRstar,0.d0)
              uRstar=0.d0
          endif

      enddo ! end iteration on Riemann problem

      do mw=1,NWAVES
         sw(mw)=lambda(mw)
         fw(1,mw)=beta(mw)*r(2,mw)
         fw(2,mw)=beta(mw)*r(3,mw)
         fw(3,mw)=beta(mw)*r(2,mw)
      enddo
      !find transverse components (ie huv jumps).
      ! MODIFIED from 5.3.1 version
      fw(3,1)=fw(3,1)*vL
      fw(3,3)=fw(3,3)*vR
      fw(3,2)= 0.d0
 
      hustar_interface = hL*uL + fw(1,1)   ! = huR - fw(1,3)
      if (hustar_interface <= 0.0d0) then
          fw(3,1) = fw(3,1) + (hR*uR*vR - hL*uL*vL - fw(3,1)- fw(3,3))
        else
          fw(3,3) = fw(3,3) + (hR*uR*vR - hL*uL*vL - fw(3,1)- fw(3,3))
        end if


      return

      end !subroutine riemann_aug_JCP-------------------------------------------------


!-----------------------------------------------------------------------


!=============================================================================
#ifdef CUDA
attributes(device) &
#endif
      subroutine riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2)

      ! TODO: do I need double precision in this if maxiter = 1 anyway?
      !determine the Riemann structure (wave-type in each family)

      implicit none

      !input
      real(CLAW_REAL) hL,hR,uL,uR

      !output 
      ! 1 and 2 represent wave speeds relavant to 1st field and 2nd field
      real(CLAW_REAL) s1m,s2m,hm
      logical rare1,rare2

      !local
      real(CLAW_REAL) u1m,u2m,um,delu
      real(CLAW_REAL) h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
      ! temporary
      real(CLAW_REAL) :: sqrtgh1,sqrtgh2
      integer iter



!     !Test for Riemann structure

      h_min=min(hR,hL)
      h_max=max(hR,hL)
      delu=uR-uL

      ! Have dry state on either side
      ! Only one rarefaction wave
      ! another shock wave has 0 jump and moves at the same speed
      ! as one edge of the rarefaction wave
      if (h_min.le.drytol) then
         hm=0.d0
         um=0.d0
         ! Eq. (54) in the JCP paper. 
         ! Either hR or hL is almost zero
         ! So the expression below corresponds to either Eq. (54a)
         ! or Eq. (54b)
         s1m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
         s2m=s1m
         if (hL.le.0.d0) then
            rare2=.true.
            rare1=.false.
         else
            rare1=.true.
            rare2=.false.
         endif

      else
         F_min= delu+2.d0*(sqrt(g*h_min)-sqrt(g*h_max))
         F_max= delu + &
             (h_max-h_min)*(sqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

         if (F_min.gt.0.d0) then !2-rarefactions

            ! (13.56) in the FVMHP book
            hm=(1.d0/(16.d0*g))* &
                max(0.d0,-delu+2.d0*(sqrt(g*hL)+sqrt(g*hR)))**2
            um=sign(1.d0,hm)*(uL+2.d0*(sqrt(g*hL)-sqrt(g*hm)))

            ! wave speed of the right edge of 1-rarefaction
            s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
            ! wave speed of the left edge of 2-rarefaction
            s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)

            rare1=.true.
            rare2=.true.

         elseif (F_max.le.0.d0) then !2 shocks
            ! Below it solves for the intersection of two Hugoniot loci
            ! to get accurate Riemann solution

!           !root finding using a Newton iteration on sqrt(h)===
            h0=h_max
            do iter=1,MAXITER
               gL=sqrt(.5d0*g*(1/h0 + 1/hL))
               gR=sqrt(.5d0*g*(1/h0 + 1/hR))
               F0=delu+(h0-hL)*gL + (h0-hR)*gR
               dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+ &
                   gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
               slope=2.d0*sqrt(h0)*dfdh
               h0=(sqrt(h0)-F0/slope)**2
            enddo
               hm=h0
               ! u1m and u2m are (13.19) and (13.20) in the FVMHP book
               ! in this case we should have u1m ~= u2m?
               u1m=uL-(hm-hL)*sqrt((.5d0*g)*(1/hm + 1/hL))
               u2m=uR+(hm-hR)*sqrt((.5d0*g)*(1/hm + 1/hR))
               um=.5d0*(u1m+u2m)

               ! QUESTION: I can't derive this.
               ! This seems to be the speeds of the two shocks
               s1m=u1m-sqrt(g*hm)
               s2m=u2m+sqrt(g*hm)
               rare1=.false.
               rare2=.false.

         else !one shock one rarefaction
            h0=h_min

            do iter=1,MAXITER
               F0=delu + 2.d0*(sqrt(g*h0)-sqrt(g*h_max)) &
                   + (h0-h_min)*sqrt(.5d0*g*(1/h0+1/h_min))
               slope=(F_max-F0)/(h_max-h_min)
               h0=h0-F0/slope
            enddo

            hm=h0
            sqrtgh2 = sqrt(g*hm)
            if (hL.gt.hR) then
                sqrtgh1 = sqrt(g*hL)
               ! ! Eq. (13.55) in the FVMHP book
               ! um=uL+2.d0*sqrt(g*hL)-2.d0*sqrt(g*hm)
               ! s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
               ! ! QUESTION: what's s2m here?
               ! s2m=uL+2.d0*sqrt(g*hL)-sqrt(g*hm)

                um=uL+2.d0*sqrtgh1-2.d0*sqrtgh2
                s1m=uL+2.d0*sqrtgh1-3.d0*sqrtgh2
                s2m=uL+2.d0*sqrtgh1-sqrtgh2

                rare1=.true.
                rare2=.false.
            else
               sqrtgh1 = sqrt(g*hR)
               ! s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)
               ! ! QUESTION: what's s1m here?
               ! s1m=uR-2.d0*sqrt(g*hR)+sqrt(g*hm)
               ! ! Eq. (13.55) in the FVMHP book
               ! um=uR-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hm)

               s1m=uR-2.d0*sqrtgh1+sqrtgh2
               s2m=uR-2.d0*sqrtgh1+3.d0*sqrtgh2
               um=uR-2.d0*sqrtgh1+2.d0*sqrtgh2
               rare2=.true.
               rare1=.false.
            endif
         endif
      endif

      return
      end ! subroutine riemanntype----------------------------------------------------------------
end module shallow_topo_device_module
