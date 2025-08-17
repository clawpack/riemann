subroutine rp1(maxm, num_eqn, num_waves, num_aux, num_ghost, &
                num_cells, ql, qr, auxl, auxr, fwave, s, amdq, apdq)

! Normal Riemann solver for the 1D SHALLOW WATER equations
!     with topography:
!     #        h_t + (hu)_x = 0                           #
!     #        (hu)_t + (hu^2 + 0.5gh^2)_x = -ghb_x      #

! This solver is based on David George's solver written for GeoClaw.
! It has been modified to be compatible with f2py (and thus PyClaw).

! waves:     2
! equations: 2

! Conserved quantities:
!       1 depth
!       2 x_momentum

! Auxiliary fields:
!       1 bathymetry

! The gravitational constant grav should be in the common block cparam.

! See http://www.clawpack.org/riemann.html for a detailed explanation
! of the Riemann solver API.

      implicit none

      real(kind=8) :: grav, g, dry_tolerance, drytol
      common /cparam/ grav, dry_tolerance

      integer, intent(in) :: maxm,num_eqn,num_aux,num_waves,num_ghost,num_cells

      real(kind=8), intent(inout) ::  ql(num_eqn, 1-num_ghost:maxm+num_ghost)
      real(kind=8), intent(inout) ::  qr(num_eqn, 1-num_ghost:maxm+num_ghost)
      real(kind=8), intent(in) ::  auxl(num_aux,1-num_ghost:maxm+num_ghost)
      real(kind=8), intent(in) ::  auxr(num_aux,1-num_ghost:maxm+num_ghost)

      real(kind=8), intent(out) :: fwave(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
      real(kind=8), intent(out) ::  s(num_waves, 1-num_ghost:maxm+num_ghost)
      real(kind=8), intent(out) ::  apdq(num_eqn,1-num_ghost:maxm+num_ghost)
      real(kind=8), intent(out) ::  amdq(num_eqn,1-num_ghost:maxm+num_ghost)

      !local only
      integer m,i,mw,maxiter
      real(kind=8) wall(2)
      real(kind=8) fw(3,3)
      real(kind=8) sw(3)

      real(kind=8) hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
      real(kind=8) bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      real(kind=8) s1m,s2m
      real(kind=8) hstar,hstartest,hstarHLL,sLtest,sRtest

      g = grav
      drytol = dry_tolerance

      !loop through Riemann problems at each grid cell
      do i=2-num_ghost,num_cells+num_ghost

      !-----------------------Initializing------------------------------
         !inform of a bad riemann problem from the start
         if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
         endif

         !Initialize Riemann problem for grid interface
         do mw=1,num_waves
             s(mw,i)=0.d0
             fwave(1,mw,i)=0.d0
             fwave(2,mw,i)=0.d0
         enddo

         !zero (small) negative values if they exist
         if (qr(1,i-1).lt.0.d0) then
               qr(1,i-1)=0.d0
               qr(2,i-1)=0.d0
         endif

         if (ql(1,i).lt.0.d0) then
               ql(1,i)=0.d0
               ql(2,i)=0.d0
         endif

         !skip problem if in a completely dry area
         if (qr(1,i-1) <= drytol .and. ql(1,i) <= drytol) then
            go to 30
         endif

         !Riemann problem variables
         hL = qr(1,i-1) 
         hR = ql(1,i) 
         huL = qr(2,i-1) 
         huR = ql(2,i) 
         bL = auxr(1,i-1)
         bR = auxl(1,i)

         hvL = 0.d0
         hvR = 0.d0

         !check for wet/dry boundary
         sE1 =  1.d99
         sE2 = -1.d99
         wall(1) = 1.d0
         wall(2) = 1.d0
         if (hR.le.drytol) then
            sLtest=min(-sqrt(g*hL),huL/hL-sqrt(g*hL)) !what would be the Einfeldt speed of wall problem
            hstartest=hL-(huL/sLtest) !what would be middle state in approx Riemann solution
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
               wall(2)=0.d0
               hR=hL
               huR=-huL
               bR=bL
               phiR=phiL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            sRtest=max(sqrt(g*hR),huR/hR+sqrt(g*hR)) !what would be the Einfeldt speed of wall
            hstartest= hR-(huR/sRtest) !what would be middle state in approx Rimeann solution
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
               wall(1)=0.d0
               hL=hR
               huL=-huR
               bL=bR
               phiL=phiR
            endif
         endif

         if (hR.gt.drytol) then
            uR=huR/hR
            vR=hvR/hR
            phiR = 0.5d0*g*hR**2 + huR**2/hR
         else
            hR = 0.d0
            huR = 0.d0
            uR = 0.d0
            vR = 0.d0
            phiR = 0.d0
            sE2 = max(sE2,huL/hL+2.d0*sqrt(g*hL))
         endif

         if (hL.gt.drytol) then
            uL=huL/hL
            vL=hvL/hL
            phiL = 0.5d0*g*hL**2 + huL**2/hL
         else
            hL=0.d0
            huL=0.d0
            uL=0.d0
            vL=0.d0
            phiL = 0.d0
            sE1 = min(sE1,huR/hR - 2.d0*sqrt(g*hR))
         endif


         !determine wave speeds
         sL=uL-sqrt(g*hL) ! 1 wave speed of left state
         sR=uR+sqrt(g*hR) ! 2 wave speed of right state

         uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
         chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
         sRoe1=uhat-chat ! Roe wave speed 1 wave
         sRoe2=uhat+chat ! Roe wave speed 2 wave

         sE1 = min(sE1,min(sL,sRoe1)) ! Eindfeldt speed 1 wave
         sE2 = max(sE2,max(sR,sRoe2)) ! Eindfeldt speed 2 wave

         !--------------------end initializing...finally----------
         !solve Riemann problem.

         maxiter = 1

         call riemann_aug_JCP(maxiter,3,3,hL,hR,huL,huR,hvL,hvR,bL,bR, &
                              uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

         s(1,i) = sw(1)*wall(1)
         s(2,i) = sw(3)*wall(2)

         do m=1,num_eqn
            fwave(m,1,i)=fw(m,1)*wall(1)
            fwave(m,2,i)=fw(m,3)*wall(2)
            if (sw(2)>0.0) then
               fwave(m,2,i) = fwave(m,2,i) + fw(m,2)*wall(2)
            else
               fwave(m,1,i) = fwave(m,1,i) + fw(m,2)*wall(1)
            endif
         enddo




 30      continue
      enddo


      do i=2-num_ghost,num_cells+num_ghost
         do m=1,num_eqn
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do  mw=1,num_waves
               if (s(mw,i).lt.0.d0) then
                  amdq(m,i) = amdq(m,i) + fwave(m,mw,i)
               elseif (s(mw,i).gt.0.d0) then
                  apdq(m,i) = apdq(m,i) + fwave(m,mw,i)
               else
                  amdq(m,i) = amdq(m,i) + .5d0*fwave(m,mw,i)
                  apdq(m,i) = apdq(m,i) + .5d0*fwave(m,mw,i)
               endif
            enddo
         enddo
      enddo

      return
      end subroutine rp1


!-----------------------------------------------------------------------
      subroutine riemann_aug_JCP(maxiter,num_eqn,num_waves,hL,hR,huL,huR, &
        hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

      ! solve shallow water equations given single left and right states
      ! This solver is described in J. Comput. Phys. (6): 3089-3113, March 2008
      ! Augmented Riemann Solvers for the Shallow Equations with Steady States and Inundation

      ! To use the original solver call with maxiter=1.

      ! This solver allows iteration when maxiter > 1. The iteration seems to help with
      ! instabilities that arise (with any solver) as flow becomes transcritical over variable topo
      ! due to loss of hyperbolicity.



      implicit none

      !input
      integer num_eqn,num_waves,maxiter
      real(kind=8) fw(num_eqn,num_waves)
      real(kind=8) sw(num_waves)
      real(kind=8) hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,sE1,sE2
      real(kind=8) hvL,hvR,vL,vR
      real(kind=8) drytol,g


      !local
      integer m,mw,k,iter
      real(kind=8) A(3,3)
      real(kind=8) r(3,3)
      real(kind=8) lambda(3)
      real(kind=8) del(3)
      real(kind=8) beta(3)

      real(kind=8) delh,delhu,delphi,delb,delnorm
      real(kind=8) rare1st,rare2st,sdelta,raremin,raremax
      real(kind=8) criticaltol,convergencetol,raretol
      real(kind=8) s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
      real(kind=8) huRstar,huLstar,uRstar,uLstar,hstarHLL
      real(kind=8) deldelh,deldelphi
      real(kind=8) s1m,s2m,hm
      real(kind=8) det1,det2,det3,determinant

      logical rare1,rare2,rarecorrector,rarecorrectortest,sonic

      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL
      delb = bR-bL
      delnorm = delh**2 + delphi**2

      call riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,1,drytol,g)


      lambda(1)= min(sE1,s2m) !Modified Einfeldt speed
      lambda(3)= max(sE2,s1m) !Modified Eindfeldt speed
      sE1=lambda(1)
      sE2=lambda(3)
      lambda(2) = 0.d0  ! ### Fix to avoid uninitialized value in loop on mw -- Correct?? ###

      
      hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve

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
            if (max(rare1st,rare2st).gt.raremin*sdelta.and. & 
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
      do mw=1,num_waves
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
      criticaltol = 1.d-6
      deldelh = -delb
      deldelphi = -g*0.5d0*(hR+hL)*delb

!     !determine a few quanitites needed for steady state wave if iterated
      hLstar=hL
      hRstar=hR
      uLstar=uL
      uRstar=uR
      huLstar=uLstar*hLstar
      huRstar=uRstar*hRstar

      !iterate to better determine the steady state wave
      convergencetol=1.d-6
      do iter=1,maxiter
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
         sonic=.false.
         if (abs(s1s2bar).le.criticaltol) sonic=.true.
         if (s1s2bar*s1s2tilde.le.criticaltol) sonic=.true.
         if (s1s2bar*sE1*sE2.le.criticaltol) sonic = .true.
         if (min(abs(sE1),abs(sE2)).lt.criticaltol) sonic=.true.
         if (sE1.lt.0.d0.and.s1m.gt.0.d0) sonic = .true.
         if (sE2.gt.0.d0.and.s2m.lt.0.d0) sonic = .true.
         if ((uL+sqrt(g*hL))*(uR+sqrt(g*hR)).lt.0.d0) sonic=.true.
         if ((uL-sqrt(g*hL))*(uR-sqrt(g*hR)).lt.0.d0) sonic=.true.

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
         do mw=1,num_waves
            if (lambda(mw).lt.0.d0) then
               hLstar= hLstar + beta(mw)*r(1,mw)
               huLstar= huLstar + beta(mw)*r(2,mw)
            endif
         enddo
         do mw=num_waves,1,-1
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

      do mw=1,num_waves
         sw(mw)=lambda(mw)
         fw(1,mw)=beta(mw)*r(2,mw)
         fw(2,mw)=beta(mw)*r(3,mw)
         fw(3,mw)=beta(mw)*r(2,mw)
      enddo
      !find transverse components (ie huv jumps).
      fw(3,1)=fw(3,1)*vL
      fw(3,3)=fw(3,3)*vR
      fw(3,2)= hR*uR*vR - hL*uL*vL - fw(3,1)- fw(3,3)

      return

      end !subroutine riemann_aug_JCP-------------------------------------------------


!-----------------------------------------------------------------------
      subroutine riemann_ssqfwave(maxiter,num_eqn,num_waves,hL,hR,huL,huR, &
          hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

      ! solve shallow water equations given single left and right states
      ! steady state wave is subtracted from delta [q,f]^T before decomposition

      implicit none

      !input
      integer num_eqn,num_waves,maxiter

      real(kind=8) hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,sE1,sE2
      real(kind=8) vL,vR,hvL,hvR
      real(kind=8) drytol,g

      !local
      integer iter

      logical sonic

      real(kind=8) delh,delhu,delphi,delb,delhdecomp,delphidecomp
      real(kind=8) s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
      real(kind=8) uRstar,uLstar,hstarHLL
      real(kind=8) deldelh,deldelphi
      real(kind=8) alpha1,alpha2,beta1,beta2,delalpha1,delalpha2
      real(kind=8) criticaltol,convergencetol
      real(kind=8) sL,sR
      real(kind=8) uhat,chat,sRoe1,sRoe2

      real(kind=8) sw(num_waves)
      real(kind=8) fw(num_eqn,num_waves)

      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL
      delb = bR-bL

      convergencetol= 1.d-16
      criticaltol = 1.d-99

      deldelh = -delb
      deldelphi = -g*0.5d0*(hR+hL)*delb

!     !if no source term, skip determining steady state wave
      if (abs(delb).gt.0.d0) then
!
         !determine a few quanitites needed for steady state wave if iterated
         hLstar=hL
         hRstar=hR
         uLstar=uL
         uRstar=uR
         hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve

         alpha1=0.d0
         alpha2=0.d0

!        !iterate to better determine Riemann problem
         do iter=1,maxiter

            !determine steady state wave (this will be subtracted from the delta vectors)
            hbar =  max(0.5d0*(hLstar+hRstar),0.d0)
            s1s2bar = 0.25d0*(uLstar+uRstar)**2 - g*hbar
            s1s2tilde= max(0.d0,uLstar*uRstar) - g*hbar


!           !find if sonic problem
            sonic=.false.
            if (abs(s1s2bar).le.criticaltol) sonic=.true.
            if (s1s2bar*s1s2tilde.le.criticaltol) sonic=.true.
            if (s1s2bar*sE1*sE2.le.criticaltol) sonic = .true.
            if (min(abs(sE1),abs(sE2)).lt.criticaltol) sonic=.true.

!           !find jump in h, deldelh
            if (sonic) then
               deldelh =  -delb
            else
               deldelh = delb*g*hbar/s1s2bar
            endif
!           !bounds in case of critical state resonance, or negative states
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

!           !find jump in phi, deldelphi
            if (sonic) then
               deldelphi = -g*hbar*delb
            else
               deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar
            endif
!           !bounds in case of critical state resonance, or negative states
            deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
            deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))

!---------determine fwaves ------------------------------------------

!           !first decomposition
            delhdecomp = delh-deldelh
            delalpha1 = (sE2*delhdecomp - delhu)/(sE2-sE1)-alpha1
            alpha1 = alpha1 + delalpha1
            delalpha2 = (delhu - sE1*delhdecomp)/(sE2-sE1)-alpha2
            alpha2 = alpha2 + delalpha2

            !second decomposition
            delphidecomp = delphi - deldelphi
            beta1 = (sE2*delhu - delphidecomp)/(sE2-sE1)
            beta2 = (delphidecomp - sE1*delhu)/(sE2-sE1)

            if ((delalpha2**2+delalpha1**2).lt.convergencetol**2) then
               exit
            endif
!
            if (sE2.gt.0.d0.and.sE1.lt.0.d0) then
               hLstar=hL+alpha1
               hRstar=hR-alpha2
!               hustar=huL+alpha1*sE1
               hustar = huL + beta1
            elseif (sE1.ge.0.d0) then
               hLstar=hL
               hustar=huL
               hRstar=hR - alpha1 - alpha2
            elseif (sE2.le.0.d0) then
               hRstar=hR
               hustar=huR
               hLstar=hL + alpha1 + alpha2
            endif
!
            if (hLstar.gt.drytol) then
               uLstar=hustar/hLstar
            else
               hLstar=max(hLstar,0.d0)
               uLstar=0.d0
            endif
!
            if (hRstar.gt.drytol) then
               uRstar=hustar/hRstar
            else
               hRstar=max(hRstar,0.d0)
               uRstar=0.d0
            endif

         enddo
      endif

      delhdecomp = delh - deldelh
      delphidecomp = delphi - deldelphi

      !first decomposition
      alpha1 = (sE2*delhdecomp - delhu)/(sE2-sE1)
      alpha2 = (delhu - sE1*delhdecomp)/(sE2-sE1)

      !second decomposition
      beta1 = (sE2*delhu - delphidecomp)/(sE2-sE1)
      beta2 = (delphidecomp - sE1*delhu)/(sE2-sE1)

      ! 1st nonlinear wave
      fw(1,1) = alpha1*sE1
      fw(2,1) = beta1*sE1
      fw(3,1) = fw(1,1)*vL
      ! 2nd nonlinear wave
      fw(1,3) = alpha2*sE2
      fw(2,3) = beta2*sE2
      fw(3,3) = fw(1,3)*vR
      ! advection of transverse wave
      fw(1,2) = 0.d0
      fw(2,2) = 0.d0
      fw(3,2) = hR*uR*vR - hL*uL*vL -fw(3,1)-fw(3,3)
      !speeds
      sw(1)=sE1
      sw(2)=0.5d0*(sE1+sE2)
      sw(3)=sE2

      return

      end subroutine !-------------------------------------------------


!-----------------------------------------------------------------------
      subroutine riemann_fwave(num_eqn,num_waves,hL,hR,huL,huR,hvL,hvR, &
                  bL,bR,uL,uR,vL,vR,phiL,phiR,s1,s2,drytol,g,sw,fw)

      ! solve shallow water equations given single left and right states
      ! solution has two waves.
      ! flux - source is decomposed.

      implicit none

      !input
      integer num_eqn,num_waves

      real(kind=8) hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,s1,s2
      real(kind=8) hvL,hvR,vL,vR
      real(kind=8) drytol,g

      real(kind=8) sw(num_waves)
      real(kind=8) fw(num_eqn,num_waves)

      !local
      real(kind=8) delh,delhu,delphi,delb,delhdecomp,delphidecomp
      real(kind=8) deldelh,deldelphi
      real(kind=8) beta1,beta2


      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL
      delb = bR-bL

      deldelphi = -g*0.5d0*(hR+hL)*delb
      delphidecomp = delphi - deldelphi

      !flux decomposition
      beta1 = (s2*delhu - delphidecomp)/(s2-s1)
      beta2 = (delphidecomp - s1*delhu)/(s2-s1)

      sw(1)=s1
      sw(2)=0.5d0*(s1+s2)
      sw(3)=s2
      ! 1st nonlinear wave
      fw(1,1) = beta1
      fw(2,1) = beta1*s1
      fw(3,1) = beta1*vL
      ! 2nd nonlinear wave
      fw(1,3) = beta2
      fw(2,3) = beta2*s2
      fw(3,3) = beta2*vR
      ! advection of transverse wave
      fw(1,2) = 0.d0
      fw(2,2) = 0.d0
      fw(3,2) = hR*uR*vR - hL*uL*vL -fw(3,1)-fw(3,3)
      return

      end !subroutine -------------------------------------------------





!=============================================================================
      subroutine riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2, &
                   maxiter,drytol,g)

      !determine the Riemann structure (wave-type in each family)


      implicit none

      !input
      real(kind=8) hL,hR,uL,uR,drytol,g
      integer maxiter

      !output
      real(kind=8) s1m,s2m
      logical rare1,rare2

      !local
      real(kind=8) hm,u1m,u2m,um,delu
      real(kind=8) h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
      integer iter



!     !Test for Riemann structure

      h_min=min(hR,hL)
      h_max=max(hR,hL)
      delu=uR-uL

      if (h_min.le.drytol) then
         hm=0.d0
         um=0.d0
         s1m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
         s2m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
         if (hL.le.0.d0) then
            rare2=.true.
            rare1=.false.
         else
            rare1=.true.
            rare2=.false.
         endif

      else
         F_min= delu+2.d0*(sqrt(g*h_min)-sqrt(g*h_max))
         F_max= delu + (h_max-h_min)*(sqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

         if (F_min.gt.0.d0) then !2-rarefactions

            hm=(1.d0/(16.d0*g))* &
                     max(0.d0,-delu+2.d0*(sqrt(g*hL)+sqrt(g*hR)))**2
            um=sign(1.d0,hm)*(uL+2.d0*(sqrt(g*hL)-sqrt(g*hm)))

            s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
            s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)

            rare1=.true.
            rare2=.true.

         elseif (F_max.le.0.d0) then !2 shocks

!           !root finding using a Newton iteration on sqrt(h)===
            h0=h_max
            do iter=1,maxiter
               gL=sqrt(.5d0*g*(1/h0 + 1/hL))
               gR=sqrt(.5d0*g*(1/h0 + 1/hR))
               F0=delu+(h0-hL)*gL + (h0-hR)*gR
               dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
               slope=2.d0*sqrt(h0)*dfdh
               h0=(sqrt(h0)-F0/slope)**2
            enddo
               hm=h0
               u1m=uL-(hm-hL)*sqrt((.5d0*g)*(1/hm + 1/hL))
               u2m=uR+(hm-hR)*sqrt((.5d0*g)*(1/hm + 1/hR))
               um=.5d0*(u1m+u2m)

               s1m=u1m-sqrt(g*hm)
               s2m=u2m+sqrt(g*hm)
               rare1=.false.
               rare2=.false.

         else !one shock one rarefaction
            h0=h_min

            do iter=1,maxiter
               F0=delu + 2.d0*(sqrt(g*h0)-sqrt(g*h_max)) &
                       + (h0-h_min)*sqrt(.5d0*g*(1/h0+1/h_min))
               slope=(F_max-F0)/(h_max-h_min)
               h0=h0-F0/slope
            enddo

            hm=h0
            if (hL.gt.hR) then
               um=uL+2.d0*sqrt(g*hL)-2.d0*sqrt(g*hm)
               s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
               s2m=uL+2.d0*sqrt(g*hL)-sqrt(g*hm)
               rare1=.true.
               rare2=.false.
            else
               s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)
               s1m=uR-2.d0*sqrt(g*hR)+sqrt(g*hm)
               um=uR-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hm)
               rare2=.true.
               rare1=.false.
            endif
         endif
      endif

      return

      end ! subroutine riemanntype----------------------------------------------------------------
