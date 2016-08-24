c======================================================================
       subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,
     &                 ql,qr,auxl,auxr,fwave,s,amdq,apdq)
c======================================================================
c
c Solves normal Riemann problems for the 2D SHALLOW WATER equations
c     with topography:
c     #        h_t + (hu)_x + (hv)_y = 0                           #
c     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
c     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #

c On input, ql contains the state vector at the left edge of each cell
c     qr contains the state vector at the right edge of each cell
c
c This data is along a slice in the x-direction if ixy=1
c     or the y-direction if ixy=2.

c  Note that the i'th Riemann problem has left state qr(i-1,:)
c     and right state ql(i,:)
c  From the basic clawpack routines, this routine is called with
c     ql = qr
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water adjoint               !
!             equations. It is modified from rpn2_geoclaw.f.                !
!                                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use geoclaw_module, only: g => grav, drytol => dry_tolerance
      use geoclaw_module, only: earth_radius, deg2rad
      use amr_module, only: mcapa

      implicit none

      !input
      integer maxm,meqn,maux,mwaves,mbc,mx,ixy

      double precision  fwave(meqn, mwaves, 1-mbc:maxm+mbc)
      double precision  s(mwaves, 1-mbc:maxm+mbc)
      double precision  ql(meqn, 1-mbc:maxm+mbc)
      double precision  qr(meqn, 1-mbc:maxm+mbc)
      double precision  apdq(meqn,1-mbc:maxm+mbc)
      double precision  amdq(meqn,1-mbc:maxm+mbc)
      double precision  auxl(maux,1-mbc:maxm+mbc)
      double precision  auxr(maux,1-mbc:maxm+mbc)

      !local only
      integer m,i,mw,maxiter,mu,nv
      double precision wall(3)
      double precision fw(3,3)
      double precision sw(3)
      double precision delta(3)

      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
      double precision bR,bL,cL,cR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision s1m,s2m
      double precision hstar,hstartest,hstarHLL,sLtest,sRtest
      double precision tw,dxdc
      double precision beta1, beta2, beta3, hRhat, hLhat

      logical rare1,rare2

      !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

!-----------------------Initializing-----------------------------------
         !inform of a bad riemann problem from the start
         if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
         endif

        !Initialize Riemann problem for grid interface
        do mw=1,mwaves
            s(mw,i)=0.d0
            fwave(1,mw,i)=0.d0
            fwave(2,mw,i)=0.d0
            fwave(3,mw,i)=0.d0

            fw(1,mw)=0.d0
            fw(2,mw)=0.d0
            fw(3,mw)=0.d0
        enddo

         !set normal direction
         if (ixy.eq.1) then
            mu=2
            nv=3
         else
            mu=3
            nv=2
         endif

         !zero (small) negative values if they exist
         if (qr(1,i-1).lt.0.d0) then
           qr(1,i-1)=0.d0
           qr(2,i-1)=0.d0
           qr(3,i-1)=0.d0
         endif

        if (ql(1,i).lt.0.d0) then
          ql(1,i)=0.d0
          ql(2,i)=0.d0
          ql(3,i)=0.d0
        endif

         !skip problem if in a completely dry area
         if (qr(1,i-1) <= drytol .and. ql(1,i) <= drytol) then
             go to 30
         endif

        !Riemann problem variables
        hL = qr(1,i-1)
        hR = ql(1,i)
        bL = auxr(1,i-1)
        bR = auxl(1,i)

        huL = qr(mu,i-1)
        huR = ql(mu,i)
        hvL=qr(nv,i-1)
        hvR=ql(nv,i)

        ! Linearizing
        hLhat = -bL
        hRhat = -bR

c       ! Check for wet/dry boundary
        if (hR <= drytol) then
            hR = 0.d0
            huR = 0.d0
            hvR = 0.d0
            uR = 0.d0
            vR = 0.d0
            hRhat = 0.d0
        else
            uR = huR/hR
            vR = hvR/hR
        endif
        if (hL <= drytol) then
            hL = 0.d0
            huL = 0.d0
            hvL = 0.d0
            uL = 0.d0
            vL = 0.d0
            hLhat = 0.d0
        else
            uL = huL/hL
            vL = hvL/hL
        endif

        wall(1) = 1.d0
        wall(2) = 1.d0
        wall(3) = 1.d0

        if (hR.le.drytol) then
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,g)
            hstartest=max(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values
c                that mirror left for wall problem
                wall(2)=0.d0
                wall(3)=0.d0
                hR=hL
                huR=-huL
                bR=bL
                uR=-uL
                vR=vL
                hRhat = hLhat  ! added since this is a linear problem
             elseif (hL+bL.lt.bR) then
                bR=hL+bL
            endif
        elseif (hL.le.drytol) then ! right surface is lower than left topo
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,g)
            hstartest=max(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values
c               that mirror right
                wall(1)=0.d0
                wall(2)=0.d0
                hL=hR
                huL=-huR
                bL=bR
                uL=-uR
                vL=vR
                hLhat = hRhat  ! added since this is a linear problem
            elseif (hR+bR.lt.bL) then
                bL=hR+bR
            endif
        endif

         !determine wave speeds
         cL=sqrt(g*hLhat) ! 1 wave speed of left state
         cR=sqrt(g*hRhat) ! 2 wave speed of right state

         !--------------------end initializing----------
         !solve Riemann problem.

c       # f-wave splitting
        delta(1) = -huR*cR**2 + huL*cL**2
        delta(2) = -(hR+bR) + (hL+bL)
        delta(3) = 0

        beta1 = (delta(1) + cR*delta(2))/(cL + cR)
        beta2 = (-delta(1) + cL*delta(2))/(cL + cR)
        beta3 = delta(3)

        if (cL + cR == 0.d0) then
            beta1 = (delta(1) + cR*delta(2))
            beta2 = (-delta(1) + cL*delta(2))
        endif

c       # Compute the waves.
        fw(1,1) = cL*beta1
        fw(mu,1) = beta1
        fw(nv,1) = 0.d0
        sw(1) = -cL

        fw(1,2) = -cR*beta2
        fw(mu,2) = beta2
        fw(nv,2) = 0.d0
        sw(2) = cR

c        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)

            fw(1,mw)=fw(1,mw)*wall(mw)
            fw(2,mw)=fw(2,mw)*wall(mw)
            fw(3,mw)=fw(3,mw)*wall(mw)
         enddo

         do mw=1,mwaves
            s(mw,i)=sw(mw)
            fwave(1,mw,i)=fw(1,mw)
            fwave(mu,mw,i)=fw(mu,mw)
            fwave(nv,mw,i)=fw(nv,mw)
         enddo

 30      continue
      enddo

c==========Capacity for mapping from latitude longitude to physical space====
        if (mcapa.gt.0) then
          do i=2-mbc,mx+mbc
            if (ixy.eq.1) then
              dxdc=(earth_radius*deg2rad)
            else
              dxdc=earth_radius*cos(auxl(3,i))*deg2rad
            endif

            do mw=1,mwaves
              s(mw,i)=dxdc*s(mw,i)
              fwave(1,mw,i)=dxdc*fwave(1,mw,i)
              fwave(2,mw,i)=dxdc*fwave(2,mw,i)
              fwave(3,mw,i)=dxdc*fwave(3,mw,i)
            enddo
          enddo
        endif

c============= compute fluctuations=============================================
         amdq(1:3,:) = 0.d0
         apdq(1:3,:) = 0.d0
         do i=2-mbc,mx+mbc
             do  mw=1,mwaves
               if (s(mw,i) < 0.d0) then
                 amdq(1:3,i) = amdq(1:3,i) + fwave(1:3,mw,i)
               else if (s(mw,i) > 0.d0) then
                 apdq(1:3,i)  = apdq(1:3,i) + fwave(1:3,mw,i)
               else
                 amdq(1:3,i) = amdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
                 apdq(1:3,i) = apdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
               endif
             enddo
         enddo

         do i=2-mbc, mx+mbc
             amdq(1:meqn,i) = fwave(1:meqn,1,i)
             apdq(1:meqn,i) = fwave(1:meqn,2,i)
         enddo

      return
      end subroutine
