! =====================================================
      subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,
     &                ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
!
!     Riemann solver in the transverse direction using 
!     Jacobian matrix from left cell (if imp==1) or right cell (if imp==2).
!
!     Note this has been modified from the version used in v5.7.x and
!     earlier, where Roe averages based on ql and qr were used, which
!     is not correct.  In addition:
!      - a bug in the second component of the eigenvectors was fixed.
!      - when s(2) is close to zero this component of flux difference
!        is split equally between bmasdq and bpasdq to improve symmetry.
!   
!     Further modified to clean up and avoid a lot of work in dry cells.

!-----------------------last modified October 2020 ----------------------

      use geoclaw_module, only: g => grav, tol => dry_tolerance
      use geoclaw_module, only: coordinate_system,earth_radius,deg2rad

      implicit none

      integer, intent(in) :: ixy,maxm,meqn,maux,mwaves,mbc,mx,imp

      real(kind=8), intent(in) ::  ql(meqn,1-mbc:maxm+mbc)
      real(kind=8), intent(in) ::  qr(meqn,1-mbc:maxm+mbc)
      real(kind=8), intent(in) ::  asdq(meqn,1-mbc:maxm+mbc)
      real(kind=8), intent(in) ::  aux1(maux,1-mbc:maxm+mbc)
      real(kind=8), intent(in) ::  aux2(maux,1-mbc:maxm+mbc)
      real(kind=8), intent(in) ::  aux3(maux,1-mbc:maxm+mbc)

      real(kind=8), intent(out) ::  bmasdq(meqn,1-mbc:maxm+mbc)
      real(kind=8), intent(out) ::  bpasdq(meqn,1-mbc:maxm+mbc)

      ! local:
      real(kind=8) ::  s(mwaves), r(meqn,mwaves), beta(mwaves)
      real(kind=8) ::  h,u,v
      real(kind=8) ::  delf1,delf2,delf3
      real(kind=8) ::  dxdcm,dxdcp,topo1,topo3,eta

      integer :: i,mw,mu,mv
      

      if (ixy == 1) then
         ! normal solve was in x-direction
         mu = 2
         mv = 3
      else
         ! normal solve was in y-direction
         mu = 3
         mv = 2
      endif

      ! initialize all components of result to 0:
      bmasdq(:,:) = 0.d0
      bpasdq(:,:) = 0.d0


      do i=2-mbc,mx+mbc

         if (imp==1) then
            h = qr(1,i-1)
         else
            h = ql(1,i)
         endif

         if (h <= tol) then
             ! fluctuation going into a dry cell, don't know how to split,
             ! so leave bmadsq(:,i)=bpasdq(:,i)=0 and go on to next i:
             cycle  
         endif

         ! compute velocities in relevant cell, and other quantities:

         if (imp==1) then
              ! fluctuation being split is left-going
              u = qr(mu,i-1)/h
              v = qr(mv,i-1)/h
              eta = h + aux2(1,i-1)
              topo1 = aux1(1,i-1)
              topo3 = aux3(1,i-1)
         else
              ! fluctuation being split is right-going
              u = ql(mu,i)/h
              v = ql(mv,i)/h
              eta = h + aux2(1,i)
              topo1 = aux1(1,i)
              topo3 = aux3(1,i)
         endif

         ! check if cell that transverse waves go into are both too high:
         ! Note: prior to v5.8.0 this checked against max rather than min
         if (eta < min(topo1,topo3)) cycle  ! go to next i

         ! if we get here, we want to do the splitting (no dry cells),
         ! so compute the necessary quantities:

         if (coordinate_system == 2) then
            ! on the sphere:
            if (ixy == 2) then
               dxdcp=(earth_radius*deg2rad)
               dxdcm = dxdcp
            else
               if (imp == 1) then
                  dxdcp = earth_radius*cos(aux3(3,i-1))*deg2rad
                  dxdcm = earth_radius*cos(aux1(3,i-1))*deg2rad
               else
                  dxdcp = earth_radius*cos(aux3(3,i))*deg2rad
                  dxdcm = earth_radius*cos(aux1(3,i))*deg2rad
               endif
            endif
         else
            ! coordinate_system == 1 means Cartesian:
            dxdcp = 1.d0
            dxdcm = 1.d0
         endif

c        Determine some speeds necessary for the Jacobian
            
         ! In v5.7.x and prior versions,
         ! we used left right states to define Roe averages,
         ! which is consistent with those used in rpn2.
         ! But now we are computing upgoing, downgoing waves either in
         ! cell on left (if imp==1) or on right (if imp==2) so we
         ! should possibly use q values in cells above/below,
         ! but these aren't available (only aux values).
         ! At any rate, there is no clear justification for using cells
         ! on the other side of the normal-solve interface.

         ! v5.8.0: modified to use left or right state alone in defining
         ! Jacobian, based on imp:

         s(1) = v - dsqrt(g*h)
         s(2) = v
         s(3) = v + dsqrt(g*h)

c        Determine asdq decomposition (beta)

         delf1 = asdq(1,i)
         delf2 = asdq(mu,i)
         delf3 = asdq(mv, i)

         ! v5.8.0: fixed bug in beta(2): u in place of s(2)=v
         beta(1) = (s(3)*delf1 - delf3) / (s(3) - s(1))
         beta(2) = -u*delf1 + delf2
         beta(3) = (delf3 - s(1)*delf1) / (s(3) - s(1))

c        Set-up eigenvectors
         r(1,1) = 1.d0
         r(2,1) = u    ! v5.8.0: fixed bug, u not s(2)=v
         r(3,1) = s(1)

         r(1,2) = 0.d0
         r(2,2) = 1.d0
         r(3,2) = 0.d0

         r(1,3) = 1.d0
         r(2,3) = u    ! v5.8.0: fixed bug, u not s(2)=v
         r(3,3) = s(3)


         ! compute fluctuations

         do  mw=1,3
            if ((s(mw) < 0.d0) .and. (eta >= topo1)) then
                 bmasdq(1,i) =bmasdq(1,i) + dxdcm*s(mw)*beta(mw)*r(1,mw)
                 bmasdq(mu,i)=bmasdq(mu,i)+ dxdcm*s(mw)*beta(mw)*r(2,mw)
                 bmasdq(mv,i)=bmasdq(mv,i)+ dxdcm*s(mw)*beta(mw)*r(3,mw)
            elseif ((s(mw) > 0.d0) .and. (eta >= topo3)) then
                 bpasdq(1,i) =bpasdq(1,i) + dxdcp*s(mw)*beta(mw)*r(1,mw)
                 bpasdq(mu,i)=bpasdq(mu,i)+ dxdcp*s(mw)*beta(mw)*r(2,mw)
                 bpasdq(mv,i)=bpasdq(mv,i)+ dxdcp*s(mw)*beta(mw)*r(3,mw)
            endif
         enddo  ! loop on mw


      enddo  ! loop on i

      return
      end
