! =====================================================
      subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,
     &                ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
      use geoclaw_module, only: g => grav, tol => dry_tolerance
      use geoclaw_module, only: coordinate_system,earth_radius,deg2rad

      implicit none
!
!     # Riemann solver in the transverse direction.

!
!     Aux1 is below, aux2 is where you are, aux3 is above

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water adjoint               !
!             equations. It is modified from rpt2_geoclaw.f.                !
!                                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer ixy,maxm,meqn,maux,mwaves,mbc,mx,imp

      double precision  ql(meqn,1-mbc:maxm+mbc)
      double precision  qr(meqn,1-mbc:maxm+mbc)
      double precision  asdq(meqn,1-mbc:maxm+mbc)
      double precision  bmasdq(meqn,1-mbc:maxm+mbc)
      double precision  bpasdq(meqn,1-mbc:maxm+mbc)
      double precision  aux1(maux,1-mbc:maxm+mbc)
      double precision  aux2(maux,1-mbc:maxm+mbc)
      double precision  aux3(maux,1-mbc:maxm+mbc)

      double precision  s(3)
      double precision  r(3,3)
      double precision  a(3)
      double precision  abs_tol
      double precision  hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
      double precision  hLhat, hRhat, hHatm, hHatp, hHat
      double precision  uhat,vhat,roe1,roe3,sm,sp,s1l,s3r,ss
      double precision  delf1,delf2,delf3,dxdcd,dxdcu
      double precision  dxdcm,dxdcp,topo1,topo3,eta

      integer i,m,mw,mu,mv, i1

      abs_tol=tol

      if (ixy.eq.1) then
        mu = 2
        mv = 3
      else
        mu = 3
        mv = 2
      endif

      do i=2-mbc,mx+mbc-1

         hl=qr(1,i-1)
         hr=ql(1,i)

!===========determine velocity from momentum===========================
       if (imp == 1) then
!            # asdq = amdq, moving to left
           i1 = i-1
       else
!            # asdq = apdq, moving to right
           i1 = i
       endif

       do mw=1,mwaves
           s(mw)=0.d0
           a(mw)=0.d0
           do m=1,meqn
               r(m,mw)=0.d0
           enddo
       enddo

       if (aux2(1,i1) >= 0.d0) go to 90

       if (hl.lt.abs_tol) then
           hl=0.d0
           hHatm = 0.d0
       else
           hHatm = -aux2(1,i1)
       endif

       if (hr.lt.abs_tol) then
           hr=0.d0
           hHatp = 0.d0
       else
           hHatp = -aux2(1,i1)
       endif

       hHat = (hHatp + hHatm)/2

       if (hl <= tol .and. hr <= tol) go to 90

*      !check and see if cell that transverse waves are going in is high and dry
       if (imp.eq.1) then
          eta = qr(1,i-1)  + aux2(1,i-1)
          topo1 = aux1(1,i-1)
          topo3 = aux3(1,i-1)
       else
          eta = ql(1,i) + aux2(1,i)
          topo1 = aux1(1,i)
          topo3 = aux3(1,i)
       endif

       if (eta.lt.max(topo1,topo3)) go to 90


c=======================Determine asdq decomposition (beta)============

         sm=dsqrt(g*hHatm)
         sp=dsqrt(g*hHatp)
         ss=dsqrt(g*hHat)

         delf1=asdq(1,i)
         delf2=asdq(mu,i)
         delf3=asdq(mv, i)

         a(1) = (delf1 + ss*delf3)/ (ss+sm)
         a(2) = (-delf1 + ss*delf3)/(sp + ss)
c======================End =================================================

c=====================Set-up eigenvectors===================================
         r(1,1) = sm
         r(3,1) = 1.d0
         r(2,1) = 0.d0

         r(1,2) = -sp
         r(3,2) = 1.d0
         r(2,2) = 0.d0

         s(1)=-sm
         s(2)=sp
         s(3)=0.d0
c============================================================================
90      continue
c============= compute fluctuations==========================================

            dxdcm = 1.d0
            dxdcp = 1.d0

       if (coordinate_system.eq.2) then
         if (ixy.eq.2) then
           dxdcp=(earth_radius*deg2rad)
           dxdcm = dxdcp
         else
           if (imp.eq.1) then
             dxdcp = earth_radius*cos(aux3(3,i-1))*deg2rad
             dxdcm = earth_radius*cos(aux1(3,i-1))*deg2rad
           else
             dxdcp = earth_radius*cos(aux3(3,i))*deg2rad
             dxdcm = earth_radius*cos(aux1(3,i))*deg2rad
           endif
         endif
       endif

            bmasdq(1,i)=0.0d0
            bpasdq(1,i)=0.0d0
            bmasdq(2,i)=0.0d0
            bpasdq(2,i)=0.0d0
            bmasdq(3,i)=0.0d0
            bpasdq(3,i)=0.0d0
            do  mw=1,3
              if (s(mw).lt.0.d0) then
                bmasdq(1,i) =bmasdq(1,i) + dxdcm*s(mw)*a(mw)*r(1,mw)
                bmasdq(mu,i)=bmasdq(mu,i)+ dxdcm*s(mw)*a(mw)*r(2,mw)
                bmasdq(mv,i)=bmasdq(mv,i)+ dxdcm*s(mw)*a(mw)*r(3,mw)
              elseif (s(mw).gt.0.d0) then
                bpasdq(1,i) =bpasdq(1,i) + dxdcp*s(mw)*a(mw)*r(1,mw)
                bpasdq(mu,i)=bpasdq(mu,i)+ dxdcp*s(mw)*a(mw)*r(2,mw)
                bpasdq(mv,i)=bpasdq(mv,i)+ dxdcp*s(mw)*a(mw)*r(3,mw)
              endif
            enddo

c========================================================================
         enddo
      return
      end
