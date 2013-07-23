!     ==================================================================
      subroutine rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr, &
                      auxl,auxr,maux,wave,s,amdq,apdq)
!     ==================================================================
!
!     # Roe-solver for the Euler equations
!      -----------------------------------------------------------
!
!     # solve Riemann problems along one slice of data.
!     # This data is along a slice in the x-direction if ixyz=1
!     #                               the y-direction if ixyz=2.
!     #                               the z-direction if ixyz=3.
!
!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell
!
!     # On output, wave contains the waves, s the speeds,
!     # and amdq, apdq the left-going and right-going flux differences,
!     # respectively.  
!
!     # Note that the i'th Riemann problem has left state qr(:,i-1)
!     #                                    and right state ql(:,i)
!     # From the basic clawpack routines, this routine is called with ql = qr
!
!
      implicit real*8(a-h,o-z)
!
      dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension amdq(meqn, 1-mbc:maxm+mbc)
      dimension apdq(meqn, 1-mbc:maxm+mbc)
      dimension auxl(maux, 1-mbc:maxm+mbc)
      dimension auxr(maux, 1-mbc:maxm+mbc)
!
      parameter (maxmrp = 1002)
      dimension delta(5)
      logical efix

      double precision u2v2w2(-1:maxmrp), &
            u(-1:maxmrp),v(-1:maxmrp),w(-1:maxmrp),enth(-1:maxmrp), &
            a(-1:maxmrp),g1a2(-1:maxmrp),euv(-1:maxmrp)
!
      common /cparam/ gamma
!
      data efix /.true./    !# use entropy fix for transonic rarefactions
!
      gamma1 = gamma - 1.d0

      if (-1.gt.1-mbc .or. maxmrp .lt. maxm+mbc) then
         write(6,*) 'need to increase maxmrp in rpn3eu'
         write(6,*) 'maxm+mbc=',maxm+mbc
         stop
         endif

      if (mwaves .ne. 3) then
         write(6,*) '*** Should have mwaves=3 for this Riemann solver'
         stop
         endif
    
!
!     # set mu to point to  the component of the system that corresponds
!     # to momentum in the direction of this slice, mv and mw to the
!     # orthogonal momentums:
!
      if(ixyz .eq. 1)then
          mu = 2
          mv = 3
          mw = 4
      else if(ixyz .eq. 2)then
          mu = 3
          mv = 4
          mw = 2
      else
          mu = 4
          mv = 2
          mw = 3
      endif
!
!
!     # note that notation for u,v, and w reflects assumption that the
!     # Riemann problems are in the x-direction with u in the normal
!     # direction and v and w in the orthogonal directions, but with the
!     # above definitions of mu, mv, and mw the routine also works with
!     # ixyz=2 and ixyz = 3
!     # and returns, for example, f0 as the Godunov flux g0 for the
!     # Riemann problems u_t + g(u)_y = 0 in the y-direction.
!
!
!     # Compute the Roe-averaged variables needed in the Roe solver.
!
      do 10 i = 2-mbc, mx+mbc
         if (qr(1,i-1) .le. 0.d0 .or. ql(1,i) .le. 0.d0) then
            write(*,*) i, mu, mv, mw
            write(*,990) (qr(j,i-1),j=1,5)
            write(*,990) (ql(j,i),j=1,5)
 990        format(5e12.4)
             if (ixyz .eq. 1) &
               write(6,*) '*** rho .le. 0 in x-sweep at ',i
             if (ixyz .eq. 2) &
               write(6,*) '*** rho .le. 0 in y-sweep at ',i
             if (ixyz .eq. 3) &
               write(6,*) '*** rho .le. 0 in z-sweep at ',i
             write(6,*) 'stopped with rho < 0...'
             stop
             endif
         rhsqrtl = dsqrt(qr(1,i-1))
         rhsqrtr = dsqrt(ql(1,i))
         pl = gamma1*(qr(5,i-1) - 0.5d0*(qr(mu,i-1)**2 + &
                qr(mv,i-1)**2 + qr(mw,i-1)**2)/qr(1,i-1))
         pr = gamma1*(ql(5,i) - 0.5d0*(ql(mu,i)**2 + &
                ql(mv,i)**2 + ql(mw,i)**2)/ql(1,i))
         rhsq2 = rhsqrtl + rhsqrtr
         u(i) = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
         v(i) = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
         w(i) = (qr(mw,i-1)/rhsqrtl + ql(mw,i)/rhsqrtr) / rhsq2
         enth(i) = (((qr(5,i-1)+pl)/rhsqrtl &
                    + (ql(5,i)+pr)/rhsqrtr)) / rhsq2
         u2v2w2(i) = u(i)**2 + v(i)**2 + w(i)**2
         a2 = gamma1*(enth(i) - .5d0*u2v2w2(i))
         if (a2 .le. 0.d0) then
             if (ixyz .eq. 1) &
               write(6,*) '*** a2 .le. 0 in x-sweep at ',i
             if (ixyz .eq. 2) &
               write(6,*) '*** a2 .le. 0 in y-sweep at ',i
             if (ixyz .eq. 3) &
               write(6,*) '*** a2 .le. 0 in z-sweep at ',i
             write(6,*) 'stopped with a2 < 0...'
             stop
             endif
         a(i) = dsqrt(a2)
         g1a2(i) = gamma1 / a2
         euv(i) = enth(i) - u2v2w2(i)
   10 continue
!
!
!     # now split the jump in q1d at each interface into waves
!
!     # find a1 thru a5, the coefficients of the 5 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
         delta(1) = ql(1,i) - qr(1,i-1)
         delta(2) = ql(mu,i) - qr(mu,i-1)
         delta(3) = ql(mv,i) - qr(mv,i-1)
         delta(4) = ql(mw,i) - qr(mw,i-1)
         delta(5) = ql(5,i) - qr(5,i-1)
         a4 = g1a2(i) * (euv(i)*delta(1) &
           + u(i)*delta(2) + v(i)*delta(3) + w(i)*delta(4) &
           - delta(5))
         a2 = delta(3) - v(i)*delta(1)
         a3 = delta(4) - w(i)*delta(1)
         a5 = (delta(2) + (a(i)-u(i))*delta(1) - a(i)*a4) / (2.d0*a(i))
         a1 = delta(1) - a4 - a5
!
!        # Compute the waves.
!        # Note that the 2-wave, 3-wave and 4-wave travel at the same speed
!        # and are lumped together in wave(.,2,.).  The 5-wave is then stored
!        # in wave(.,3,.).
!
         wave(1,1,i)  = a1
         wave(mu,1,i) = a1*(u(i)-a(i))
         wave(mv,1,i) = a1*v(i)
         wave(mw,1,i) = a1*w(i)
         wave(5,1,i)  = a1*(enth(i) - u(i)*a(i))
         s(1,i) = u(i)-a(i)
!
         wave(1,2,i)  = a4
         wave(mu,2,i) = a4*u(i)
         wave(mv,2,i) = a4*v(i) + a2
         wave(mw,2,i) = a4*w(i) + a3
         wave(5,2,i)  = a4*0.5d0*u2v2w2(i)  + a2*v(i) + a3*w(i)
         s(2,i) = u(i)
!
         wave(1,3,i)  = a5
         wave(mu,3,i) = a5*(u(i)+a(i))
         wave(mv,3,i) = a5*v(i)
         wave(mw,3,i) = a5*w(i)
         wave(5,3,i)  = a5*(enth(i)+u(i)*a(i))
         s(3,i) = u(i)+a(i)
   20    continue
!
!
!    # compute flux differences amdq and apdq.
!    ---------------------------------------
!
      if (efix) go to 110
!
!     # no entropy fix
!     ----------------
!
!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves
!
      do 100 m=1,meqn
         do 100 i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mws=1,mwaves
                if (s(m,iws) .lt. 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mws,i)*wave(m,mws,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mws,i)*wave(m,mws,i)
                endif
   90          continue
  100       continue
      go to 900
!
!-----------------------------------------------------
!
  110 continue
!
!     # With entropy fix
!     ------------------
!
!    # compute flux differences amdq and apdq.
!    # First compute amdq as sum of s*wave for left going waves.
!    # Incorporate entropy fix by adding a modified fraction of wave
!    # if s should change sign.
!
      do 200 i = 2-mbc, mx+mbc
!
!        # check 1-wave:
!        ---------------
!
         rhoim1 = qr(1,i-1)
         pim1 = gamma1*(qr(5,i-1) - 0.5d0*(qr(mu,i-1)**2 &
                + qr(mv,i-1)**2 + qr(mw,i-1)**2) / rhoim1)
         cim1 = dsqrt(gamma*pim1/rhoim1)
         s0 = qr(mu,i-1)/rhoim1 - cim1     !# u-c in left state (cell i-1)
!
!
!        # check for fully supersonic case:
         if (s0.ge.0.d0 .and. s(1,i).gt.0.d0)then
!            # everything is right-going
              do 60 m=1,meqn
                amdq(m,i) = 0.d0
   60           continue
              go to 200
              endif
!
         rho1 = qr(1,i-1) + wave(1,1,i)
         rhou1 = qr(mu,i-1) + wave(mu,1,i)
         rhov1 = qr(mv,i-1) + wave(mv,1,i)
         rhow1 = qr(mw,i-1) + wave(mw,1,i)
         en1 = qr(5,i-1) + wave(5,1,i)
         p1 = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2 + &
                     rhow1**2)/rho1)
         c1 = dsqrt(gamma*p1/rho1)
         s1 = rhou1/rho1 - c1  !# u-c to right of 1-wave
         if (s0.lt.0.d0 .and. s1.gt.0.d0) then
!            # transonic rarefaction in the 1-wave
           sfract = s0 * (s1-s(1,i)) / (s1-s0)
         else if (s(1,i) .lt. 0.d0) then
!           # 1-wave is leftgoing
           sfract = s(1,i)
         else
!           # 1-wave is rightgoing
             sfract = 0.d0   !# this shouldn't happen since s0 < 0
         endif
       do 120 m=1,meqn
          amdq(m,i) = sfract*wave(m,1,i)
  120       continue
!
!        # check 2-wave:
!        ---------------
!
         if (s(2,i) .ge. 0.d0) go to 200  !# 2-,3- and 4- waves are rightgoing
       do 140 m=1,meqn
          amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
  140       continue
!
!        # check 3-wave:
!        ---------------
!
       rhoi = ql(1,i)
       pi = gamma1*(ql(5,i) - 0.5d0*(ql(mu,i)**2 &
                + ql(mv,i)**2 + ql(mw,i)**2) / rhoi)
       ci = dsqrt(gamma*pi/rhoi)
       s3 = ql(mu,i)/rhoi + ci     !# u+c in right state  (cell i)
!
         rho2 = ql(1,i) - wave(1,3,i)
         rhou2 = ql(mu,i) - wave(mu,3,i)
         rhov2 = ql(mv,i) - wave(mv,3,i)
         rhow2 = ql(mw,i) - wave(mw,3,i)
         en2 = ql(5,i) - wave(5,3,i)
         p2 = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2 + &
                     rhow2**2)/rho2)
         c2 = dsqrt(gamma*p2/rho2)
         s2 = rhou2/rho2 + c2   !# u+c to left of 3-wave
         if (s2 .lt. 0.d0 .and. s3.gt.0.d0 ) then
!            # transonic rarefaction in the 3-wave
           sfract = s2 * (s3-s(3,i)) / (s3-s2)
         else if (s(3,i) .lt. 0.d0) then
!            # 3-wave is leftgoing
           sfract = s(3,i)
         else
!            # 3-wave is rightgoing
           go to 200
         endif
!
       do 160 m=1,5
          amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
  160       continue
  200    continue
!
!     # compute the rightgoing flux differences:
!     # df = SUM s*wave   is the total flux difference and apdq = df - amdq
!
      do 220 m=1,meqn
       do 220 i = 2-mbc, mx+mbc
          df = 0.d0
          do 210 mws=1,mwaves
             df = df + s(mws,i)*wave(m,mws,i)
  210          continue
          apdq(m,i) = df - amdq(m,i)
  220       continue
!
  900 continue
      return
      end

