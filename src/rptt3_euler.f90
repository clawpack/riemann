!     ==================================================================
      subroutine rptt3  (ixyz,icoor,maxm,meqn,mwaves,mbc,mx, &
                       ql,qr,aux1,aux2,aux3,maux,ilr,impt,asdq, &
                       bmasdq,bpasdq)
!     ==================================================================
!
!     # Riemann solver in the transverse direction for the
!     # Euler equations.
!     #
!     # On input,
!
!     #    ql,qr is the data along some one-dimensional slice, as in rpn3
!     #         This slice is
!     #             in the x-direction if ixyz=1,
!     #             in the y-direction if ixyz=2, or
!     #             in the z-direction if ixyz=3.
!
!     #    bsasdq is an array of flux differences that result from a
!     #         transverse splitting (a previous call to rpt3).  
!     #         This stands for B^* A^* \Dq but could represent any of 
!     #         6 possibilities, e.g.  C^* B^* \Dq, as specified by ixyz
!     #         and icoor (see below).
!     #         Moreover, each * represents either + or -, as specified by
!     #         imp and impt.
!
!     #    ixyz indicates the direction of the original Riemann solve,
!     #         called the x-like direction in the table below:
!
!     #               x-like direction   y-like direction   z-like direction
!     #      ixyz=1:        x                  y                  z         
!     #      ixyz=2:        y                  z                  x         
!     #      ixyz=3:        z                  x                  y         
!
!     #    icoor indicates direction in which the transverse solve should 
!     #         be performed.
!     #      icoor=2: split in the y-like direction.
!     #      icoor=3: split in the z-like direction.
!
!     #    For example,
!     #        ixyz=1, icoor=3 means bsasdq=B^*A^*\Dq, and should be 
!     #                        split in z into 
!     #                           cmbsasdq = C^-B^*A^*\Dq,
!     #                           cpbsasdq = C^+B^*A^*\Dq.
!     #
!     #        ixyz=2, icoor=3 means bsasdq=C^*B^*\Dq, and should be
!     #                        split in x into 
!     #                           cmbsasdq = A^-C^*B^*\Dq,
!     #                           cpbsasdq = A^+C^*B^*\Dq.
!
!     #    The parameters imp and impt are generally needed only if aux
!     #    arrays are being used, in order to access the appropriate
!     #    variable coefficients.
!
      implicit real*8(a-h,o-z)
      dimension     ql(1-mbc:maxm+mbc, meqn)
      dimension     qr(1-mbc:maxm+mbc, meqn)
      dimension   asdq(1-mbc:maxm+mbc, meqn)
      dimension bmasdq(1-mbc:maxm+mbc, meqn)
      dimension bpasdq(1-mbc:maxm+mbc, meqn)
      dimension   aux1(1-mbc:maxm+mbc, maux, 3)
      dimension   aux2(1-mbc:maxm+mbc, maux, 3)
      dimension   aux3(1-mbc:maxm+mbc, maux, 3)
!
      dimension waveb(5,3),sb(3)
      parameter (maxmrp = 1002)
      integer tflag

      double precision u2v2w2(-1:maxmrp), &
            u(-1:maxmrp),v(-1:maxmrp),w(-1:maxmrp),enth(-1:maxmrp), &
            a(-1:maxmrp),g1a2(-1:maxmrp),euv(-1:maxmrp)
!
      if (-3.gt.1-mbc .or. maxmrp .lt. maxm+mbc) then
      write(6,*) 'need to increase maxmrp in rp3t'
      stop
      endif
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
!     # Solve Riemann problem in the second coordinate direction
!
      if( icoor .eq. 2 )then

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

      do 20 i = 2-mbc, mx+mbc
            a4 = g1a2(i) * (euv(i)*asdq(1,i) &
                  + u(i)*asdq(mu,i) + v(i)*asdq(mv,i) &
                  + w(i)*asdq(mw,i) - asdq(5,i))
         a2 = asdq(mu,i) - u(i)*asdq(1,i)
            a3 = asdq(mw,i) - w(i)*asdq(1,i)
         a5 = (asdq(mv,i) + (a(i)-v(i))*asdq(1,i) - a(i)*a4) &
                   / (2.d0*a(i))
         a1 = asdq(1,i) - a4 - a5
!
            waveb(1,1)  = a1
            waveb(mu,1) = a1*u(i)
            waveb(mv,1) = a1*(v(i)-a(i))
            waveb(mw,1) = a1*w(i)
            waveb(5,1)  = a1*(enth(i) - v(i)*a(i))
         sb(1) = v(i) - a(i)
!
            waveb(1,2)  = a4
            waveb(mu,2) = a2 + u(i)*a4
            waveb(mv,2) = v(i)*a4
            waveb(mw,2) = a3 + w(i)*a4
            waveb(5,2)  = a4*0.5d0*u2v2w2(i) + a2*u(i) + a3*w(i)
         sb(2) = v(i)
!
            waveb(1,3)  = a5
            waveb(mu,3) = a5*u(i)
            waveb(mv,3) = a5*(v(i)+a(i))
            waveb(mw,3) = a5*w(i)
            waveb(5,3)  = a5*(enth(i)+v(i)*a(i))
         sb(3) = v(i) + a(i)
!
      do 25 m=1,meqn
         bmasdq(m,i) = 0.d0
         bpasdq(m,i) = 0.d0
         do 25 mws=1,mwaves
            bmasdq(m,i) = bmasdq(m,i) &
                       + dmin1(sb(mws), 0.d0) * waveb(m,mws)
            bpasdq(m,i) = bpasdq(m,i) &
                       + dmax1(sb(mws), 0.d0) * waveb(m,mws)
 25         continue
!
   20    continue
!
      else
!
!        # Solve Riemann problem in the third coordinate direction
!
      do 30 i = 2-mbc, mx+mbc
            a4 = g1a2(i) * (euv(i)*asdq(1,i) &
                  + u(i)*asdq(mu,i) + v(i)*asdq(mv,i) &
                  + w(i)*asdq(mw,i) - asdq(5,i))
         a2 = asdq(mu,i) - u(i)*asdq(1,i)
            a3 = asdq(mv,i) - v(i)*asdq(1,i)
         a5 = (asdq(mw,i) + (a(i)-w(i))*asdq(1,i) - a(i)*a4) &
                   / (2.d0*a(i))
         a1 = asdq(1,i) - a4 - a5
!
            waveb(1,1)  = a1
            waveb(mu,1) = a1*u(i)
            waveb(mv,1) = a1*v(i)
            waveb(mw,1) = a1*(w(i) - a(i))
            waveb(5,1)  = a1*(enth(i) - w(i)*a(i))
         sb(1) = w(i) - a(i)
!
            waveb(1,2)  = a4
            waveb(mu,2) = a2 + u(i)*a4
            waveb(mv,2) = a3 + v(i)*a4
            waveb(mw,2) = w(i)*a4
            waveb(5,2)  = a4*0.5d0*u2v2w2(i) + a2*u(i) + a3*v(i)
         sb(2) = w(i)
!
            waveb(1,3)  = a5
            waveb(mu,3) = a5*u(i)
            waveb(mv,3) = a5*v(i)
            waveb(mw,3) = a5*(w(i)+a(i))
            waveb(5,3)  = a5*(enth(i)+w(i)*a(i))
         sb(3) = w(i) + a(i)
!
      do 35 m=1,meqn
         bmasdq(m,i) = 0.d0
         bpasdq(m,i) = 0.d0
         do 35 mws=1,mwaves
            bmasdq(m,i) = bmasdq(m,i) &
                       + dmin1(sb(mws), 0.d0) * waveb(m,mws)
            bpasdq(m,i) = bpasdq(m,i) &
                       + dmax1(sb(mws), 0.d0) * waveb(m,mws)
 35         continue
!
 30      continue
!
      endif
!
      return
      end

