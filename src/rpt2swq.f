c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &          imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the shallow water
c     # equations  on a quadrilateral grid.
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # Uses Roe averages and other quantities which were
c     # computed in rpn2sh and stored in the common block comroe.
c
      dimension     ql(meqn, 1-mbc:maxm+mbc)
      dimension     qr(meqn, 1-mbc:maxm+mbc)
      dimension   asdq(meqn, 1-mbc:maxm+mbc)
      dimension bmasdq(meqn, 1-mbc:maxm+mbc)
      dimension bpasdq(meqn, 1-mbc:maxm+mbc)
      dimension   aux1(*, 1-mbc:maxm+mbc)
      dimension   aux2(*, 1-mbc:maxm+mbc)
      dimension   aux3(*, 1-mbc:maxm+mbc)
c
      parameter (maxm2 = 1002)  !# assumes at most 1000x1000 grid with mbc=2
      dimension u(-1:maxm2),v(-1:maxm2),a(-1:maxm2),h(-1:maxm2)
      dimension wave(4, 3, -1:maxm2)
      dimension    s(3, -1:maxm2)
      dimension enx(-1:maxm2), eny(-1:maxm2), enz(-1:maxm2)
      dimension etx(-1:maxm2), ety(-1:maxm2), etz(-1:maxm2)
      dimension gamma(-1:maxm2)
     
      dimension delta(4)
      common /sw/  g
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

c
      if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
	 write(6,*) 'need to increase maxm2 in rpB'
	 stop
      endif

      if(ixy.eq.1) then
        dx = dxcom
      else
        dx = dycom
      endif 
c
c
      if(ixy.eq.1) then
        ioff = 7
      else
        ioff = 1
      endif

c
c        # imp is used to flag whether wave is going to left or right,
c        # since states and grid are different on each side
c
         if (imp.eq.1) then
c            # asdq = amdq, moving to left
             ix1 = 2-mbc
	     ixm1 = mx+mbc
           else
c            # asdq = apdq, moving to right
             ix1 = 1-mbc
	     ixm1 = mx+mbc
           endif
c
c        --------------
c        # up-going:
c        --------------
c

c       # determine rotation matrix for interface above cell, using aux3
c               [ alf  beta ]
c               [-beta  alf ]
c
        do i=ix1,ixm1
c
         if (imp.eq.1) then
             i1 = i-1
           else
             i1 = i
           endif
c
           enx(i) =   aux3(ioff+1,i1)
           eny(i) =   aux3(ioff+2,i1)
           enz(i) =   aux3(ioff+3,i1)
           etx(i) =   aux3(ioff+4,i1)
           ety(i) =   aux3(ioff+5,i1)
           etz(i) =   aux3(ioff+6,i1)
           gamma(i) = dsqrt(etx(i)**2 + ety(i)**2 + etz(i)**2)
           etx(i) =   etx(i) / gamma(i)
           ety(i) =   ety(i) / gamma(i)
           etz(i) =   etz(i) / gamma(i)

           h(i) = ql(1,i1)
           u(i) = (enx(i)*ql(2,i1)+eny(i)*ql(3,i1)+enz(i)*ql(4,i1)) 
     &                / h(i)
           v(i) = (etx(i)*ql(2,i1)+ety(i)*ql(3,i1)+etz(i)*ql(4,i1)) 
     &                / h(i)
	   a(i) = dsqrt(g*h(i))
           enddo
c
c
c
c     # now split asdq into waves:
c
c     # find a1 thru a3, the coefficients of the 3 eigenvectors:
      do 20 i = ix1,ixm1
         delta(1) = asdq(1,i) 
         delta(2) = enx(i)*asdq(2,i)+eny(i)*asdq(3,i)+enz(i)*asdq(4,i)
         delta(3) = etx(i)*asdq(2,i)+ety(i)*asdq(3,i)+etz(i)*asdq(4,i)
         a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
         a2 = -v(i)*delta(1) + delta(3)
         a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))
c
c        # Compute the waves.
c
         wave(1,1,i) = a1
         wave(2,1,i) = enx(i)*a1*(u(i)-a(i)) + etx(i)*a1*v(i)
         wave(3,1,i) = eny(i)*a1*(u(i)-a(i)) + ety(i)*a1*v(i)
         wave(4,1,i) = enz(i)*a1*(u(i)-a(i)) + etz(i)*a1*v(i) 
         s(1,i) = (u(i)-a(i))*gamma(i)/dx
c
         wave(1,2,i) = 0.0d0
         wave(2,2,i) = etx(i)*a2
         wave(3,2,i) = ety(i)*a2
         wave(4,2,i) = etz(i)*a2
         s(2,i) = u(i) * gamma(i)/dx
c
         wave(1,3,i) = a3
         wave(2,3,i) = enx(i)*a3*(u(i)+a(i)) + etx(i)*a3*v(i)
         wave(3,3,i) = eny(i)*a3*(u(i)+a(i)) + ety(i)*a3*v(i)
         wave(4,3,i) = enz(i)*a3*(u(i)+a(i)) + etz(i)*a3*v(i)
         s(3,i) = (u(i)+a(i)) * gamma(i)/dx
   20    continue
c
c
c    # compute flux difference bpasdq
c    --------------------------------
c
      do 40 m=1,meqn
         do 40 i=ix1,ixm1
	    bpasdq(m,i) = 0.d0
	    do 30 mw=1,mwaves
	       bpasdq(m,i) = bpasdq(m,i) 
     &	           + dmax1(s(mw,i),0.d0)*wave(m,mw,i)
   30          continue
   40       continue
c
c     # project momentum component of bpasdq to tangent plane:
      do i=ix1, ixm1
        if (imp.eq.1) then
          i1 = i-1
        else
          i1 = i
        endif
        erx = aux3(14,i1)
        ery = aux3(15,i1)
        erz = aux3(16,i1)
        bn = erx*bpasdq(2,i)+ery*bpasdq(3,i)+erz*bpasdq(4,i)
        bpasdq(2,i) = bpasdq(2,i) - bn*erx
        bpasdq(3,i) = bpasdq(3,i) - bn*ery
        bpasdq(4,i) = bpasdq(4,i) - bn*erz

      enddo
 
c
c
c        --------------
c        # down-going:
c        --------------
c

c       # determine rotation matrix for interface below cell, using aux2
c               [ alf  beta ]
c               [-beta  alf ]
c
        do i=ix1,ixm1
c
         if (imp.eq.1) then
             i1 = i-1
           else
             i1 = i
           endif
c
           enx(i) =   aux2(ioff+1,i1)
           eny(i) =   aux2(ioff+2,i1)
           enz(i) =   aux2(ioff+3,i1)
           etx(i) =   aux2(ioff+4,i1)
           ety(i) =   aux2(ioff+5,i1)
           etz(i) =   aux2(ioff+6,i1)
           gamma(i) = dsqrt(etx(i)**2 + ety(i)**2 + etz(i)**2)
           etx(i) =   etx(i) / gamma(i)
           ety(i) =   ety(i) / gamma(i)
           etz(i) =   etz(i) / gamma(i)
           u(i) = (enx(i)*ql(2,i1)+eny(i)*ql(3,i1)+enz(i)*ql(4,i1)) 
     &               / h(i)
           v(i) = (etx(i)*ql(2,i1)+ety(i)*ql(3,i1)+etz(i)*ql(4,i1))
     &               / h(i)
           enddo
c
c
c
c     # now split asdq into waves:
c
c     # find a1 thru a3, the coefficients of the 3 eigenvectors:
      do 80 i = ix1,ixm1
         delta(1) = asdq(1,i) 
         delta(2) = enx(i)*asdq(2,i)+eny(i)*asdq(3,i)+enz(i)*asdq(4,i)
         delta(3) = etx(i)*asdq(2,i)+ety(i)*asdq(3,i)+etz(i)*asdq(4,i)

         a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
         a2 = -v(i)*delta(1) + delta(3)
         a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))
c
c        # Compute the waves.
c
         wave(1,1,i) = a1
         wave(2,1,i) = enx(i)*a1*(u(i)-a(i)) + etx(i)*a1*v(i)
         wave(3,1,i) = eny(i)*a1*(u(i)-a(i)) + ety(i)*a1*v(i)
         wave(4,1,i) = enz(i)*a1*(u(i)-a(i)) + etz(i)*a1*v(i)
         s(1,i) = (u(i)-a(i)) * gamma(i)/dx
c
         wave(1,2,i) = 0.0d0
         wave(2,2,i) = etx(i)*a2
         wave(3,2,i) = ety(i)*a2
         wave(4,2,i) = etz(i)*a2
         s(2,i) = u(i) * gamma(i)/dx
c
         wave(1,3,i) = a3
         wave(2,3,i) = enx(i)*a3*(u(i)+a(i)) + etx(i)*a3*v(i)
         wave(3,3,i) = eny(i)*a3*(u(i)+a(i)) + ety(i)*a3*v(i)
         wave(4,3,i) = enz(i)*a3*(u(i)+a(i)) + etz(i)*a3*v(i)
         s(3,i) = (u(i)+a(i)) * gamma(i)/dx
   80    continue
c
c
c    # compute flux difference bmasdq
c    --------------------------------
c
      do 100 m=1,meqn
         do 100 i=ix1,ixm1
	    bmasdq(m,i) = 0.d0
	    do 90 mw=1,mwaves
	       bmasdq(m,i) = bmasdq(m,i) 
     &	                    + dmin1(s(mw,i), 0.d0)*wave(m,mw,i)
   90          continue
  100       continue
c     # project momentum component of bmasdq to tangent plane:
      do i=ix1, ixm1
        if (imp.eq.1) then
          i1 = i-1
        else
          i1 = i
        endif
        erx = aux1(14,i1)
        ery = aux1(15,i1)
        erz = aux1(16,i1)
        bn = erx*bmasdq(2,i)+ery*bmasdq(3,i)+erz*bmasdq(4,i)
        bmasdq(2,i) = bmasdq(2,i) - bn*erx
        bmasdq(3,i) = bmasdq(3,i) - bn*ery
        bmasdq(4,i) = bmasdq(4,i) - bn*erz

      enddo
c
c
    
      return
      end
