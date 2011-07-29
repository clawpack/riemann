c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &			imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision(a-h,o-z)
c
c     # Riemann solver in the transverse direction for the advection equation.
c
      dimension     ql(meqn, 1-mbc:maxm+mbc)
      dimension     qr(meqn, 1-mbc:maxm+mbc)
      dimension   asdq(meqn, 1-mbc:maxm+mbc)
      dimension bmasdq(meqn, 1-mbc:maxm+mbc)
      dimension bpasdq(meqn, 1-mbc:maxm+mbc)
      dimension   aux1(2,    1-mbc:maxm+mbc)
      dimension   aux2(2,    1-mbc:maxm+mbc)
      dimension   aux3(2,    1-mbc:maxm+mbc)
c
c
      kv = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1
      do 10 i=2-mbc,mx+mbc
	 i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
	 bmasdq(i,1) = dmin1(aux2(i1,kv), 0.d0) * asdq(i,1)
	 bpasdq(i,1) = dmax1(aux3(i1,kv), 0.d0) * asdq(i,1)
   10    continue
c
      return
      end
