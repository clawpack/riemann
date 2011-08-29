c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &                   aux1,aux2,aux3,imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the acoustics equations
c     # with varying material properties rho and kappa
c
c     # auxN(1,i) holds rho
c     # auxN(2,i) holds c
c     #  N = 1 for row below
c     #      2 for this row
c     #      3 for row above
c
c     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.
c
c     # imp=1  means  asdq=amdq,    imp=2 means asdq=apdq
c
      dimension    ql(meqn, 1-mbc:maxm+mbc)
      dimension    qr(meqn, 1-mbc:maxm+mbc)
      dimension    asdq(meqn, 1-mbc:maxm+mbc)
      dimension bmasdq(meqn, 1-mbc:maxm+mbc)
      dimension bpasdq(meqn, 1-mbc:maxm+mbc)
      dimension   aux1(2, 1-mbc:maxm+mbc)
      dimension   aux2(2, 1-mbc:maxm+mbc)
      dimension   aux3(2, 1-mbc:maxm+mbc)
c
c
c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
c
      do 20 i = 2-mbc, mx+mbc
c
c        # imp is used to flag whether wave is going to left or right,
c        # since material properties are different on the two sides
c
         if (imp.eq.1) then 
c            # asdq = amdq, moving to left
             i1 = i-1
           else
c            # asdq = apdq, moving to right
             i1 = i
           endif
c
c        # The flux difference asdq is split into downward moving part
c        # traveling at speed -c relative to the medium below and
c        # an upward moving part traveling
c        # at speed +c relative to the medium above.
c
c        # Note that the sum of these parts does not give all of asdq
c        # since there is also reflection at the interfaces which decreases
c        # the flux.
c
c        # sound speed in each row of cells:
         cm = aux1(2,i1)
         c = aux2(2,i1)
         cp = aux3(2,i1)
c
c        # impedances:
         zm = aux1(1,i1)*aux1(2,i1)
         zz = aux2(1,i1)*aux2(2,i1)
         zp = aux3(1,i1)*aux3(2,i1)

c        # transmitted part of down-going wave:
         a1 = (-asdq(1,i) + asdq(mv,i)*zz) /
     &         (zm + zz)

c        # transmitted part of up-going wave:
         a2 = (asdq(1,i) + asdq(mv,i)*zz) /
     &         (zz + zp)
c
c        # The down-going flux difference bmasdq is the product  -c * wave
c
         bmasdq(1,i) = cm * a1*zm
         bmasdq(mu,i) = 0.d0
         bmasdq(mv,i) = -cm * a1
c
c        # The up-going flux difference bpasdq is the product  c * wave
c
         bpasdq(1,i) = cp * a2*zp
         bpasdq(mu,i) = 0.d0
         bpasdq(mv,i) = cp * a2
c
   20    continue
c
      return
      end
